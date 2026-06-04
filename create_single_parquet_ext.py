#!/usr/bin/env python3
"""
create_single_parquet_ext.py

Converts a directory of BigWig files into a single Parquet file using a
memory-efficient, two-stage (map-reduce) approach. This is designed to handle
very large datasets on systems with limited memory, such as HPC nodes.

The final Parquet file is sorted by chromosome and start position to enable
efficient querying of genomic regions.

Stage 1 (Map):
- Each BigWig file is processed in parallel.
- The data for each file is transformed and written to a temporary Parquet file
  on disk, releasing memory.

Stage 2 (Reduce):
- DuckDB reads all the intermediate Parquet files.
- The data is sorted and aggregated efficiently, chromosome by chromosome.
- The final single Parquet file is written to the destination.

Usage:
    python create_single_parquet_ext.py <input_dir> <output_path> [--num_processes N] [--chrom CHROM] [--pos_start N] [--pos_end N]

Example:
    python create_single_parquet_ext.py data_qc_passed coverage_data.parquet --num_processes 8 --chrom 22 --pos_start 16000000 --pos_end 51000000
"""
import sys
import argparse
import os
import glob
from multiprocessing import Pool, cpu_count
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyBigWig
from tqdm import tqdm
import shutil
import duckdb

# Define the target schema to ensure consistent data types and column order
TARGET_SCHEMA = pa.schema([
    pa.field('seqnames', pa.string()),
    pa.field('start', pa.int32()),
    pa.field('end', pa.int32()),
    pa.field('strand', pa.string()),
    pa.field('score', pa.float64()),
    pa.field('sample_id', pa.string())
])


def sql_string_literal(value):
    """Return a DuckDB SQL string literal for paths and simple filter values."""
    return "'" + str(value).replace("'", "''") + "'"


def stream_parquet_to_writer(parquet_path, writer, batch_size):
    """
    Append a Parquet file to an existing writer in bounded-size Arrow batches.
    """
    parquet_file = pq.ParquetFile(parquet_path)
    rows_written = 0

    for batch in parquet_file.iter_batches(
        batch_size=batch_size,
        columns=TARGET_SCHEMA.names,
    ):
        if batch.num_rows == 0:
            continue

        table = pa.Table.from_batches([batch])
        if table.schema != TARGET_SCHEMA:
            table = table.cast(TARGET_SCHEMA)

        writer.write_table(table)
        rows_written += table.num_rows

    return rows_written


def process_and_save_bigwig(task):
    """
    Wrapper function for parallel processing. Reads a BigWig file, converts it
    to a DataFrame with the target schema, and saves it to a temporary Parquet file.

    Args:
        task (tuple): A tuple containing task parameters.
    """
    filepath, output_path, allowed_chroms, pos_start, pos_end = task

    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)

    bw = None
    try:
        bw = pyBigWig.open(filepath)
        if bw is None:
            raise RuntimeError("pyBigWig.open returned None")

        data_frames = []
        for chrom_name, _ in bw.chroms().items():
            clean_chrom = chrom_name.replace('chr', '')
            if clean_chrom in allowed_chroms:
                intervals = bw.intervals(chrom_name)
                if intervals:
                    df = pd.DataFrame(intervals, columns=['start', 'end', 'score'])
                    # BigWig intervals are 0-based, half-open: [start, end).
                    # Convert to 1-based, inclusive coordinates (GRanges-style): [start+1, end].
                    df['start'] = df['start'] + 1
                    df = df[df['end'] >= df['start']]

                    if pos_start is not None:
                        df = df[df['end'] >= pos_start]
                    if pos_end is not None:
                        df = df[df['start'] <= pos_end]

                    if df.empty:
                        continue

                    df['seqnames'] = clean_chrom
                    data_frames.append(df)

        rows_written = 0
        if data_frames:
            full_df = pd.concat(data_frames, ignore_index=True)
            sample_id = os.path.basename(filepath).split('.')[0]
            full_df['sample_id'] = sample_id
            full_df['strand'] = '*'

            # Reorder columns and ensure schema is correct before saving
            full_df = full_df[['seqnames', 'start', 'end', 'strand', 'score', 'sample_id']]
            rows_written = len(full_df)

            table = pa.Table.from_pandas(full_df, schema=TARGET_SCHEMA, preserve_index=False)
            pq.write_table(table, output_path)

        return {
            'filepath': filepath,
            'ok': True,
            'rows_written': rows_written,
        }
    except Exception as e:
        return {
            'filepath': filepath,
            'ok': False,
            'error': str(e),
        }
    finally:
        if bw is not None:
            bw.close()


def extract_bigwigs(input_dir, output_dir, allowed_chroms, pos_start, pos_end, num_processes):
    """
    Convert BigWig files to intermediate Parquet parts.
    """
    bw_files = glob.glob(os.path.join(input_dir, '*.bw')) + \
               glob.glob(os.path.join(input_dir, '*.bigwig'))

    if not bw_files:
        print(f"Error: No .bw or .bigwig files found in '{input_dir}'", file=sys.stderr)
        sys.exit(1)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print(f"Found {len(bw_files)} BigWig files. Processing in parallel with {num_processes} processes...")
    print(f"Intermediate files will be stored in '{output_dir}'")

    tasks = [
        (bw_file, os.path.join(output_dir, f"part_{i}.parquet"), allowed_chroms, pos_start, pos_end)
        for i, bw_file in enumerate(bw_files)
    ]

    with Pool(num_processes) as pool:
        results = list(tqdm(
            pool.imap_unordered(process_and_save_bigwig, tasks),
            total=len(tasks),
            desc="[Stage 1/2] Processing BigWig files",
        ))

    failures = [result for result in results if not result['ok']]
    if failures:
        print(
            f"\nError: {len(failures)} BigWig file(s) could not be opened or processed.",
            file=sys.stderr,
        )
        for failure in sorted(failures, key=lambda item: item['filepath']):
            print(f"  {failure['filepath']}: {failure['error']}", file=sys.stderr)
        print("\nAborting without writing final Parquet output.", file=sys.stderr)
        sys.exit(1)

    no_data = [result for result in results if result['rows_written'] == 0]
    if no_data:
        print(
            f"Warning: {len(no_data)} BigWig file(s) had no intervals after filtering.",
            file=sys.stderr,
        )

    return output_dir


def make_unique_temp_dir(base_temp_dir):
    temp_dir = base_temp_dir
    counter = 1
    while os.path.exists(temp_dir):
        temp_dir = f"{base_temp_dir}_{counter}"
        counter += 1
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir


def merge_parquet_parts(parts_dir, output_path, selected_chrom, args):
    """
    Merge intermediate Parquet parts into the final sorted Parquet file.
    """
    part_files = glob.glob(os.path.join(parts_dir, '*.parquet'))
    if not part_files:
        print(f"Warning: no intermediate Parquet parts found in '{parts_dir}'. Writing empty output.")
        with pq.ParquetWriter(output_path, TARGET_SCHEMA):
            pass
        return 0

    all_parts_path = os.path.join(parts_dir, '*.parquet')
    chromosomes = [selected_chrom] if selected_chrom is not None else [str(i) for i in range(1, 23)]
    merge_temp_dir = make_unique_temp_dir("merge_parquet_tmp")
    chrom_temp_dir = os.path.join(merge_temp_dir, 'chromosome_parquet_parts')
    duckdb_temp_dir = args.duckdb_temp_dir or os.path.join(merge_temp_dir, 'duckdb_tmp')
    duckdb_threads = args.duckdb_threads or min(args.num_processes, 4)
    os.makedirs(chrom_temp_dir, exist_ok=True)
    os.makedirs(duckdb_temp_dir, exist_ok=True)

    con = None
    rows_written = 0
    try:
        con = duckdb.connect()
        con.execute(f"SET temp_directory = {sql_string_literal(duckdb_temp_dir)}")
        con.execute(f"SET threads = {duckdb_threads}")
        con.execute("SET preserve_insertion_order = false")
        if args.duckdb_memory_limit:
            con.execute(f"SET memory_limit = {sql_string_literal(args.duckdb_memory_limit)}")

        print("[Stage 2/2] Writing final dataset by chromosome...")
        print(f"DuckDB temp directory: {duckdb_temp_dir}")
        print(f"DuckDB memory limit: {args.duckdb_memory_limit}")
        print(f"DuckDB threads: {duckdb_threads}")

        with pq.ParquetWriter(output_path, TARGET_SCHEMA) as writer:
            for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
                chrom_output_path = os.path.join(chrom_temp_dir, f"chrom_{chrom}.parquet")
                where_clauses = [f"seqnames = {sql_string_literal(chrom)}"]
                if args.pos_start is not None:
                    where_clauses.append(f"\"end\" >= {args.pos_start}")
                if args.pos_end is not None:
                    where_clauses.append(f"start <= {args.pos_end}")

                copy_query = f"""
                    COPY (
                        SELECT
                            CAST(seqnames AS VARCHAR) AS seqnames,
                            CAST(start AS INTEGER) AS start,
                            CAST("end" AS INTEGER) AS "end",
                            CAST(strand AS VARCHAR) AS strand,
                            CAST(score AS DOUBLE) AS score,
                            CAST(sample_id AS VARCHAR) AS sample_id
                        FROM read_parquet({sql_string_literal(all_parts_path)})
                        WHERE {' AND '.join(where_clauses)}
                        ORDER BY start
                    )
                    TO {sql_string_literal(chrom_output_path)}
                    (FORMAT PARQUET, COMPRESSION 'SNAPPY')
                """
                con.execute(copy_query)

                if os.path.exists(chrom_output_path):
                    rows_written += stream_parquet_to_writer(
                        chrom_output_path,
                        writer,
                        args.stream_batch_size,
                    )
                    os.remove(chrom_output_path)
    finally:
        if con is not None:
            con.close()
        if os.path.exists(merge_temp_dir):
            shutil.rmtree(merge_temp_dir)

    return rows_written


def main():
    """
    Main function to orchestrate the conversion of BigWig files to a
    single Parquet file using a memory-efficient strategy.
    """
    parser = argparse.ArgumentParser(
        description="Convert a directory of BigWig files to a single, sorted Parquet file."
    )
    parser.add_argument(
        'input_dir',
        help='Input directory containing BigWig (.bw or .bigwig) files.'
    )
    parser.add_argument(
        'output_path',
        help='Output path for the single Parquet file, or parts directory in extract mode.'
    )
    parser.add_argument(
        '--mode',
        choices=['full', 'extract', 'merge'],
        default='full',
        help='Run both stages, only BigWig extraction, or only final merge (default: full).'
    )
    parser.add_argument(
        '--num_processes', '-p',
        type=int,
        default=cpu_count(),
        help='Number of parallel processes to use (default: number of CPU cores).'
    )
    parser.add_argument(
        '--chrom',
        type=str,
        default=None,
        help="Optional chromosome to keep (e.g. 22, X, chr22). Default: 1-22 (autosomes)."
    )
    parser.add_argument(
        '--pos_start',
        type=int,
        default=None,
        help='Optional 1-based inclusive start position filter (requires --chrom).'
    )
    parser.add_argument(
        '--pos_end',
        type=int,
        default=None,
        help='Optional 1-based inclusive end position filter (requires --chrom).'
    )
    parser.add_argument(
        '--duckdb_memory_limit',
        default='64GB',
        help='DuckDB memory limit for stage 2 sorting before spilling to disk (default: 64GB).'
    )
    parser.add_argument(
        '--duckdb_threads',
        type=int,
        default=None,
        help='DuckDB thread count for stage 2 (default: min(num_processes, 4)).'
    )
    parser.add_argument(
        '--duckdb_temp_dir',
        default=None,
        help='Directory for DuckDB temporary spill files (default: inside the stage 1 temp directory).'
    )
    parser.add_argument(
        '--stream_batch_size',
        type=int,
        default=65536,
        help='Rows per Arrow batch when streaming chromosome Parquet files into the final output.'
    )
    args = parser.parse_args()

    valid_chroms = {str(i) for i in range(1, 23)} | {'X'}
    autosomes = {str(i) for i in range(1, 23)}
    selected_chrom = None
    if args.chrom is not None:
        selected_chrom = args.chrom.replace('chr', '')
        if selected_chrom not in valid_chroms:
            parser.error("--chrom must be one of 1-22 or X (with or without 'chr').")

    if (args.pos_start is not None or args.pos_end is not None) and selected_chrom is None:
        parser.error("--pos_start and --pos_end require --chrom.")

    if args.pos_start is not None and args.pos_start < 1:
        parser.error("--pos_start must be >= 1.")
    if args.pos_end is not None and args.pos_end < 1:
        parser.error("--pos_end must be >= 1.")
    if args.pos_start is not None and args.pos_end is not None and args.pos_end < args.pos_start:
        parser.error("--pos_end must be >= --pos_start.")
    if args.duckdb_threads is not None and args.duckdb_threads < 1:
        parser.error("--duckdb_threads must be >= 1.")
    if args.stream_batch_size < 1:
        parser.error("--stream_batch_size must be >= 1.")

    allowed_chroms = {selected_chrom} if selected_chrom is not None else autosomes

    if selected_chrom is not None:
        print(f"Chromosome filter: {selected_chrom}")
    if args.pos_start is not None or args.pos_end is not None:
        print(f"Position filter: start={args.pos_start}, end={args.pos_end}")

    if args.mode == 'extract':
        extract_bigwigs(
            args.input_dir,
            args.output_path,
            allowed_chroms,
            args.pos_start,
            args.pos_end,
            args.num_processes,
        )
        print("\nDone.")
        print(f"Successfully created intermediate Parquet parts at: {args.output_path}")
        return

    if args.mode == 'merge':
        print("Merging intermediate Parquet parts into final output...")
        merge_parquet_parts(args.input_dir, args.output_path, selected_chrom, args)
        print("\nDone.")
        print(f"Successfully created single Parquet file at: {args.output_path}")
        return

    temp_dir = make_unique_temp_dir("temp_parquet_parts")
    try:
        extract_bigwigs(
            args.input_dir,
            temp_dir,
            allowed_chroms,
            args.pos_start,
            args.pos_end,
            args.num_processes,
        )
        print("\nAll BigWig files processed. Merging and writing final Parquet file...")
        merge_parquet_parts(temp_dir, args.output_path, selected_chrom, args)
    finally:
        if os.path.exists(temp_dir):
            print(f"Cleaning up temporary directory '{temp_dir}'...")
            shutil.rmtree(temp_dir)

    print("\nDone.")
    print(f"Successfully created single Parquet file at: {args.output_path}")


if __name__ == '__main__':
    main()
