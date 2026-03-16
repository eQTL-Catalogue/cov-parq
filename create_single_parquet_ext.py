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

    try:
        bw = pyBigWig.open(filepath)
    except Exception as e:
        print(f"Error opening {filepath}: {e}", file=sys.stderr)
        return

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

    if data_frames:
        full_df = pd.concat(data_frames, ignore_index=True)
        sample_id = os.path.basename(filepath).split('.')[0]
        full_df['sample_id'] = sample_id
        full_df['strand'] = '*'

        # Reorder columns and ensure schema is correct before saving
        full_df = full_df[['seqnames', 'start', 'end', 'strand', 'score', 'sample_id']]

        table = pa.Table.from_pandas(full_df, schema=TARGET_SCHEMA, preserve_index=False)
        pq.write_table(table, output_path)

    bw.close()


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
        help='Output path for the single Parquet file.'
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

    allowed_chroms = {selected_chrom} if selected_chrom is not None else autosomes

    bw_files = glob.glob(os.path.join(args.input_dir, '*.bw')) + \
               glob.glob(os.path.join(args.input_dir, '*.bigwig'))

    if not bw_files:
        print(f"Error: No .bw or .bigwig files found in '{args.input_dir}'", file=sys.stderr)
        sys.exit(1)

    base_temp_dir = "temp_parquet_parts"
    temp_dir = base_temp_dir
    counter = 1
    while os.path.exists(temp_dir):
        temp_dir = f"{base_temp_dir}_{counter}"
        counter += 1

    print(f"Found {len(bw_files)} BigWig files. Processing in parallel with {args.num_processes} processes...")
    print(f"Intermediate files will be stored in '{temp_dir}'")
    if selected_chrom is not None:
        print(f"Chromosome filter: {selected_chrom}")
    if args.pos_start is not None or args.pos_end is not None:
        print(f"Position filter: start={args.pos_start}, end={args.pos_end}")

    tasks = [
        (bw_file, os.path.join(temp_dir, f"part_{i}.parquet"), allowed_chroms, args.pos_start, args.pos_end)
        for i, bw_file in enumerate(bw_files)
    ]

    with Pool(args.num_processes) as pool:
        list(tqdm(pool.imap_unordered(process_and_save_bigwig, tasks), total=len(tasks), desc="[Stage 1/2] Processing BigWig files"))

    print("\nAll BigWig files processed. Merging and writing final Parquet file...")

    all_parts_path = os.path.join(temp_dir, '*.parquet')
    chromosomes = [selected_chrom] if selected_chrom is not None else [str(i) for i in range(1, 23)]

    try:
        con = duckdb.connect()
        print("[Stage 2/2] Writing final dataset by chromosome...")

        with pq.ParquetWriter(args.output_path, TARGET_SCHEMA) as writer:
            for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
                where_clauses = [f"seqnames = '{chrom}'"]
                if args.pos_start is not None:
                    where_clauses.append(f"\"end\" >= {args.pos_start}")
                if args.pos_end is not None:
                    where_clauses.append(f"start <= {args.pos_end}")

                query = f"""
                    SELECT seqnames, start, "end", strand, score, sample_id
                    FROM read_parquet('{all_parts_path}')
                    WHERE {' AND '.join(where_clauses)}
                    ORDER BY start
                """
                table = con.execute(query).fetch_arrow_table()

                if table.num_rows > 0:
                    writer.write_table(table)
    finally:
        con.close()
        print(f"Cleaning up temporary directory '{temp_dir}'...")
        shutil.rmtree(temp_dir)

    print("\nDone.")
    print(f"Successfully created single Parquet file at: {args.output_path}")


if __name__ == '__main__':
    main()
