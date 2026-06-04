import argparse
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


DEFAULT_PARQUET = (
    "/gpfs/helios/projects/eQTLCatalogue/coverage_parquet/INTERVAL/blood.parquet"
)


def format_size(num_bytes):
    gib = num_bytes / (1024**3)
    gb = num_bytes / 1_000_000_000
    return f"{gib:.3f} GiB ({gb:.3f} GB)"


def print_basic_metadata(path, parquet_file):
    metadata = parquet_file.metadata

    print(f"File: {path}")
    print(f"Size: {format_size(path.stat().st_size)}")
    print(f"Rows: {metadata.num_rows}")
    print(f"Columns: {metadata.num_columns}")
    print(f"Row groups: {parquet_file.num_row_groups}")
    print("\nSchema:")
    print(parquet_file.schema)


def check_metadata_null_counts(parquet_file):
    null_counts = {}
    missing_null_stats = 0

    for row_group_index in range(parquet_file.num_row_groups):
        row_group = parquet_file.metadata.row_group(row_group_index)

        for column_index in range(row_group.num_columns):
            column = row_group.column(column_index)
            stats = column.statistics

            if stats is None or stats.null_count is None:
                missing_null_stats += 1
                continue

            if stats.null_count:
                name = column.path_in_schema
                null_counts[name] = null_counts.get(name, 0) + stats.null_count

    print("\nMetadata Null Check:")
    if null_counts:
        print("Columns with nulls recorded in parquet metadata:")
        for name, count in sorted(null_counts.items()):
            print(f"  {name}: {count}")
    else:
        print("No nulls recorded in parquet metadata.")

    if missing_null_stats:
        print(
            "Warning: "
            f"{missing_null_stats} column chunks did not include null-count statistics."
        )


def count_batch_nans(array):
    if pa.types.is_floating(array.type):
        return pc.sum(pc.is_nan(array)).as_py() or 0
    return 0


def scan_batches(parquet_file, batch_size):
    null_counts = {}
    nan_counts = {}
    rows_scanned = 0
    batches_scanned = 0

    print(f"\nBatch Scan: batch_size={batch_size}")

    for batch in parquet_file.iter_batches(batch_size=batch_size):
        batches_scanned += 1
        rows_scanned += batch.num_rows

        for name, array in zip(batch.schema.names, batch.columns):
            null_count = array.null_count
            if null_count:
                null_counts[name] = null_counts.get(name, 0) + null_count

            nan_count = count_batch_nans(array)
            if nan_count:
                nan_counts[name] = nan_counts.get(name, 0) + nan_count

        if batches_scanned % 100 == 0:
            print(
                f"  scanned {batches_scanned} batches / {rows_scanned} rows...",
                flush=True,
            )

    print(f"Scanned rows: {rows_scanned}")
    print(f"Scanned batches: {batches_scanned}")

    if null_counts:
        print("\nColumns with nulls from data scan:")
        for name, count in sorted(null_counts.items()):
            print(f"  {name}: {count}")
    else:
        print("\nNo nulls found in data scan.")

    if nan_counts:
        print("\nFloating-point columns with NaN values:")
        for name, count in sorted(nan_counts.items()):
            print(f"  {name}: {count}")
    else:
        print("No NaN values found in floating-point columns.")


def sanity_check(path, batch_size, metadata_only):
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    parquet_file = pq.ParquetFile(path)
    print_basic_metadata(path, parquet_file)
    check_metadata_null_counts(parquet_file)

    if metadata_only:
        print("\nSkipped data scan because --metadata-only was set.")
        return

    scan_batches(parquet_file, batch_size)


def main():
    parser = argparse.ArgumentParser(
        description="Low-memory sanity check for a generated parquet file."
    )
    parser.add_argument(
        "file_path",
        nargs="?",
        default=DEFAULT_PARQUET,
        help=f"Parquet file to check. Defaults to {DEFAULT_PARQUET}",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10_000,
        help="Rows per Arrow batch for the data scan. Default: 10000",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Only inspect parquet metadata and skip scanning the data.",
    )
    args = parser.parse_args()

    sanity_check(Path(args.file_path), args.batch_size, args.metadata_only)


if __name__ == "__main__":
    main()
