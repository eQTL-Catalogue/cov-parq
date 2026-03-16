#!/usr/bin/env python3
"""
Validate cov-parq output parquet files for basic integrity and expected shape.

Example:
  python scripts/test_output.py /path/to/CIL_6.parquet
  python scripts/test_output.py /path/to/CIL_6_chr22_18500000_20000000.parquet \
      --expected-chrom 22 --expected-pos-start 18500000 --expected-pos-end 20000000
"""

from __future__ import annotations

import argparse
import os
import sys

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq


EXPECTED_COLUMNS = ["seqnames", "start", "end", "strand", "score", "sample_id"]
EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("seqnames", pa.string()),
        pa.field("start", pa.int32()),
        pa.field("end", pa.int32()),
        pa.field("strand", pa.string()),
        pa.field("score", pa.float64()),
        pa.field("sample_id", pa.string()),
    ]
)


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


def normalize_chrom(chrom: str | None) -> str | None:
    if chrom is None:
        return None
    return chrom.replace("chr", "")


def chrom_rank(seqname: str) -> int:
    try:
        value = int(seqname)
        if 1 <= value <= 22:
            return value
    except Exception:
        pass
    if seqname == "X":
        return 23
    return 1000


def validate_parquet_schema(parquet_path: str) -> None:
    try:
        pf = pq.ParquetFile(parquet_path)
    except Exception as exc:
        fail(f"cannot open parquet file: {exc}")

    schema = pf.schema_arrow
    columns = schema.names
    if columns != EXPECTED_COLUMNS:
        fail(f"unexpected columns/order: got {columns}, expected {EXPECTED_COLUMNS}")

    for name in EXPECTED_COLUMNS:
        got_type = schema.field(name).type
        expected_type = EXPECTED_SCHEMA.field(name).type
        if got_type != expected_type:
            fail(f"column '{name}' has type {got_type}, expected {expected_type}")

    print("OK: schema and column order match expected output.")


def run_integrity_checks(
    parquet_path: str,
    expected_chrom: str | None,
    expected_pos_start: int | None,
    expected_pos_end: int | None,
) -> None:
    con = duckdb.connect()
    try:
        parquet_path_escaped = parquet_path.replace("'", "''")
        rel_name = "parq"
        con.execute(f"CREATE VIEW {rel_name} AS SELECT * FROM read_parquet('{parquet_path_escaped}')")

        row_count = con.execute(f"SELECT COUNT(*) FROM {rel_name}").fetchone()[0]
        if row_count <= 0:
            fail("file has zero rows")
        print(f"OK: row count = {row_count}")

        null_count = con.execute(
            f"""
            SELECT COUNT(*)
            FROM {rel_name}
            WHERE seqnames IS NULL
               OR start IS NULL
               OR "end" IS NULL
               OR strand IS NULL
               OR score IS NULL
               OR sample_id IS NULL
            """
        ).fetchone()[0]
        if null_count > 0:
            fail(f"found {null_count} rows with NULLs in required columns")
        print("OK: no NULLs in required columns.")

        bad_interval_count = con.execute(
            f"""SELECT COUNT(*) FROM {rel_name} WHERE start < 1 OR "end" < start"""
        ).fetchone()[0]
        if bad_interval_count > 0:
            fail(f"found {bad_interval_count} invalid intervals (start < 1 or end < start)")
        print("OK: genomic intervals are valid.")

        strand_bad = con.execute(f"SELECT COUNT(*) FROM {rel_name} WHERE strand <> '*'").fetchone()[0]
        if strand_bad > 0:
            fail(f"found {strand_bad} rows where strand is not '*'")
        print("OK: strand values are as expected ('*').")

        chrom_counts = con.execute(
            f"""
            SELECT seqnames, COUNT(*) AS n
            FROM {rel_name}
            GROUP BY seqnames
            ORDER BY
                CASE
                    WHEN TRY_CAST(seqnames AS INTEGER) IS NOT NULL THEN TRY_CAST(seqnames AS INTEGER)
                    WHEN seqnames = 'X' THEN 23
                    ELSE 1000
                END,
                seqnames
            """
        ).fetchall()
        print(f"INFO: chromosomes present: {', '.join([f'{c}:{n}' for c, n in chrom_counts])}")

        if expected_chrom is not None:
            mismatch = con.execute(
                f"SELECT COUNT(*) FROM {rel_name} WHERE seqnames <> ?",
                [expected_chrom],
            ).fetchone()[0]
            if mismatch > 0:
                fail(f"found {mismatch} rows not on expected chromosome {expected_chrom}")
            print(f"OK: all rows are on expected chromosome {expected_chrom}.")

        if expected_pos_start is not None:
            mismatch = con.execute(
                f"""SELECT COUNT(*) FROM {rel_name} WHERE "end" < ?""",
                [expected_pos_start],
            ).fetchone()[0]
            if mismatch > 0:
                fail(f"found {mismatch} rows ending before expected_pos_start={expected_pos_start}")
            print(f"OK: all rows satisfy end >= {expected_pos_start}.")

        if expected_pos_end is not None:
            mismatch = con.execute(
                f"SELECT COUNT(*) FROM {rel_name} WHERE start > ?",
                [expected_pos_end],
            ).fetchone()[0]
            if mismatch > 0:
                fail(f"found {mismatch} rows starting after expected_pos_end={expected_pos_end}")
            print(f"OK: all rows satisfy start <= {expected_pos_end}.")

        sample_count = con.execute(f"SELECT COUNT(DISTINCT sample_id) FROM {rel_name}").fetchone()[0]
        print(f"INFO: distinct sample_id count = {sample_count}")

    finally:
        con.close()


def validate_file_row_order(parquet_path: str, batch_size: int = 250000) -> None:
    """
    Memory-safe streaming validation of physical row order in the parquet file.
    Ensures chromosome rank is nondecreasing and start is nondecreasing within each rank.
    """
    pf = pq.ParquetFile(parquet_path)
    prev_rank = None
    prev_start = None
    checked_rows = 0

    for batch in pf.iter_batches(columns=["seqnames", "start"], batch_size=batch_size):
        seqnames = batch.column(0).to_pylist()
        starts = batch.column(1).to_pylist()

        for seqname, start in zip(seqnames, starts):
            rank = chrom_rank(seqname)
            if prev_rank is not None:
                if rank < prev_rank:
                    fail(
                        f"ordering violation at row {checked_rows + 1}: "
                        f"chromosome rank decreased ({rank} < {prev_rank})"
                    )
                if rank == prev_rank and start < prev_start:
                    fail(
                        f"ordering violation at row {checked_rows + 1}: "
                        f"start decreased ({start} < {prev_start}) within chromosome {seqname}"
                    )
            prev_rank = rank
            prev_start = start
            checked_rows += 1

    print(f"OK: row order is nondecreasing by chromosome block then start ({checked_rows} rows checked).")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check output parquet integrity and expected content.")
    parser.add_argument("parquet_path", help="Path to parquet file to validate.")
    parser.add_argument(
        "--expected-chrom",
        default=None,
        help="Optional expected chromosome (e.g. 22 or chr22).",
    )
    parser.add_argument(
        "--expected-pos-start",
        type=int,
        default=None,
        help="Optional expected inclusive start filter bound.",
    )
    parser.add_argument(
        "--expected-pos-end",
        type=int,
        default=None,
        help="Optional expected inclusive end filter bound.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    parquet_path = os.path.abspath(args.parquet_path)
    if not os.path.exists(parquet_path):
        fail(f"file does not exist: {parquet_path}")

    expected_chrom = normalize_chrom(args.expected_chrom)
    if (args.expected_pos_start is not None or args.expected_pos_end is not None) and expected_chrom is None:
        fail("--expected-pos-start/--expected-pos-end require --expected-chrom")

    print(f"Checking parquet: {parquet_path}")
    validate_parquet_schema(parquet_path)
    run_integrity_checks(
        parquet_path=parquet_path,
        expected_chrom=expected_chrom,
        expected_pos_start=args.expected_pos_start,
        expected_pos_end=args.expected_pos_end,
    )
    validate_file_row_order(parquet_path)
    print("PASS: parquet file looks intact and consistent with expected output.")


if __name__ == "__main__":
    main()
