#!/usr/bin/env python3
"""
Open a BigWig file and print basic metadata plus the first few intervals.

Example:
  python scripts/check_bigwig.py /path/to/sample.bigwig
  python scripts/check_bigwig.py /path/to/sample.bigwig --chrom chr1 --n 20
  python scripts/check_bigwig.py /path/to/bigwig_dir
  python scripts/check_bigwig.py /path/to/bigwig_dir --recursive
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pyBigWig


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    sys.exit(1)


def format_interval(interval: tuple[int, int, float]) -> str:
    start, end, score = interval
    return f"{start}\t{end}\t{score}"


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    clean = chrom.removeprefix("chr")
    try:
        value = int(clean)
        return (value, chrom)
    except ValueError:
        if clean == "X":
            return (23, chrom)
        if clean == "Y":
            return (24, chrom)
        return (1000, chrom)


def find_bigwig_files(path: Path, recursive: bool) -> list[Path]:
    patterns = ("*.bw", "*.bigwig")
    files: list[Path] = []
    for pattern in patterns:
        files.extend(path.rglob(pattern) if recursive else path.glob(pattern))
    return sorted(files)


def check_bigwig_readable(bigwig_path: Path) -> tuple[bool, str]:
    try:
        bw = pyBigWig.open(str(bigwig_path))
    except Exception as exc:
        return False, str(exc)

    try:
        chroms = bw.chroms()
        if not chroms:
            return False, "opened, but no chromosomes were found"
        return True, f"{len(chroms)} chromosomes"
    except Exception as exc:
        return False, str(exc)
    finally:
        bw.close()


def report_directory(directory: Path, recursive: bool) -> None:
    if not directory.is_dir():
        fail(f"not a directory: {directory}")

    bigwig_files = find_bigwig_files(directory, recursive=recursive)
    if not bigwig_files:
        fail(f"no .bw or .bigwig files found in {directory}")

    bad_files: list[tuple[Path, str]] = []
    print(f"Checking {len(bigwig_files)} BigWig files in {directory}")

    for index, bigwig_path in enumerate(bigwig_files, start=1):
        ok, message = check_bigwig_readable(bigwig_path)
        if not ok:
            bad_files.append((bigwig_path, message))
            print(f"BAD\t{index}\t{bigwig_path.name}\t{message}")

    good_count = len(bigwig_files) - len(bad_files)
    print("\nSummary")
    print(f"  Total files: {len(bigwig_files)}")
    print(f"  Readable: {good_count}")
    print(f"  Unreadable: {len(bad_files)}")

    if bad_files:
        print("\nUnreadable files:")
        for bigwig_path, message in bad_files:
            print(f"  {bigwig_path}\t{message}")
        sys.exit(1)

    print("\nPASS: all BigWig files opened successfully.")


def print_bigwig_head(
    bigwig_path: Path,
    chrom: str | None,
    start: int,
    end: int | None,
    n: int,
) -> None:
    if not bigwig_path.exists():
        fail(f"file does not exist: {bigwig_path}")

    try:
        bw = pyBigWig.open(str(bigwig_path))
    except Exception as exc:
        fail(f"could not open BigWig: {exc}")

    try:
        chroms = bw.chroms()
        if not chroms:
            fail("BigWig opened, but no chromosomes were found")

        print(f"BigWig: {bigwig_path}")
        print(f"File size: {bigwig_path.stat().st_size} bytes")
        print(f"Chromosomes: {len(chroms)}")
        print("First chromosomes:")
        for chrom_name, chrom_size in sorted(chroms.items(), key=lambda item: chrom_sort_key(item[0]))[:10]:
            print(f"  {chrom_name}\t{chrom_size}")

        chroms_to_check = [chrom] if chrom is not None else [
            chrom_name for chrom_name, _ in sorted(chroms.items(), key=lambda item: chrom_sort_key(item[0]))
        ]

        for chrom_name in chroms_to_check:
            if chrom_name not in chroms:
                fail(f"chromosome '{chrom_name}' is not present in this BigWig")

            query_end = end if end is not None else chroms[chrom_name]
            query_end = min(query_end, chroms[chrom_name])
            if start >= query_end:
                fail(f"invalid range for {chrom_name}: start={start}, end={query_end}")

            intervals = bw.intervals(chrom_name, start, query_end)
            if not intervals:
                print(f"\nNo intervals found in {chrom_name}:{start}-{query_end}")
                continue

            print(f"\nHead intervals from {chrom_name}:{start}-{query_end}")
            print("start\tend\tscore")
            for interval in intervals[:n]:
                print(format_interval(interval))
            return

        fail("no intervals found in any chromosome")
    finally:
        bw.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Print the head of a BigWig file, or report unreadable BigWigs in a directory."
    )
    parser.add_argument("path", help="Path to a .bw/.bigwig file or a directory containing BigWigs.")
    parser.add_argument("--chrom", default=None, help="Optional chromosome to inspect, e.g. chr1 or 1.")
    parser.add_argument("--start", type=int, default=0, help="0-based query start position. Default: 0.")
    parser.add_argument("--end", type=int, default=None, help="Optional query end position.")
    parser.add_argument("--n", type=int, default=10, help="Number of intervals to print. Default: 10.")
    parser.add_argument("--recursive", action="store_true", help="Recursively scan directories for BigWigs.")
    args = parser.parse_args()

    if args.start < 0:
        parser.error("--start must be >= 0")
    if args.end is not None and args.end <= args.start:
        parser.error("--end must be greater than --start")
    if args.n < 1:
        parser.error("--n must be >= 1")

    return args


def main() -> None:
    args = parse_args()
    path = Path(args.path).expanduser().resolve()

    if path.is_dir():
        report_directory(path, recursive=args.recursive)
        return

    print_bigwig_head(
        bigwig_path=path,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        n=args.n,
    )


if __name__ == "__main__":
    main()
