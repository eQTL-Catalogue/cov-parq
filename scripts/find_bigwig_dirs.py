#!/usr/bin/env python3

import argparse
import csv
import sys
from pathlib import Path


HEADER = ("study_name", "qtl_group", "bigwig_dir")


def find_bigwig_dirs(root: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []

    for study_dir in sorted((p for p in root.iterdir() if p.is_dir()), key=lambda p: p.name):
        for qtl_group_dir in sorted(
            (p for p in study_dir.iterdir() if p.is_dir()), key=lambda p: p.name
        ):
            bigwig_dir = qtl_group_dir / "bigwig"
            if bigwig_dir.is_dir():
                rows.append((study_dir.name, qtl_group_dir.name, str(bigwig_dir)))

    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Scan an rnaseq root for exact {study_name}/{qtl_group}/bigwig directories "
            "and write a studies.tsv-style file."
        ),
        epilog=(
            "Example:\n"
            "  python scripts/find_bigwig_dirs.py "
            "--root /gpfs/helios/projects/eQTLCatalogue/r8_run_folders/rnaseq "
            "--out studies.auto.tsv"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Root directory to scan (e.g. .../r8_run_folders/rnaseq).",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV path (new file or overwrite existing file).",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    out_path = Path(args.out).expanduser().resolve()

    if not root.exists():
        print(f"Error: root path does not exist: {root}", file=sys.stderr)
        return 1
    if not root.is_dir():
        print(f"Error: root path is not a directory: {root}", file=sys.stderr)
        return 1

    rows = find_bigwig_dirs(root)
    if not rows:
        print(
            "Error: no matching directories found for pattern "
            "{root}/{study_name}/{qtl_group}/bigwig",
            file=sys.stderr,
        )
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADER)
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
