#!/usr/bin/env python3
"""validate a mosaic samplesheet csv and echo it (with defaults filled) to stdout."""

import argparse
import csv
import sys
from pathlib import Path

REQUIRED = ["sample_name", "tumor_bam", "normal_bam"]
OPTIONAL = ["tumor_type"]


def die(msg: str) -> None:
    """print error to stderr and exit non-zero."""
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(1)


def has_bai(bam: Path) -> bool:
    """check both <bam>.bai and <bam>.bai (replacing .bam) conventions."""
    return bam.with_suffix(bam.suffix + ".bai").exists() or bam.with_suffix(".bai").exists()


def validate(path: Path) -> list[dict]:
    """read csv, validate required fields and bam/bai existence, return rows."""
    with path.open() as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            die(f"empty samplesheet: {path}")
        missing = [c for c in REQUIRED if c not in reader.fieldnames]
        if missing:
            die(f"missing required columns: {', '.join(missing)}")
        rows = list(reader)

    if not rows:
        die("samplesheet has no data rows")

    seen: set[str] = set()
    for i, row in enumerate(rows, start=2):  # line numbers include header
        name = (row.get("sample_name") or "").strip()
        if not name:
            die(f"line {i}: empty sample_name")
        if name in seen:
            die(f"line {i}: duplicate sample_name '{name}'")
        seen.add(name)

        for col in ("tumor_bam", "normal_bam"):
            bam_path = Path((row.get(col) or "").strip())
            if not bam_path:
                die(f"line {i}: empty {col}")
            if not bam_path.exists():
                die(f"line {i}: BAM not found: {bam_path} ({col})")
            if not has_bai(bam_path):
                die(f"line {i}: BAI index not found for {bam_path}")

        if "tumor_type" not in row or not (row.get("tumor_type") or "").strip():
            row["tumor_type"] = "UNKNOWN"

    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="validate a mosaic samplesheet")
    ap.add_argument("samplesheet", type=Path)
    args = ap.parse_args()

    if not args.samplesheet.exists():
        die(f"samplesheet not found: {args.samplesheet}")

    rows = validate(args.samplesheet)
    cols = REQUIRED + OPTIONAL
    writer = csv.DictWriter(sys.stdout, fieldnames=cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
