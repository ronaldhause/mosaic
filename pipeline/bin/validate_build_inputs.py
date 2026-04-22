#!/usr/bin/env python3
"""sanity-check inputs for build_reference.nf before the pipeline fires.

checks:
- reference fasta exists and is non-empty
- matching .fai exists (samtools faidx) and lists at least one autosome
- capture bait bed exists, is non-empty, and its chrom naming style
  (chr-prefixed vs not) matches the fasta

fails fast with a readable message. exits 0 on success.
"""

from __future__ import annotations

import argparse
import os
import sys


def _nonempty(path: str) -> None:
    if not os.path.isfile(path):
        raise SystemExit(f"file not found: {path}")
    if os.path.getsize(path) == 0:
        raise SystemExit(f"file is empty: {path}")


def _chrom_style(names: list[str]) -> str:
    """return 'chr' if majority have a chr prefix, 'bare' otherwise"""
    if not names:
        return "bare"
    n_chr = sum(1 for n in names if n.startswith("chr"))
    return "chr" if n_chr > len(names) / 2 else "bare"


def _fai_contigs(fai: str) -> list[str]:
    names: list[str] = []
    with open(fai) as fh:
        for line in fh:
            parts = line.split("\t", 1)
            if parts:
                names.append(parts[0])
    return names


def _bed_contigs(bed: str, sample_n: int = 200) -> list[str]:
    names: list[str] = []
    with open(bed) as fh:
        for i, line in enumerate(fh):
            if i >= sample_n:
                break
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            names.append(line.split("\t", 1)[0])
    return names


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True, help="hg19 fasta")
    p.add_argument("--fai", required=True, help=".fai index")
    p.add_argument("--capture_bait_bed", required=True, help="capture bait bed")
    args = p.parse_args()

    _nonempty(args.reference)
    _nonempty(args.fai)
    _nonempty(args.capture_bait_bed)

    fai_names = _fai_contigs(args.fai)
    if not fai_names:
        raise SystemExit(f"no contigs in {args.fai}")
    autosomes_in_fai = {n for n in fai_names if n.lstrip("chr").isdigit() and 1 <= int(n.lstrip("chr")) <= 22}
    if not autosomes_in_fai:
        raise SystemExit(f"no autosomal contigs (1-22) in {args.fai}; is this hg19?")

    bed_names = _bed_contigs(args.capture_bait_bed)
    if not bed_names:
        raise SystemExit(f"no data rows in {args.capture_bait_bed}")

    fai_style = _chrom_style(fai_names)
    bed_style = _chrom_style(bed_names)
    if fai_style != bed_style:
        raise SystemExit(
            f"chrom naming mismatch: fasta is '{fai_style}'-style "
            f"(e.g. {fai_names[0]}), bait bed is '{bed_style}'-style "
            f"(e.g. {bed_names[0]}). rename one to match before running."
        )

    print(
        f"ok: {len(autosomes_in_fai)} autosomes in fasta, "
        f"{len(bed_names)} sampled bait rows; naming = {fai_style}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
