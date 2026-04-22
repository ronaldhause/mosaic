#!/usr/bin/env python3
"""join merged-loci metadata with bait-filtered loci (and optional annovar).

produces the final locus_annotation.tsv used by downstream analyses —
one row per locus retained after the bait-proximity filter.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

COLS = [
    "locus_id", "chrom", "start", "end", "repeat_type",
    "repeat_subunits", "subunit_lengths", "n_repeats",
    "is_compound", "is_complex", "n_members", "members",
    "genomic_class", "gene",
]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--merged_tsv", required=True, help="merge_compound_complex.py tsv output")
    p.add_argument("--ms_bed", required=True, help="microsatellites.bed (post bait filter)")
    p.add_argument("--annovar_tsv", default=None,
                   help="optional annovar tsv with chrom,start,end,genomic_class,gene")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    merged = pd.read_csv(args.merged_tsv, sep="\t")
    kept = pd.read_csv(
        args.ms_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "locus_id", "score", "strand"],
    )
    ann = merged.merge(kept[["locus_id"]], on="locus_id", how="inner")

    if args.annovar_tsv:
        annovar = pd.read_csv(args.annovar_tsv, sep="\t")
        ann = ann.merge(
            annovar[["chrom", "start", "end", "genomic_class", "gene"]],
            on=["chrom", "start", "end"], how="left",
        )
    else:
        ann["genomic_class"] = pd.NA
        ann["gene"] = pd.NA

    ann[COLS].to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(ann)} annotated loci to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
