#!/usr/bin/env python3
"""write build_summary.tsv with counts for every stage of build_reference.nf."""

from __future__ import annotations

import argparse
import sys

import pandas as pd


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--loci_raw", required=True, help="misa_to_bed.py output tsv")
    p.add_argument("--merged_tsv", required=True, help="merge_compound_complex.py tsv output")
    p.add_argument("--ms_bed", required=True, help="microsatellites.bed after bait filter")
    p.add_argument("--defb_locus", default="8:7679723-7679741",
                   help="defb locus id to probe for in the final bed")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    raw = pd.read_csv(args.loci_raw, sep="\t")
    merged = pd.read_csv(args.merged_tsv, sep="\t")
    kept = pd.read_csv(
        args.ms_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "locus_id", "score", "strand"],
    )
    autosome_names = {str(i) for i in range(1, 23)}
    raw_chroms = raw["chrom"].astype(str).str.replace("^chr", "", regex=True)
    defb_present = bool(
        (kept["locus_id"] == args.defb_locus).any()
        or (kept["locus_id"].str.replace("^chr", "", regex=True) == args.defb_locus).any()
    )

    rows = [
        ("n_loci_total_raw", len(raw)),
        ("n_loci_autosomal_raw", int(raw_chroms.isin(autosome_names).sum())),
        ("n_loci_merged", len(merged)),
        ("n_loci_after_bait_filter", len(kept)),
        ("n_compound", int(merged["is_compound"].sum()) if "is_compound" in merged else 0),
        ("n_complex", int(merged["is_complex"].sum()) if "is_complex" in merged else 0),
        ("defb_locus_present", defb_present),
    ]
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(args.output, sep="\t", index=False)
    for k, v in rows:
        print(f"{k}\t{v}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
