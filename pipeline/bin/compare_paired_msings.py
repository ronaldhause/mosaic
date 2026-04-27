#!/usr/bin/env python3
"""compute peak_diff per locus from paired tumor/normal msings outputs.

implements the paper's "high-sensitivity" approach:

1. for each locus, per sample, compute relative abundance of each tract
   length as count / max(count at that locus).
2. drop lengths with relative abundance < min_rel_abundance (default 5%).
3. count remaining distinct lengths (peaks).
4. peak_diff = peaks_tumor - peaks_normal. positive = unstable under this
   definition.

input format (tall tsv, one row per (locus, length), tab-delimited):
    locus<TAB>length<TAB>count
where locus is e.g. "8:7679723-7679741", length is the offset in bp
(int), count is the supporting read count (int).

convert msings analyzer output to this form upstream. the wide/legacy
msings format ("offset:count" pairs in one column) is also accepted
when the file has two columns: locus<TAB>offset_count_string (with
pairs separated by ';' or ',').

output: 3-column csv: sample, locus, peak_diff
(directly consumed by pipeline/bin/compute_features.R)
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd


def _parse_wide_line(locus: str, blob: str) -> list[tuple[str, int, int]]:
    rows = []
    for item in blob.replace(",", ";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) < 2:
            continue
        try:
            length = int(parts[0])
            count = int(float(parts[1]))
        except ValueError:
            continue
        rows.append((locus, length, count))
    return rows


def read_msings_tall(path: str) -> pd.DataFrame:
    """accept either canonical 3-col (locus, length, count) or 2-col wide"""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#")
    if df.shape[1] >= 3:
        df = df.iloc[:, :3].copy()
        df.columns = ["locus", "length", "count"]
        df["length"] = df["length"].astype(int)
        df["count"] = df["count"].astype(float).astype(int)
        return df
    if df.shape[1] == 2:
        rows: list[tuple[str, int, int]] = []
        for locus, blob in df.itertuples(index=False):
            rows.extend(_parse_wide_line(str(locus), str(blob)))
        return pd.DataFrame(rows, columns=["locus", "length", "count"])
    raise ValueError(f"{path}: expected 2 or 3 tab-delimited columns")


def count_peaks(dist: pd.DataFrame, min_rel_abundance: float) -> pd.Series:
    """return one row per locus: the number of lengths above the relabund cutoff"""
    if dist.empty:
        return pd.Series(dtype=int, name="peaks")
    peak_max = dist.groupby("locus")["count"].transform("max")
    keep = dist["count"] >= peak_max * min_rel_abundance
    kept = dist[keep & (dist["count"] > 0)]
    return kept.groupby("locus")["length"].nunique().rename("peaks")


def call_instability(
    tumor: pd.DataFrame,
    normal: pd.DataFrame,
    min_rel_abundance: float = 0.05,
) -> pd.DataFrame:
    """compute per-locus peak_diff between tumor and matched normal.

    drops tract lengths with relative abundance < min_rel_abundance at each
    locus (normalized to the most frequent length at that locus). peak_diff
    is (tumor peaks) - (normal peaks) across the union of loci; loci absent
    in either sample contribute a peak count of 0.

    returns columns: locus, peak_diff.
    """
    t_peaks = count_peaks(tumor, min_rel_abundance)
    n_peaks = count_peaks(normal, min_rel_abundance)
    all_loci = t_peaks.index.union(n_peaks.index)
    out = pd.DataFrame(index=all_loci)
    out["tumor_peaks"] = t_peaks.reindex(all_loci).fillna(0).astype(int)
    out["normal_peaks"] = n_peaks.reindex(all_loci).fillna(0).astype(int)
    out["peak_diff"] = out["tumor_peaks"] - out["normal_peaks"]
    out.index.name = "locus"
    return out.reset_index()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tumor", required=True, help="tumor per-locus length distribution tsv")
    p.add_argument("--normal", required=True, help="normal per-locus length distribution tsv")
    p.add_argument("--sample_name", required=True, help="sample id for the output")
    p.add_argument("--min_rel_abundance", type=float, default=0.05,
                   help="drop lengths below this relative abundance [%(default)s]")
    p.add_argument("--output", required=True,
                   help="output tsv: sample, locus, peak_diff (consumed by compute_features.R)")
    args = p.parse_args()

    tumor = read_msings_tall(args.tumor)
    normal = read_msings_tall(args.normal)
    calls = call_instability(tumor, normal, min_rel_abundance=args.min_rel_abundance)
    calls.insert(0, "sample", args.sample_name)
    calls[["sample", "locus", "peak_diff"]].to_csv(
        args.output, sep="\t", index=False
    )
    n_unstable = int((calls["peak_diff"] > 0).sum())
    print(
        f"{args.sample_name}: {len(calls)} loci, {n_unstable} unstable "
        f"(peak_diff > 0)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
