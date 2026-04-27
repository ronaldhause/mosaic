#!/usr/bin/env python3
"""parse misa output into a padded, autosome-only bed.

misa (http://pgrc.ipk-gatersleben.de/misa/) emits a tab-delimited file
with columns: id, ssr_nr, ssr_type, ssr, size, start, end. with
`-interrupts=0` every row is a single perfect repeat; compound/complex
merging happens downstream in merge_compound_complex.py.

input coords are 1-based inclusive; bed is 0-based half-open. we also
pad +/- pad_bp and clamp to chrom bounds (read from the fai).

emits a tsv with columns:
    chrom, start, end, locus_id, repeat_subunit, subunit_length, n_repeats
where locus_id = "chrom:start-end" using the padded coordinates with
chrom stripped of any leading "chr" (matches the paper's convention).
"""

from __future__ import annotations

import argparse
import re
import sys

import pandas as pd

SSR_RE = re.compile(r"^\(([ACGTN]+)\)(\d+)$", re.IGNORECASE)
AUTOSOMES = {str(i) for i in range(1, 23)} | {f"chr{i}" for i in range(1, 23)}


def read_chrom_sizes(fai_path: str) -> dict[str, int]:
    """load contig lengths from a samtools .fai index"""
    sizes: dict[str, int] = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes


def parse_ssr(ssr: str) -> tuple[str, int] | None:
    """split a misa ssr column like '(AT)6' into (subunit, n_repeats)"""
    m = SSR_RE.match(ssr.strip())
    if not m:
        return None
    return m.group(1).upper(), int(m.group(2))


def misa_to_bed(
    misa_path: str,
    chrom_sizes: dict[str, int],
    pad: int = 5,
) -> pd.DataFrame:
    """parse misa output into a padded bed dataframe.

    keeps autosomes only; clamps padded coords to chrom bounds. drops any
    row whose ssr doesn't match the simple (unit)count form (i.e. anything
    misa produced with interrupts > 0 that slipped through).
    """
    df = pd.read_csv(
        misa_path,
        sep="\t",
        comment="#",
        header=0,
        dtype=str,
    )
    df.columns = [c.strip().lower().replace(".", "").replace(" ", "_") for c in df.columns]
    # misa's default header is: "ID  SSR nr.  SSR type  SSR  size  start  end"
    # which normalises to: id, ssr_nr, ssr_type, ssr, size, start, end
    needed = {"id", "ssr", "start", "end"}
    if not needed.issubset(df.columns):
        raise ValueError(
            f"missa output missing columns {needed - set(df.columns)}; "
            f"saw: {list(df.columns)}"
        )

    df = df[df["id"].isin(AUTOSOMES)].copy()
    parsed = df["ssr"].map(parse_ssr)
    keep = parsed.notna()
    df = df[keep].copy()
    df["repeat_subunit"] = [p[0] for p in parsed[keep]]
    df["n_repeats"] = [p[1] for p in parsed[keep]]
    df["subunit_length"] = df["repeat_subunit"].str.len()

    # coord conversion: misa 1-based inclusive -> bed 0-based half-open,
    # then pad.
    df["start"] = df["start"].astype(int) - 1 - pad
    df["end"] = df["end"].astype(int) + pad
    df["start"] = df["start"].clip(lower=0)

    # clamp end to contig length (strip leading chr for lookup if needed)
    def _cap(row):
        chrom = row["id"]
        size = chrom_sizes.get(chrom) or chrom_sizes.get(f"chr{chrom}") or chrom_sizes.get(chrom.replace("chr", ""))
        return min(row["end"], size) if size else row["end"]

    df["end"] = df.apply(_cap, axis=1)
    df = df[df["end"] > df["start"]]

    # normalise chrom naming for the locus_id (paper convention: no "chr")
    df["chrom"] = df["id"]
    df["id_bare"] = df["chrom"].str.replace(r"^chr", "", regex=True)
    df["locus_id"] = df["id_bare"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)

    out = df[
        ["chrom", "start", "end", "locus_id", "repeat_subunit", "subunit_length", "n_repeats"]
    ].sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--misa", required=True, help="misa output file (*.misa)")
    p.add_argument("--fai", required=True, help="reference fasta index (.fai)")
    p.add_argument("--pad", type=int, default=5, help="pad +/- bp on each side [%(default)s]")
    p.add_argument("--output", required=True, help="output tsv")
    args = p.parse_args()

    sizes = read_chrom_sizes(args.fai)
    out = misa_to_bed(args.misa, sizes, pad=args.pad)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(out)} autosomal loci to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
