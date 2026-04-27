#!/usr/bin/env python3
"""merge adjacent microsatellite loci into compound/complex entries.

paper rule: if two perfect repeats are separated by <= max_gap bp, treat
them as a single locus. if the constituent subunit lengths are all equal
it's "compound" (c); otherwise "complex" (c*).

input: tsv from misa_to_bed.py with columns
    chrom, start, end, locus_id, repeat_subunit, subunit_length, n_repeats

output: tsv with columns
    chrom, start, end, locus_id, repeat_type, repeat_subunits, subunit_lengths,
    n_repeats, is_compound, is_complex, n_members, members
and a matching 6-column bed (chrom, start, end, locus_id, score, strand).
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd


def merge_adjacent(loci: pd.DataFrame, max_gap: int = 10) -> pd.DataFrame:
    """merge adjacent microsatellite loci separated by <= max_gap bp.

    expects columns: chrom, start, end, locus_id, repeat_subunit,
    subunit_length, n_repeats.

    returns one row per merged locus with repeat_type in {p1..p5, c, c*},
    flags for is_compound and is_complex, and an audit 'members' column
    listing the original locus_ids that merged.
    """
    if loci.empty:
        return loci.assign(
            repeat_type=pd.Series(dtype=str),
            repeat_subunits=pd.Series(dtype=str),
            subunit_lengths=pd.Series(dtype=str),
            is_compound=pd.Series(dtype=bool),
            is_complex=pd.Series(dtype=bool),
            n_members=pd.Series(dtype=int),
            members=pd.Series(dtype=str),
        )

    df = loci.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    out_rows: list[dict] = []
    cur: dict | None = None

    def flush(row: dict) -> None:
        subs = row["_subunits"]
        lens = row["_lengths"]
        n_members = len(subs)
        is_compound = n_members > 1 and len(set(lens)) == 1
        is_complex = n_members > 1 and len(set(lens)) > 1
        if n_members == 1:
            repeat_type = f"p{lens[0]}"
        elif is_compound:
            repeat_type = "c"
        else:
            repeat_type = "c*"
        bare = row["chrom"].replace("chr", "", 1) if row["chrom"].startswith("chr") else row["chrom"]
        out_rows.append({
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "locus_id": f"{bare}:{row['start']}-{row['end']}",
            "repeat_type": repeat_type,
            "repeat_subunits": ",".join(subs),
            "subunit_lengths": ",".join(str(x) for x in lens),
            "n_repeats": ",".join(str(n) for n in row["_n_repeats"]),
            "is_compound": is_compound,
            "is_complex": is_complex,
            "n_members": n_members,
            "members": ",".join(row["_members"]),
        })

    for _, r in df.iterrows():
        if cur is None or r["chrom"] != cur["chrom"] or r["start"] - cur["end"] > max_gap:
            if cur is not None:
                flush(cur)
            cur = {
                "chrom": r["chrom"],
                "start": int(r["start"]),
                "end": int(r["end"]),
                "_subunits": [r["repeat_subunit"]],
                "_lengths": [int(r["subunit_length"])],
                "_n_repeats": [int(r["n_repeats"])],
                "_members": [r["locus_id"]],
            }
        else:
            cur["end"] = max(cur["end"], int(r["end"]))
            cur["_subunits"].append(r["repeat_subunit"])
            cur["_lengths"].append(int(r["subunit_length"]))
            cur["_n_repeats"].append(int(r["n_repeats"]))
            cur["_members"].append(r["locus_id"])

    if cur is not None:
        flush(cur)

    return pd.DataFrame(out_rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="misa_to_bed.py output tsv")
    p.add_argument("--max_gap", type=int, default=10,
                   help="max gap (bp) between adjacent loci to merge [%(default)s]")
    p.add_argument("--output_tsv", required=True, help="merged loci tsv")
    p.add_argument("--output_bed", required=True, help="merged loci bed (6-col)")
    args = p.parse_args()

    loci = pd.read_csv(args.input, sep="\t")
    merged = merge_adjacent(loci, max_gap=args.max_gap)
    merged.to_csv(args.output_tsv, sep="\t", index=False)

    bed = merged[["chrom", "start", "end", "locus_id"]].copy()
    bed["score"] = 0
    bed["strand"] = "+"
    bed.to_csv(args.output_bed, sep="\t", index=False, header=False)

    n_total = len(merged)
    n_compound = int(merged["is_compound"].sum()) if "is_compound" in merged else 0
    n_complex = int(merged["is_complex"].sum()) if "is_complex" in merged else 0
    print(
        f"merged {len(loci)} input loci into {n_total} (compound={n_compound}, complex={n_complex})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
