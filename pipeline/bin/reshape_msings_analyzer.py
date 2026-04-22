#!/usr/bin/env python3
"""reshape a msings 'msi analyzer' output into the wide 2-col form
consumed by compare_paired_msings.py.

msings versions differ on column names/ordering. heuristic: data rows
start with a coordinate string like '1:10234-10245' or 'chr1:10234-10245';
the last tab-delimited column carries the 'offset:count[:...]' pairs.

outputs: locus<TAB>offset_count_blob per locus.
"""

from __future__ import annotations

import argparse
import sys


def reshape(path_in: str, path_out: str) -> int:
    n = 0
    with open(path_in) as fh, open(path_out, "w") as out:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # header rows use words; data rows start with a locus coord
            if not (parts and ":" in parts[0] and "-" in parts[0]):
                continue
            locus = parts[0]
            blob = parts[-1]
            out.write(f"{locus}\t{blob}\n")
            n += 1
    return n


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    n = reshape(args.input, args.output)
    print(f"reshaped {n} loci from {args.input} -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
