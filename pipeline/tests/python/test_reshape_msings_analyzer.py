"""tests for bin/reshape_msings_analyzer.py"""

import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
spec = importlib.util.spec_from_file_location("reshape_msings_analyzer", BIN / "reshape_msings_analyzer.py")
mod = importlib.util.module_from_spec(spec)
sys.modules["reshape_msings_analyzer"] = mod
spec.loader.exec_module(mod)


def test_reshape_skips_headers_and_comments(tmp_path: Path):
    src = tmp_path / "analyzer.txt"
    dst = tmp_path / "out.tsv"
    src.write_text(
        "# msings header comment\n"
        "chrom\tstart\tend\tcounts\n"               # header
        "1:10234-10245\textra\t10:42;11:3\n"         # locus row
        "chr1:20000-20010\tfoo\t9:10;10:5\n"         # chr-prefixed
        "\n"                                          # blank
        "not_a_locus\t1\t2\n"                         # skipped: no ':' + '-'
    )
    n = mod.reshape(str(src), str(dst))
    assert n == 2
    lines = dst.read_text().splitlines()
    assert lines == ["1:10234-10245\t10:42;11:3", "chr1:20000-20010\t9:10;10:5"]


def test_reshape_empty_input(tmp_path: Path):
    src = tmp_path / "a.txt"
    dst = tmp_path / "o.tsv"
    src.write_text("# only comments\n\n")
    n = mod.reshape(str(src), str(dst))
    assert n == 0
    assert dst.read_text() == ""
