"""tests for bin/validate_samplesheet.py"""

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"
spec = importlib.util.spec_from_file_location("validate_samplesheet", BIN / "validate_samplesheet.py")
mod = importlib.util.module_from_spec(spec)
sys.modules["validate_samplesheet"] = mod
spec.loader.exec_module(mod)


def _write_csv(path: Path, rows: list[dict], cols: list[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def _touch_bam(dirpath: Path, name: str) -> Path:
    bam = dirpath / name
    bam.write_bytes(b"")
    bam.with_suffix(bam.suffix + ".bai").write_bytes(b"")
    return bam


def test_valid_samplesheet_fills_default_tumor_type(tmp_path: Path):
    t = _touch_bam(tmp_path, "t.bam")
    n = _touch_bam(tmp_path, "n.bam")
    sheet = tmp_path / "s.csv"
    _write_csv(sheet, [{"sample_name": "s1", "tumor_bam": str(t), "normal_bam": str(n)}],
               ["sample_name", "tumor_bam", "normal_bam"])
    rows = mod.validate(sheet)
    assert len(rows) == 1
    assert rows[0]["tumor_type"] == "UNKNOWN"


def test_missing_required_column_exits(tmp_path: Path):
    sheet = tmp_path / "s.csv"
    _write_csv(sheet, [{"sample_name": "s1"}], ["sample_name"])
    with pytest.raises(SystemExit):
        mod.validate(sheet)


def test_missing_bai_exits(tmp_path: Path):
    t = tmp_path / "t.bam"
    t.write_bytes(b"")  # no .bai
    n = _touch_bam(tmp_path, "n.bam")
    sheet = tmp_path / "s.csv"
    _write_csv(sheet, [{"sample_name": "s1", "tumor_bam": str(t), "normal_bam": str(n)}],
               ["sample_name", "tumor_bam", "normal_bam"])
    with pytest.raises(SystemExit):
        mod.validate(sheet)


def test_duplicate_sample_name_exits(tmp_path: Path):
    t = _touch_bam(tmp_path, "t.bam")
    n = _touch_bam(tmp_path, "n.bam")
    sheet = tmp_path / "s.csv"
    _write_csv(sheet, [
        {"sample_name": "s1", "tumor_bam": str(t), "normal_bam": str(n)},
        {"sample_name": "s1", "tumor_bam": str(t), "normal_bam": str(n)},
    ], ["sample_name", "tumor_bam", "normal_bam"])
    with pytest.raises(SystemExit):
        mod.validate(sheet)


def test_accepts_bai_without_bam_suffix(tmp_path: Path):
    # some pipelines produce <stem>.bai instead of <stem>.bam.bai
    t = tmp_path / "t.bam"
    t.write_bytes(b"")
    t.with_suffix(".bai").write_bytes(b"")
    n = _touch_bam(tmp_path, "n.bam")
    sheet = tmp_path / "s.csv"
    _write_csv(sheet, [{"sample_name": "s1", "tumor_bam": str(t), "normal_bam": str(n)}],
               ["sample_name", "tumor_bam", "normal_bam"])
    rows = mod.validate(sheet)
    assert len(rows) == 1
