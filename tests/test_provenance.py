"""Provenance gate tests.

The load-bearing one is `test_column_scan_cannot_detect_first_progression`:
without it the gate could be satisfied by a scan of `assay_platform` and the
manifest would be decoration. Every other test asserts the gate bites.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@pytest.fixture
def prov(tmp_path, monkeypatch):
    """Fresh provenance module bound to an isolated manifest."""
    mod = importlib.import_module("data.provenance")
    manifest = tmp_path / "provenance.json"
    monkeypatch.setattr(mod, "MANIFEST", manifest)
    return mod


def test_absent_manifest_makes_everything_unknown(prov):
    assert not prov.MANIFEST.exists()
    assert prov.table_status("first_progression.parquet") == "unknown"


def test_undeclared_table_is_not_promotable(prov):
    """Fail closed: no declaration means no promotion, not an implicit pass."""
    with pytest.raises(prov.ProvenanceError, match="never declared"):
        prov.require_real(["first_progression.parquet"])


def test_synthetic_table_blocks(prov):
    prov.stamp("t.parquet", "synthetic", "make_synthetic.py")
    with pytest.raises(prov.ProvenanceError, match="not-real|non-real"):
        prov.require_real(["t.parquet"])


def test_real_table_passes_and_returns_source(prov):
    prov.stamp("study_table.parquet", "real", "esr1_strata.csv PMID 39706915")
    out = prov.require_real(["study_table.parquet"])
    assert out == {"study_table.parquet": "real"}


def test_allow_synthetic_reports_status_instead_of_raising(prov):
    prov.stamp("t.parquet", "synthetic", "make_synthetic.py")
    out = prov.require_real(["t.parquet"], allow_synthetic=True)
    assert out["t.parquet"] == "synthetic"


def test_multiple_blockers_are_all_named(prov):
    """A gate that reports one offender hides the others."""
    prov.stamp("a.parquet", "synthetic", "gen")
    prov.stamp("b.parquet", "real", "real source")
    with pytest.raises(prov.ProvenanceError) as exc:
        prov.require_real(["a.parquet", "b.parquet", "c.parquet"])
    msg = str(exc.value)
    assert "a.parquet" in msg and "c.parquet" in msg
    assert "b.parquet" not in msg


def test_stamp_rejects_empty_source(prov):
    with pytest.raises(ValueError, match="empty source"):
        prov.stamp("t.parquet", "real", "   ")


def test_stamp_rejects_unknown_status(prov):
    with pytest.raises(ValueError, match="must be 'real' or 'synthetic'"):
        prov.stamp("t.parquet", "unknown", "src")  # type: ignore[arg-type]


def test_stamp_persists_to_disk(prov):
    """The declaration must survive the process, not live in memory."""
    import json

    prov.stamp("t.parquet", "real", "https://example.test/d")
    on_disk = json.loads(prov.MANIFEST.read_text(encoding="utf-8"))
    assert on_disk["tables"]["t.parquet"]["status"] == "real"
    assert on_disk["tables"]["t.parquet"]["source"] == "https://example.test/d"


# --- the negative control -------------------------------------------------


def test_column_scan_cannot_detect_first_progression():
    """`frame_is_synthetic` must FAIL on the one table that matters most.

    first_progression.parquet is written by make_synthetic.py with a hard-coded
    3x hazard ratio and carries no marker column, so this returns False. If it
    ever returned True the manifest would be redundant and this test should be
    revisited -- but until first_progression gains provenance the gate has to
    consult the manifest, not the frame.
    """
    frame = pd.DataFrame(
        {
            "patient_id": ["P000", "P001"],
            "esr1_status": ["mut", "wt"],
            "time_days": [100, 200],
            "event": ["brain_first", "nonCNS_first"],
            "detection_method": ["tumor_tissue_NGS"] * 2,
            "subtype_scope": ["ER+/HER2-"] * 2,
        }
    )
    assert "assay_platform" not in frame.columns
    assert prov_module().frame_is_synthetic(frame) is False


def test_column_scan_does_detect_marked_tables():
    frame = pd.DataFrame({"assay_platform": ["synthetic"]})
    assert prov_module().frame_is_synthetic(frame) is True


def prov_module():
    return importlib.import_module("data.provenance")