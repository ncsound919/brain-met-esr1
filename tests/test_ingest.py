"""Ingest bridge tests: hermetic fixture + live ecosystem integration."""

import importlib
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ECO_DEFAULT = ROOT.parents[1] / "SCIENCE" / "Oncology Ecosystem"

STRATA = """subtype,site,n,esr1_mut,denominator_unit
all,brain,822,53,samples
all,local,11988,442,samples
all,non-CNS,15516,3220,samples
"""


@pytest.fixture
def eco_root(tmp_path):
    eco = tmp_path / "eco"
    (eco / "brain-met-esr1" / "data").mkdir(parents=True)
    (eco / "brain-met-esr1" / "data" / "esr1_strata.csv").write_text(STRATA)
    return eco


def test_ingest_builds_real_row(eco_root, tmp_path, monkeypatch):
    mod = importlib.import_module("data.ingest_ecosystem")
    fake_root = tmp_path / "scaffold"
    monkeypatch.setattr(mod, "ROOT", fake_root)
    mod.main(eco_root)
    study = pd.read_parquet(fake_root / "data" / "raw" / "study_table.parquet")
    assert len(study) == 1
    row = study.iloc[0]
    assert row["study_id"] == "Giannoudis2024"
    assert (row["n_brain"], row["esr1_brain_mut"]) == (822, 53)
    assert (row["n_nonCNS"], row["esr1_nonCNS_mut"]) == (15516, 3220)
    assert row["denominator_unit"] == "samples"
    anchor = json.loads((fake_root / "data" / "raw" / "meta_anchor.json").read_text())
    assert anchor["annals_106P"]["n_studies"] == 10
    assert anchor["annals_106P"]["or"] == pytest.approx(0.35)


def test_ingest_missing_eco_fails(tmp_path, monkeypatch):
    mod = importlib.import_module("data.ingest_ecosystem")
    monkeypatch.setattr(mod, "ROOT", tmp_path / "scaffold")
    with pytest.raises(SystemExit):
        mod.main(tmp_path / "nope")


@pytest.mark.skipif(
    not (ECO_DEFAULT / "brain-met-esr1" / "data" / "esr1_strata.csv").exists(),
    reason="live ecosystem not mounted",
)
def test_live_pooled_or_is_0263():
    df = pd.read_csv(ECO_DEFAULT / "brain-met-esr1" / "data" / "esr1_strata.csv")
    pooled = df[df["subtype"] == "all"].set_index("site")
    brain = pooled.loc["brain"]
    cns = pooled.loc["non-CNS"]
    a, b = brain["esr1_mut"], brain["n"] - brain["esr1_mut"]
    c, d = cns["esr1_mut"], cns["n"] - cns["esr1_mut"]
    assert (a * d) / (b * c) == pytest.approx(0.263183, abs=1e-4)
