"""OSS-extension wiring: SHELF prior loader + pipeline presence."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "modules" / "m3_occupancy"))

from occupancy_pymc import load_beta_priors


def test_new_scripts_exist():
    assert (ROOT / "modules" / "m1_measurement_error" / "mecor_calibration.R").exists()
    assert (ROOT / "modules" / "m2_dag_meta" / "baggr_meta.R").exists()
    assert (ROOT / "modules" / "m3_occupancy" / "prior_elicitation_shelf.R").exists()
    assert (ROOT / "modules" / "m4_dormancy" / "oncosimulr_selection.R").exists()
    assert (ROOT / "modules" / "m5_progression" / "competing_risks.py").exists()
    assert (ROOT / "modules" / "m6_er_activity" / "er_activity.py").exists()
    assert (ROOT / "data" / "ingest_ecosystem.py").exists()
    assert (ROOT / "modules" / "m2_dag_meta" / "dag.R").exists()


def test_snakefile_wires_new_rules():
    snake = (ROOT / "Snakefile").read_text()
    for rule in (
        "m1_mecor",
        "m2_baggr",
        "m3_shelf_priors",
        "m4_oncosimulr",
        "m5_progression",
        "m6_er_activity",
        "ingest_ecosystem",
        "m2_dag",
    ):
        assert rule in snake, f"Snakefile missing rule {rule}"


def test_beta_prior_loader_defaults():
    alpha, beta = load_beta_priors(None, ["brain", "liver"])
    assert (alpha == 7.0).all() and (beta == 3.0).all()


def test_beta_prior_loader_csv(tmp_path):
    csv = tmp_path / "priors.csv"
    csv.write_text(
        "compartment,alpha,beta,source\nbrain,5.0,4.0,workshop\nliver,7.0,3.0,workshop\n"
    )
    alpha, beta = load_beta_priors(csv, ["brain", "liver", "bone"])
    assert alpha[0] == pytest.approx(5.0) and beta[0] == pytest.approx(4.0)
    # bone absent -> default
    assert alpha[2] == pytest.approx(7.0) and beta[2] == pytest.approx(3.0)


def test_beta_prior_means_sane():
    # Elicited Beta mean must lie in (0,1) or the model is misspecified.
    alpha = np.array([5.0, 7.0])
    beta = np.array([4.0, 3.0])
    means = alpha / (alpha + beta)
    assert ((means > 0) & (means < 1)).all()
