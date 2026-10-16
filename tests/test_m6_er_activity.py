"""M6: MLM activity math. Pure NumPy, no decoupler needed."""

import sys
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "modules" / "m6_er_activity"))
sys.path.insert(0, str(ROOT))

from er_activity import mlm_scores

from data.schema import ExpressionRow, ExpressionSample


def test_mlm_recovers_known_beta():
    rng = np.random.default_rng(5)
    w = rng.normal(2.0, 0.5, size=100)
    expr = np.column_stack([2.0 * w + rng.normal(0, 0.1, size=100) for _ in range(5)])
    betas, tstats = mlm_scores(expr, w)
    assert (np.abs(betas - 2.0) < 0.1).all()
    assert (tstats > 5).all()  # strong positive activity


def test_mlm_sign_follows_direction():
    rng = np.random.default_rng(6)
    w = rng.normal(0, 1.0, size=80)
    up = mlm_scores((1.5 * w + rng.normal(0, 0.2, 80)).reshape(-1, 1), w)
    down = mlm_scores((-1.5 * w + rng.normal(0, 0.2, 80)).reshape(-1, 1), w)
    assert up[1][0] > 0 and down[1][0] < 0


def test_mlm_zero_weights_safe():
    y = np.ones((10, 2))
    betas, tstats = mlm_scores(y, np.zeros(10))
    assert np.isfinite(betas).all() and np.isfinite(tstats).all()


def test_expression_schema():
    s = ExpressionSample(
        sample_id="S1",
        patient_id="P1",
        compartment="brain",
        esr1_status="wt",
        assay_platform="x",
    )
    ExpressionRow(sample_id=s.sample_id, gene="ESR1", log_norm_count=3.2)
    with pytest.raises(ValidationError):
        ExpressionSample(
            sample_id="S2",
            patient_id="P2",
            compartment="mars",
            esr1_status="wt",
            assay_platform="x",
        )
