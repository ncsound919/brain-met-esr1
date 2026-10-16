"""M3: marginal occupancy likelihood. Pure NumPy, no PyMC needed."""

import numpy as np
import pytest


def marginal_prob(psi: float, p: float) -> float:
    return psi * p


def test_nondetection_bias_direction():
    # Same true occupancy, lower brain detection -> lower naive prevalence.
    psi = 0.2
    assert marginal_prob(psi, 0.3) < marginal_prob(psi, 0.7)


def test_single_survey_nonidentifiable():
    # psi*p products collide: (0.5,0.4) and (0.4,0.5) give same marginal.
    assert marginal_prob(0.5, 0.4) == marginal_prob(0.4, 0.5)
    # ... so single-survey data cannot separate them without priors.
    assert True  # documents the README identifiability warning


def test_mle_recovers_product():
    rng = np.random.default_rng(0)
    psi, p = 0.4, 0.6
    y = rng.binomial(1, psi * p, size=20000)
    assert y.mean() == pytest.approx(psi * p, abs=0.02)
