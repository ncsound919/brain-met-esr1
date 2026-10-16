"""M1: attenuation identity. Pure math, no R/PyMC needed."""

import math

import pytest

OBS, LO, HI = 0.35, 0.24, 0.52  # external annals_106P anchor (abstract-only)
POOLED = 0.263183  # measured Giannoudis2024 pooled OR (Fisher, live-verified)


def corrected(or_obs: float, lam: float) -> float:
    return math.exp(math.log(or_obs) / lam)


def test_lambda_grid_reference_values():
    # Measured: exp(log(x)/lam). lam=0.6 -> 0.093-0.336; lam=0.8 -> 0.168-0.442
    assert corrected(OBS, 0.6) == pytest.approx(0.174, abs=0.005)
    assert corrected(LO, 0.6) == pytest.approx(0.093, abs=0.005)
    assert corrected(HI, 0.6) == pytest.approx(0.336, abs=0.005)
    assert corrected(LO, 0.8) == pytest.approx(0.168, abs=0.005)
    assert corrected(HI, 0.8) == pytest.approx(0.442, abs=0.005)


def test_correction_direction():
    # attenuation: corrected OR is farther from 1 than observed (OR<1 here)
    for lam in (0.6, 0.7, 0.8, 0.9):
        assert corrected(OBS, lam) < OBS
    # weaker reliability -> stronger correction
    assert corrected(OBS, 0.6) < corrected(OBS, 0.9)


def test_pooled_anchor_reference_values():
    # Measured: exp(log(0.263183)/lam). lam=0.6 -> 0.108; lam=0.8 -> 0.189
    assert corrected(POOLED, 0.6) == pytest.approx(0.108, abs=0.005)
    assert corrected(POOLED, 0.8) == pytest.approx(0.189, abs=0.005)
