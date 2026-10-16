"""Property tests: invariants that must hold for ALL inputs, not just examples."""

import math
import sys
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from pydantic import ValidationError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.schema import DetectionHistory, StudyRecord

# Lean profile: invariants here are cheap, 20 examples each keeps CI snappy.
settings.register_profile("lean", max_examples=20)
settings.load_profile("lean")

COMPARTMENTS = ["brain", "liver", "bone", "lung", "lymph_node", "other"]
CTDNA_METHODS = ["ctDNA_NGS", "ctDNA_ddPCR"]
TISSUE_METHODS = ["tumor_tissue_NGS", "WES", "IHC", "other"]


def det_row(**over):
    row = {
        "patient_id": "P001",
        "compartment": "liver",
        "timepoint_days": 0,
        "detection_method": "ctDNA_NGS",
        "assay_platform": "x",
        "limit_of_detection": 0.01,
        "detected": True,
        "ctdna_fraction": 0.02,
        "denominator_unit": "patients",
    }
    row.update(over)
    return row


@given(method=st.sampled_from(CTDNA_METHODS))
def test_ctdna_always_requires_fraction(method):
    with pytest.raises(ValidationError):
        DetectionHistory(**det_row(detection_method=method, ctdna_fraction=None))


@given(method=st.sampled_from(TISSUE_METHODS))
def test_tissue_never_requires_fraction(method):
    DetectionHistory(**det_row(detection_method=method, ctdna_fraction=None))


@given(
    n_brain=st.integers(0, 500),
    n_cns=st.integers(0, 500),
    over=st.integers(1, 50),
)
def test_counts_above_n_always_fail(n_brain, n_cns, over):
    base = {
        "study_id": "S",
        "citation": "c",
        "subtype_scope": "ER+/HER2-",
        "detection_method": "WES",
        "assay_platform": "x",
        "denominator_unit": "patients",
    }
    with pytest.raises(ValidationError):
        StudyRecord(
            **base,
            n_brain=n_brain,
            n_nonCNS=n_cns,
            esr1_brain_mut=n_brain + over,
            esr1_nonCNS_mut=0,
        )


@given(days=st.integers(min_value=0, max_value=36500))
def test_nonnegative_timepoint_always_valid(days):
    DetectionHistory(
        **det_row(detection_method="tumor_tissue_NGS", ctdna_fraction=None, timepoint_days=days)
    )


@given(days=st.integers(max_value=-1))
def test_negative_timepoint_always_fails(days):
    with pytest.raises(ValidationError):
        DetectionHistory(
            **det_row(detection_method="tumor_tissue_NGS", ctdna_fraction=None, timepoint_days=days)
        )


@given(
    lam1=st.floats(0.1, 0.99),
    lam2=st.floats(0.1, 0.99),
)
def test_attenuation_monotonic_in_lambda(lam1, lam2):
    # Weaker reliability -> correction farther from 1 (OR_obs < 1 here).
    c1 = math.exp(math.log(0.35) / lam1)
    c2 = math.exp(math.log(0.35) / lam2)
    if lam1 < lam2:
        assert c1 < c2
    elif lam2 < lam1:
        assert c2 < c1
    else:
        assert c1 == c2


@given(lam=st.floats(0.1, 1.0))
def test_corrected_interval_brackets_corrected_or(lam):
    lo = math.exp(math.log(0.24) / lam)
    mid = math.exp(math.log(0.35) / lam)
    hi = math.exp(math.log(0.52) / lam)
    assert lo < mid < hi
