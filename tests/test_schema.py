"""Schema tests: no heavy deps beyond pydantic/pandas."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.schema import DetectionHistory, StudyRecord, check_duplicates


def base_det(**over):
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


def test_ctdna_requires_fraction():
    with pytest.raises(ValueError, match="ctdna_fraction"):
        DetectionHistory(**base_det(ctdna_fraction=None))
    # tissue rows may omit it
    DetectionHistory(**base_det(detection_method="tumor_tissue_NGS", ctdna_fraction=None))


def test_study_counts():
    StudyRecord(
        study_id="S1",
        citation="c",
        subtype_scope="ER+/HER2-",
        n_brain=10,
        n_nonCNS=20,
        esr1_brain_mut=1,
        esr1_nonCNS_mut=4,
        detection_method="WES",
        assay_platform="x",
        denominator_unit="patients",
    )
    with pytest.raises(ValueError, match="exceeds n_brain"):
        StudyRecord(
            study_id="S2",
            citation="c",
            subtype_scope="ER+/HER2-",
            n_brain=10,
            n_nonCNS=20,
            esr1_brain_mut=11,
            esr1_nonCNS_mut=4,
            detection_method="WES",
            assay_platform="x",
            denominator_unit="patients",
        )


def test_contradictory_duplicates_fail():
    rows = [
        base_det(
            patient_id="P1",
            compartment="brain",
            detection_method="tumor_tissue_NGS",
            ctdna_fraction=None,
            detected=True,
        ),
        base_det(
            patient_id="P1",
            compartment="brain",
            detection_method="tumor_tissue_NGS",
            ctdna_fraction=None,
            detected=False,
        ),
    ]
    with pytest.raises(ValueError, match="Contradictory duplicate"):
        check_duplicates(rows)
