"""M5: Aalen-Johansen math. Pure NumPy/pandas, no comprisk needed."""

import sys
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "modules" / "m5_progression"))
sys.path.insert(0, str(ROOT))

from competing_risks import CODE, aalen_johansen

from data.schema import FirstProgression


def test_hand_checkable_case():
    # times [5,5,10], events [brain, censored, brain]:
    # t=5: n=3, d_brain=1 -> cif=1/3, S=2/3; t=10: n=1 -> cif=1/3+2/3=1.0
    times = np.array([5.0, 5.0, 10.0])
    events = np.array([CODE["brain_first"], CODE["censored"], CODE["brain_first"]])
    curve = aalen_johansen(times, events, CODE["brain_first"])
    assert curve["cif"].iloc[0] == pytest.approx(1 / 3)
    assert curve["cif"].iloc[-1] == pytest.approx(1.0)


def test_cif_bounded_and_monotone():
    rng = np.random.default_rng(3)
    times = rng.uniform(1, 100, size=300)
    events = rng.choice([0, 1, 2], size=300, p=[0.4, 0.3, 0.3])
    for cause in (1, 2):
        cif = aalen_johansen(times, events, cause)["cif"].to_numpy()
        assert ((cif >= 0) & (cif <= 1)).all()
        assert (np.diff(cif) >= -1e-12).all()


def test_competing_causes_sum_to_one_minus_survival():
    # No censoring: CIF_1 + CIF_2 == 1 at the end, and == 1 - S throughout.
    times = np.array([1.0, 2.0, 3.0, 4.0])
    events = np.array([1, 2, 1, 2])
    c1 = aalen_johansen(times, events, 1)["cif"].iloc[-1]
    c2 = aalen_johansen(times, events, 2)["cif"].iloc[-1]
    assert c1 + c2 == pytest.approx(1.0)


def test_progression_schema():
    FirstProgression(
        patient_id="P1",
        esr1_status="mut",
        time_days=100,
        event="brain_first",
        detection_method="WES",
    )
    with pytest.raises(ValidationError):
        FirstProgression(
            patient_id="P1",
            esr1_status="mut",
            time_days=-5,
            event="brain_first",
            detection_method="WES",
        )
