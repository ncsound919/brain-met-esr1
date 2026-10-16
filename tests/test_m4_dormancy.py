"""M4: dormancy fallback simulator. No gillespy2 needed."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "m4_dormancy"))

from simulate import fallback_tauleap, summarize


def test_fallback_runs_and_nonnegative():
    a = fallback_tauleap(0.01, 0.0, seed=1, t_end=50)
    assert isinstance(a, int) and a >= 0


def test_pause_targeting_does_not_create_cells_from_nothing():
    # With no switching and no seed actives, actives stay zero regardless.
    assert fallback_tauleap(0.0, 0.1, seed=2, t_end=50) == 0


def test_summarize_groups_correctly():
    import pandas as pd

    df = pd.DataFrame(
        {
            "switching_rate": [0.01, 0.01, 0.02],
            "pause_target_rate": [0.0, 0.0, 0.0],
            "final_active": [10, 20, 5],
        }
    )
    summ = summarize(df)
    assert len(summ) == 2
    assert summ.loc[summ["switching_rate"] == 0.01, "mean_final_active"].iloc[0] == 15.0


def test_model_builder_api():
    pytest.importorskip("gillespy2")
    from model import build_model

    m = build_model(switching_rate=0.01, pause_target_rate=0.05)
    names = (
        [r.name for r in m.listOfReactions.values()]
        if hasattr(m, "listOfReactions")
        else [r.name for r in m.reactions]
    )
    assert "pause_target" in names
