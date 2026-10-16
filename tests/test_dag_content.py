"""M2 DAG content guards: edge/node regression tests (no pydagitty needed)."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "modules" / "m2_dag_meta"))

from dag import G1, G2

G1_EDGES = ["M -> E", "M -> S", "C -> E", "C -> S", "T -> E", "T -> S", "E -> S"]


def test_g1_has_all_edges_and_roles():
    for edge in G1_EDGES:
        assert edge in G1, f"G1 missing edge {edge}"
    assert "E [exposure]" in G1
    assert "S [outcome]" in G1


def test_g2_adds_selection_collider_only():
    for edge in G1_EDGES:
        assert edge in G2, f"G2 missing base edge {edge}"
    assert "E -> I" in G2
    assert "S -> I" in G2
    assert "I -> E" not in G2 and "I -> S" not in G2  # I is a consequence, not a cause


def test_graphs_parse_if_pydagitty_present():
    pydagitty = pytest.importorskip("pydagitty")
    for g in (G1, G2):
        dag = pydagitty.DAG(g)  # must not raise on our syntax
        assert dag is not None
