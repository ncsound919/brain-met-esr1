"""Encode the DAG, derive adjustment sets, enumerate collider paths.

DEPRECATED in favor of dag.R (R dagitty, jtextor/dagitty, verified live:
G1 adjustment set {C, M, T}). Kept because the graph strings are the
specification both implementations share, and for the day a verifiable
pydagitty materializes. No pydagitty API source could be found, so every
call below is trust-me code with loud-degradation guards.
"""

from __future__ import annotations

G1 = """dag {
E [exposure]
S [outcome]
M -> E
M -> S
C -> E
C -> S
T -> E
T -> S
E -> S
}"""

G2 = """dag {
E [exposure]
S [outcome]
M -> E
M -> S
C -> E
C -> S
T -> E
T -> S
E -> S
E -> I
S -> I
}"""


def report(dag_string: str, label: str) -> None:
    from pydagitty import DAG

    dag = DAG(dag_string)
    print(f"--- {label} ---")
    try:
        print("Adjustment sets:", dag.adjustment_sets())
    except Exception as exc:  # keep API drift visible, not silent
        print(f"adjustment_sets() unavailable: {exc}")
    # Path-enumeration helpers differ across pydagitty versions, and the 0.1.x
    # API is unverified in this repo (pydagitty not installed here). Use
    # whichever helper exists; never fail silently, never claim more.
    for method in ("collider_paths", "dconnected_paths", "paths"):
        if hasattr(dag, method):
            print(f"({method} helper present; see DAG plot for enumeration)")
            break
    else:
        print("(no path-enumeration helper; inspect graph manually)")


def main() -> None:
    report(G1, "G1 confounding graph")
    report(G2, "G2 selection-collider graph (I descendant of E, S)")


if __name__ == "__main__":
    main()
