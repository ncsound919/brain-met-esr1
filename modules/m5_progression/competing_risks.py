"""M5: brain-first vs non-CNS-first progression as competing risks.

Math: cause-specific cumulative incidence CIF_k(t) = P(T <= t, cause = k),
Aalen-Johansen estimator: CIF_k += S(t-) * d_kj / n_j summed over distinct
event times, with S the overall (any-event) survival. Treating the other
site as censored (naive KM) overestimates absolute risk -- this module
exists to not do that.

Two layers:
  1. Hand-rolled AJ below (transparent, tested, zero extra deps). Always runs.
  2. comprisk (PyPI, sklearn API) for Gray's K-sample test, Fine-Gray
     regression and CompetingRiskForest when installed. Preprint-stage
     (July 2026): pin the version that works, expect API churn.

Pre-registered: Gray's test on brain-first CIF, ESR1-mut vs wt. A higher
brain-first CIF in wt with overlapping CIs challenges L1; identical CIFs
after adjusting for detection method (Fine-Gray + method covariate)
supports L3 over L1.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.provenance import ProvenanceError, require_real, stamp  # noqa: E402

CENS, BRAIN, NONCNS = "censored", "brain_first", "nonCNS_first"
CODE = {CENS: 0, BRAIN: 1, NONCNS: 2}
SOURCE_TABLE = "first_progression.parquet"


def aalen_johansen(times: np.ndarray, events: np.ndarray, cause: int) -> pd.DataFrame:
    """AJ estimate of CIF for one cause. events: 0=censored, 1..K=causes.

    Returns DataFrame(time, cif, n_risk, n_events_cause).
    """
    order = np.argsort(times, kind="stable")
    times, events = times[order], events[order]
    rows = []
    surv, cif = 1.0, 0.0
    for t in np.unique(times[events > 0]):
        at = times == t
        n = int((times >= t).sum())
        d_all = int((at & (events > 0)).sum())
        d_k = int((at & (events == cause)).sum())
        if n > 0:
            cif += surv * d_k / n
            surv *= 1 - d_all / n
        rows.append({"time": float(t), "cif": cif, "n_risk": n, "n_events_cause": d_k})
    return pd.DataFrame(rows)


def cif_by_group(df: pd.DataFrame, group_col: str, cause: int) -> dict[str, pd.DataFrame]:
    return {
        str(g): aalen_johansen(
            gdf["time_days"].to_numpy(dtype=float),
            gdf["event"].map(CODE).to_numpy(dtype=int),
            cause,
        )
        for g, gdf in df.groupby(group_col)
    }


def try_comprisk(df: pd.DataFrame) -> None:
    """Gray's test via comprisk (verified 0.8.0 API). Skips cleanly if missing."""
    try:
        from comprisk import gray_test
    except ImportError:
        print("comprisk not installed (`pip install comprisk`); AJ estimates above stand alone.")
        return

    res = gray_test(
        df["time_days"].to_numpy(dtype=float),
        df["event"].map(CODE).to_numpy(dtype=int),
        df["esr1_status"].to_numpy(),
        cause=CODE[BRAIN],
    )
    print("Gray's test (brain-first CIF, mut vs wt groups):", res)
    print("Next: Fine-Gray with detection_method covariate; CompetingRiskForest + TreeSHAP.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cause", choices=[BRAIN, NONCNS], default=BRAIN)
    ap.add_argument(
        "--allow-synthetic",
        action="store_true",
        help="plumbing check on placeholder data; output is marked non-promotable",
    )
    args = ap.parse_args()

    # Fail closed. first_progression.parquet has no marker column, so the
    # generator that sets h_brain=0.001 (mut) vs 0.003 (wt) is invisible to any
    # scan of the data -- the manifest is the only thing that knows. Running this
    # on it reproduces the hazard ratio it was built from.
    try:
        statuses = require_real([SOURCE_TABLE], allow_synthetic=args.allow_synthetic)
    except ProvenanceError as exc:
        raise SystemExit(f"PROVENANCE GATE: {exc}")
    synthetic = statuses[SOURCE_TABLE] != "real"

    df = pd.read_parquet(ROOT / "data" / "raw" / SOURCE_TABLE)
    curves = cif_by_group(df, "esr1_status", CODE[args.cause])
    out = ROOT / "results" / "m5"
    out.mkdir(parents=True, exist_ok=True)
    for group, curve in curves.items():
        path = out / f"cif_{args.cause}_{group}.csv"
        curve.to_csv(path, index=False)
        if synthetic:
            stamp(
                path.name,
                "synthetic",
                f"derived from synthetic {SOURCE_TABLE} via --allow-synthetic",
                note=f"cause={args.cause}",
            )
        final = curve["cif"].iloc[-1] if len(curve) else 0.0
        print(f"{args.cause} | {group}: final CIF = {final:.3f} (n={len(curve)} timepoints)")
    if synthetic:
        print(
            "PROVENANCE: synthetic input. These CIFs restate the generator's "
            "hard-coded hazard and are NOT evidence for or against L1."
        )
    try_comprisk(df)


if __name__ == "__main__":
    main()
