"""M6: ER-pathway activity by ESR1 status (PROGENy Estrogen footprint + MLM).

Math: for each sample, regress expression y (genes) on footprint weights W:
  beta = (W'W)^-1 W'y,  t = beta / se  (se from residual variance).
Positive t = pathway active, negative = inactive (decoupler MLM logic).
Implemented by hand below (transparent, tested); decoupler cross-checks
when installed.

Pre-registered, both directions:
  ESR1-mut brain mets with LOW Estrogen activity = decoupled ER signaling,
    a mechanistic bridge for L1 (mutation escapes estrogen dependence).
  ESR1-mut with HIGH activity = mutation does not blunt ER output,
    challenging the simple L1 story (then look to L2-L4 + M5 timing).

License note: decoupler's VIPER method is non-commercial; this module uses
mlm/wmean logic only.
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

SOURCES = [
    "expression_samples.parquet",
    "expression.parquet",
    "footprint_estrogen.csv",
]


def mlm_scores(expr: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """MLM activity scores. expr: genes x samples; weights: genes.

    Returns (betas, t-stats) per sample.
    """
    w = weights.reshape(-1, 1)
    # pinv (not inv): footprint weights can be near-collinear; zero weights
    # must degrade to beta=0, not a LinAlgError.
    wt_w_inv = np.linalg.pinv(w.T @ w)
    betas, tstats = np.zeros(expr.shape[1]), np.zeros(expr.shape[1])
    dof = max(expr.shape[0] - 1, 1)
    for s in range(expr.shape[1]):
        y = expr[:, s]
        beta = (wt_w_inv @ w.T @ y)[0]
        resid = y - w[:, 0] * beta
        se = float(np.sqrt((resid @ resid / dof) * wt_w_inv[0, 0]))
        betas[s] = beta
        tstats[s] = beta / se if se > 0 else 0.0
    return betas, tstats


def try_decoupler(matrix: pd.DataFrame, net: pd.DataFrame) -> None:
    """Cross-check with decoupler 2.x (dc.mt.mlm). Skips cleanly otherwise."""
    try:
        import decoupler as dc
    except ImportError:
        print("decoupler not installed (`pip install decoupler`); hand MLM stands alone.")
        return
    try:
        import anndata as ad

        adata = ad.AnnData(matrix.T)
        # 2.x convention: net long-format with source/target[/weight] columns,
        # no column-name kwargs. Estimates land in .obsm; the call returns None.
        dc.mt.mlm(adata, net=net, verbose=False)
        print("decoupler mlm cross-check (obsm keys:", list(adata.obsm.keys()), "):")
        key = "score_mlm" if "score_mlm" in adata.obsm.keys() else None
        if key is not None:
            print(adata.obsm[key].describe().to_string())
        else:
            print("no score key found in obsm")
    except Exception as exc:  # API drift guard: report, don't crash
        print(f"decoupler cross-check skipped ({type(exc).__name__}: {exc})")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=500, help="top-N footprint genes (docs only)")
    ap.add_argument(
        "--allow-synthetic",
        action="store_true",
        help="plumbing check on placeholder data; output is marked non-promotable",
    )
    args = ap.parse_args()

    # Fail closed. The placeholder generator writes shift=1.5 for ESR1-mut, i.e.
    # it encodes the *falsifying* case (mut brain mets with ER activity as high
    # as wt) by construction, and expression.parquet carries no marker column.
    try:
        statuses = require_real(SOURCES, allow_synthetic=args.allow_synthetic)
    except ProvenanceError as exc:
        raise SystemExit(f"PROVENANCE GATE: {exc}")
    synthetic = any(s != "real" for s in statuses.values())

    samples = pd.read_parquet(ROOT / "data" / "raw" / "expression_samples.parquet")
    expr = pd.read_parquet(ROOT / "data" / "raw" / "expression.parquet")
    foot = pd.read_csv(ROOT / "data" / "raw" / "footprint_estrogen.csv")
    if (foot["source"] == "PLACEHOLDER-synthetic").any():
        print("WARNING: footprint weights are synthetic placeholders.")

    matrix = expr.pivot(index="gene", columns="sample_id", values="log_norm_count")
    genes = [g for g in foot["gene"] if g in matrix.index]
    w = foot.set_index("gene").loc[genes, "weight"].to_numpy(dtype=float)
    betas, tstats = mlm_scores(matrix.loc[genes].to_numpy(dtype=float), w)

    res = samples.copy()
    res["er_beta"] = betas
    res["er_t"] = tstats
    out = ROOT / "results" / "m6"
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "er_activity.csv", index=False)
    if synthetic:
        stamp(
            "er_activity.csv",
            "synthetic",
            f"derived from synthetic {[s for s, v in statuses.items() if v != 'real']} "
            "via --allow-synthetic",
        )
        print(
            "PROVENANCE: synthetic input. These ER-activity scores are generated by "
            "the same placeholder that hard-codes the mut/wt shift, so they cannot "
            "test M6 in either direction."
        )
    print(res.groupby(["esr1_status", "compartment"])[["er_beta", "er_t"]].mean().to_string())

    net = pd.DataFrame({"source": "Estrogen", "target": foot["gene"], "weight": foot["weight"]})
    try_decoupler(matrix, net)


if __name__ == "__main__":
    main()
