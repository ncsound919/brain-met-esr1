"""PyMC occupancy model with compartment-specific detection.

Marginal likelihood (no discrete latent sampling):
  y_i ~ Bernoulli(psi[comp] * p[comp])

psi = true occupancy per compartment; p = per-survey detection.
Informative Beta priors on p should come from assay sensitivity
validation data. Single-survey estimates are prior-sensitive by
construction -- report the prior alongside the posterior.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.provenance import ProvenanceError, frame_is_synthetic  # noqa: E402
COMPS = ["brain", "liver", "bone", "lung", "lymph_node", "other"]


def load(path: Path) -> pd.DataFrame:
    det = pd.read_parquet(path)
    det = det[det["compartment"].isin(COMPS)].copy()
    return det


def load_beta_priors(
    path: Path | None,
    comps: list[str],
    default: tuple[float, float] = (7.0, 3.0),
    allow_synthetic: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Elicited Beta(a,b) per compartment, aligned to `comps`.

    CSV columns: compartment,alpha,beta[,source]. Compartments missing from
    the CSV fall back to `default`. PLACEHOLDER sources are refused unless
    `allow_synthetic`, which makes the returned priors non-promotable.
    """
    alpha = np.full(len(comps), default[0])
    beta = np.full(len(comps), default[1])
    if path is None:
        return alpha, beta
    priors = pd.read_csv(path)
    if "source" in priors.columns and (priors["source"] == "PLACEHOLDER").any():
        if not allow_synthetic:
            raise ProvenanceError(
                f"elicited priors at {path} contain PLACEHOLDER judgements; pass "
                "allow_synthetic=True for a plumbing check only"
            )
        print("WARNING: elicited priors contain PLACEHOLDER judgements; output is non-promotable.")
    for i, c in enumerate(comps):
        hit = priors[priors["compartment"] == c]
        if len(hit):
            alpha[i], beta[i] = float(hit.iloc[0]["alpha"]), float(hit.iloc[0]["beta"])
    return alpha, beta


def main(
    draws: int = 1000,
    tune: int = 1000,
    beta_priors: Path | None = None,
    sampler: str = "pymc",
    allow_synthetic: bool = False,
) -> None:
    import pymc as pm

    if sampler == "nutpie":
        try:
            import nutpie  # noqa: F401
        except ImportError as exc:
            raise SystemExit(
                "nutpie not installed. Run `pip install brain-met-esr1[fast]` "
                "or use --sampler pymc."
            ) from exc

    # Fail closed on placeholder detection histories. This table DOES carry
    # assay_platform == "synthetic", so the frame check is a genuine second line
    # of defence here rather than the only one.
    det = load(ROOT / "data" / "raw" / "detection_history.parquet")
    if frame_is_synthetic(det) and not allow_synthetic:
        raise SystemExit(
            "PROVENANCE GATE: detection_history.parquet is synthetic. psi/p estimates "
            "from it describe the generator, not assay nondetection. Pass "
            "--allow-synthetic for a plumbing check only."
        )
    comps = sorted(det["compartment"].unique())
    idx = det["compartment"].map({c: i for i, c in enumerate(comps)}).to_numpy()
    y = det["detected"].to_numpy(dtype=int)
    alpha, beta = load_beta_priors(beta_priors, comps, allow_synthetic=allow_synthetic)

    with pm.Model():
        psi = pm.Beta("psi", alpha=2.0, beta=2.0, shape=len(comps))
        # Elicited (SHELF) priors where available, else weakly informative.
        p = pm.Beta("p", alpha=alpha, beta=beta, shape=len(comps))
        pm.Bernoulli("y", p=psi[idx] * p[idx], observed=y)
        idata = pm.sample(draws, tune=tune, target_accept=0.95, nuts_sampler=sampler)

    out = ROOT / "results" / "m3"
    out.mkdir(parents=True, exist_ok=True)
    summ = {
        "compartments": comps,
        "psi_median": np.median(
            idata.posterior["psi"].values.reshape(-1, len(comps)), axis=0
        ).tolist(),
        "p_median": np.median(idata.posterior["p"].values.reshape(-1, len(comps)), axis=0).tolist(),
    }
    pd.DataFrame(summ).to_csv(out / "occupancy_pymc_summary.csv", index=False)
    print(summ)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--tune", type=int, default=1000)
    ap.add_argument("--sampler", choices=["pymc", "nutpie"], default="pymc")
    ap.add_argument(
        "--beta-priors", type=Path, default=None, help="CSV from prior_elicitation_shelf.R"
    )
    ap.add_argument(
        "--allow-synthetic",
        action="store_true",
        help="plumbing check on placeholder data; output is non-promotable",
    )
    args = ap.parse_args()
    main(args.draws, args.tune, args.beta_priors, args.sampler, args.allow_synthetic)
