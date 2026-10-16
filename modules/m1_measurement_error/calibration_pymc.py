"""Bayesian regression calibration for lambda uncertainty.

Reads data/raw/study_table.parquet. Produces a posterior on OR_true.
Use instead of (not alongside) the SIMEX grid when you want a posterior
on lambda propagated into the corrected estimate.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load_effects(path: Path) -> tuple[np.ndarray, np.ndarray]:
    study = pd.read_parquet(path)
    a = study["esr1_brain_mut"] + 0.5
    b = study["n_brain"] - study["esr1_brain_mut"] + 0.5
    c = study["esr1_nonCNS_mut"] + 0.5
    d = study["n_nonCNS"] - study["esr1_nonCNS_mut"] + 0.5
    log_or = np.log((a * d) / (b * c)).to_numpy(dtype=float)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d).to_numpy(dtype=float)
    return log_or, se


def main(draws: int = 2000, tune: int = 1000, sampler: str = "pymc") -> None:
    import arviz as az
    import pymc as pm

    if sampler == "nutpie":
        try:
            import nutpie  # noqa: F401
        except ImportError as exc:
            raise SystemExit(
                "nutpie not installed. Run `pip install brain-met-esr1[fast]` "
                "or use --sampler pymc."
            ) from exc

    log_or, se = load_effects(ROOT / "data" / "raw" / "study_table.parquet")
    with pm.Model():
        lam = pm.Beta("lambda", alpha=2.0, beta=1.0)
        true_log_or = pm.Normal("true_log_or", mu=0.0, sigma=1.0)
        pm.Normal("obs", mu=lam * true_log_or, sigma=se, observed=log_or)
        # nutpie = Rust NUTS via nuts-rs; faster ESS/grad, same posterior.
        # Requires a fully continuous model -- holds here (all free RVs continuous).
        idata = pm.sample(draws, tune=tune, target_accept=0.95, nuts_sampler=sampler)

    posterior_or = np.exp(idata.posterior["true_log_or"].values)
    print(f"Posterior median OR_true: {np.median(posterior_or):.3f}")
    print(f"94% HDI: {az.hdi(posterior_or, hdi_prob=0.94)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--tune", type=int, default=1000)
    ap.add_argument("--sampler", choices=["pymc", "nutpie"], default="pymc")
    args = ap.parse_args()
    main(args.draws, args.tune, args.sampler)
