"""Generate schema-valid synthetic placeholders so the pipeline runs pre-extraction.

Real data replaces these files 1:1. Do NOT tune methods on this file and
report the numbers as findings.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
rng = np.random.default_rng(7)

sys.path.insert(0, str(ROOT))
from data.provenance import stamp  # noqa: E402

SYNTHETIC_SOURCE = (
    "data/make_synthetic.py (numpy default_rng(7)); first_progression and the "
    "expression matrix encode the hypothesis by construction"
)


def make_detection_history(n_patients: int = 60) -> pd.DataFrame:
    rows = []
    for p in range(n_patients):
        for comp in ("brain", "liver", "bone"):
            method = "ctDNA_NGS" if comp != "brain" else "tumor_tissue_NGS"
            detected = bool(rng.random() < (0.08 if comp == "brain" else 0.22))
            rows.append(
                {
                    "patient_id": f"P{p:03d}",
                    "compartment": comp,
                    "timepoint_days": 0,
                    "detection_method": method,
                    "assay_platform": "synthetic",
                    "limit_of_detection": 0.01,
                    "detected": detected,
                    "ctdna_fraction": (
                        float(rng.uniform(0.001, 0.05)) if "ctDNA" in method else None
                    ),
                    "denominator_unit": "patients",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    # NOTE: study_table.parquet is owned by ingest_ecosystem.py (real rows).
    # This generator covers the still-synthetic tables only.
    make_detection_history().to_parquet(RAW / "detection_history.parquet", index=False)
    make_first_progression().to_parquet(RAW / "first_progression.parquet", index=False)
    samples, expr = make_expression()
    samples.to_parquet(RAW / "expression_samples.parquet", index=False)
    expr.to_parquet(RAW / "expression.parquet", index=False)
    make_footprint().to_csv(RAW / "footprint_estrogen.csv", index=False)
    # Declare provenance. first_progression.parquet and expression.parquet carry
    # no marker column, so a scan of the data cannot detect them -- only this can.
    for table in (
        "detection_history.parquet",
        "first_progression.parquet",
        "expression_samples.parquet",
        "expression.parquet",
        "footprint_estrogen.csv",
    ):
        stamp(table, "synthetic", SYNTHETIC_SOURCE)
    print(f"Wrote synthetic placeholders to {RAW}")


def make_first_progression(n_patients: int = 120) -> pd.DataFrame:
    """Synthetic competing events. ESR1-mut has lower brain-first hazard,
    mirroring the measured pooled OR 0.26 anchor. PLACEHOLDER."""
    rows = []
    for p in range(n_patients):
        esr1 = "mut" if rng.random() < 0.25 else "wt"
        # Cause-specific hazards (per-day, toy scale).
        h_brain = 0.001 if esr1 == "mut" else 0.003
        h_cns = 0.004
        t_brain = rng.exponential(1 / h_brain)
        t_cns = rng.exponential(1 / h_cns)
        t_cens = rng.exponential(1 / 0.002)
        t = min(t_brain, t_cns, t_cens)
        event = "censored" if t == t_cens else ("brain_first" if t == t_brain else "nonCNS_first")
        rows.append(
            {
                "patient_id": f"P{p:03d}",
                "esr1_status": esr1,
                "time_days": int(min(t, 1825)),
                "event": event if t <= 1825 else "censored",
                "detection_method": "tumor_tissue_NGS",
                "subtype_scope": "ER+/HER2-",
            }
        )
    return pd.DataFrame(rows)


ESTROGEN_GENES = [f"ESRTARGET{i:03d}" for i in range(40)]


def make_expression(n_genes: int = 200, n_per_group: int = 6) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Tiny synthetic RNA: ESR1-mut samples overexpress estrogen-footprint genes
    (constitutive ER activity). Real data replaces 1:1 at full transcriptome scale."""
    genes = ESTROGEN_GENES + [f"GENE{i:04d}" for i in range(n_genes - len(ESTROGEN_GENES))]
    samples, expr = [], []
    sid = 0
    for esr1 in ("mut", "wt"):
        for _ in range(n_per_group):
            sid += 1
            sample_id = f"S{sid:03d}"
            samples.append(
                {
                    "sample_id": sample_id,
                    "patient_id": f"P{sid:03d}",
                    "compartment": "brain" if sid % 2 else "liver",
                    "esr1_status": esr1,
                    "assay_platform": "synthetic",
                }
            )
            shift = 1.5 if esr1 == "mut" else 0.0
            for g in genes:
                mu = shift if g in ESTROGEN_GENES else 0.0
                expr.append(
                    {
                        "sample_id": sample_id,
                        "gene": g,
                        "log_norm_count": float(rng.normal(mu, 1.0)),
                    }
                )
    return pd.DataFrame(samples), pd.DataFrame(expr)


def make_footprint() -> pd.DataFrame:
    """Synthetic PROGENy-style Estrogen footprint. Production: vendor the real
    weights via decoupler.get_progeny(organism='human', top=500) once, store here."""
    return pd.DataFrame(
        {
            "gene": ESTROGEN_GENES,
            "weight": [float(rng.normal(2.0, 0.5)) for _ in ESTROGEN_GENES],
            "p_value": [1e-6] * len(ESTROGEN_GENES),
            "source": ["PLACEHOLDER-synthetic"] * len(ESTROGEN_GENES),
        }
    )


if __name__ == "__main__":
    main()
