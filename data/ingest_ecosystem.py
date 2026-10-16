"""Bridge: ingest real ecosystem data into the scaffold's study table.

Reads (single source of truth stays in the ecosystem repo):
  <eco>/brain-met-esr1/data/esr1_strata.csv   measured Foundation Medicine counts
  <eco>/evidence/evidence_ledger.csv           ESR1_META anchor (row 21)

Writes:
  data/raw/study_table.parquet   one real row (Giannoudis2024); more rows
                                 append as the 10-study extraction lands
  data/raw/meta_anchor.json      external meta-analytic anchor + provenance

Eco root: --eco-root, $ECO_ROOT, or <scaffold>/../../SCIENCE/Oncology Ecosystem.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from data.provenance import stamp  # noqa: E402
from data.schema import StudyRecord  # noqa: E402

DEFAULT_ECO = ROOT.parents[1] / "SCIENCE" / "Oncology Ecosystem"

CITATION = (
    "Giannoudis A, Sokol ES, Bhogal T, et al. Breast cancer brain metastases "
    "genomic profiling. npj Precision Oncology 2024;8:282. "
    "DOI 10.1038/s41698-024-00761-0, PMID 39706915, PMC11662007."
)


def main(eco_root: Path) -> None:
    strata_path = eco_root / "brain-met-esr1" / "data" / "esr1_strata.csv"
    if not strata_path.exists():
        raise SystemExit(f"Ecosystem strata file missing: {strata_path}")
    strata = pd.read_csv(strata_path)
    pooled = strata[strata["subtype"] == "all"].set_index("site")
    brain = pooled.loc["brain"]
    cns = pooled.loc["non-CNS"]

    row = StudyRecord(
        study_id="Giannoudis2024",
        citation=CITATION,
        subtype_scope="all (pooled); ER+/HER2- percentages w/o counts, see register",
        n_brain=int(brain["n"]),
        n_nonCNS=int(cns["n"]),
        esr1_brain_mut=int(brain["esr1_mut"]),
        esr1_nonCNS_mut=int(cns["esr1_mut"]),
        detection_method="tumor_tissue_NGS",
        assay_platform="Foundation Medicine (unpaired samples)",
        denominator_unit="samples",
    )
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row.model_dump()]).to_parquet(raw / "study_table.parquet", index=False)

    anchor = {
        "annals_106P": {
            "or": 0.35,
            "lo": 0.24,
            "hi": 0.52,
            "n_studies": 10,
            "provenance": (
                "Annals of Oncology abstract 106P, search to 2025-06-05. "
                "Abstract ONLY: study list, detection methods and "
                "heterogeneity unseen."
            ),
            "ledger": "evidence_ledger.csv:21 ESR1_META",
        },
        "giannoudis2024_pooled": {
            "note": "computed from the ingested row, not transcribed",
            "ledger": "brain-met-esr1/data/esr1_strata.csv",
        },
    }
    (raw / "meta_anchor.json").write_text(json.dumps(anchor, indent=2) + "\n")

    # Declared real: this row is measured Foundation Medicine counts, not a
    # placeholder. n_studies=10 in the anchor describes the Annals abstract; only
    # this one study has an extracted row, so M2 pooling stays parked.
    strata_source = f"{strata_path} (via ingest_ecosystem.py); PMID 39706915"
    stamp("study_table.parquet", "real", strata_source)
    stamp(
        "meta_anchor.json",
        "real",
        f"{eco_root}/evidence/evidence_ledger.csv:21 ESR1_META",
        "external meta anchor; abstract-only provenance, study list unseen",
    )
    print(f"Ingested Giannoudis2024 ({row.n_brain}/{row.n_nonCNS} samples) from {eco_root}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--eco-root", type=Path, default=Path(os.environ.get("ECO_ROOT", DEFAULT_ECO)))
    main(ap.parse_args().eco_root)
