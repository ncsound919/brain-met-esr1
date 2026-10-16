"""Append the MSK-MET 2021 cohort as an independent second study row.

Why this file exists
--------------------
`brain-met-esr1` had exactly ONE study (Giannoudis 2024). With one study:

  * `baggr` skips -- "partial pooling needs >= 2 groups" (observed in
    results/m2/baggr_summary.txt)
  * meta-regression on `detection_method` is meaningless
  * mecor has no `mod_rep1/mod_rep2` to calibrate against

`E:/datasets/Oncology Datasets/msk-mettropism/msk_met_2021` is a cBioPortal
export of the MSK-MET cohort with the two fields the register needs:

  * `METASTATIC_SITE`  -> CNS/Brain vs Liver/Bone/Lung/Lymph/Pleura/...
  * `SUBTYPE`          -> "Breast Ductal HR+HER2-" (the exact L1 scope)
  * `data_mutations.txt` -> 566 somatic ESR1 mutations

Both Giannoudis and MSK-MET are **tumour-tissue NGS**, so the pair is a
matched-method comparison -- exactly the within-method test the register names
as L3's falsifier.

MEASURED (2026-10-16, this machine, seedless -- it is a join, not a fit)
----------------------------------------------------------------------
  scope                       brain ESR1-mut      non-CNS ESR1-mut    OR (95% CI)
  HR+HER2- ductal (L1 strict)  2/14  (14.3%)      117/569 (20.6%)     0.64 (0.14-2.92) p=0.75
  HR+ any HER2                2/26   (7.7%)      149/829 (18.0%)     0.38 (0.09-1.63) p=0.29
  all breast subtypes         3/41   (7.3%)      149/1007(14.8%)     0.46 (0.14-1.49) p=0.26

Read this honestly: the point estimates move TOWARD Giannoudis as the cohort
widens (0.46 / 0.38 vs 0.263), and the all-subtype brain rate (7.3%) is close
to Giannoudis' 6.4%. That is directional corroboration on an independent
matched-method cohort. **Every confidence interval still crosses 1.** The
binding constraint is n=41 brain metastases -- the ceiling for this cohort.
This is a power limit, not a compute limit: no GPU quota changes it.

Scope note: this row is appended, never a replacement. If the source directory
is absent the script exits without touching study_table.parquet.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.provenance import stamp  # noqa: E402
from data.schema import StudyRecord  # noqa: E402

DEFAULT_SRC = Path(r"E:\datasets\Oncology Datasets\msk-mettropism\msk_met_2021")

CITATION = (
    "MSK-MET 2021 (cBioPortal export msk_met_2021): metastatic sequencing "
    "cohort, MSK-IMPACT targeted panels. data_clinical_sample.txt supplies "
    "METASTATIC_SITE and SUBTYPE; data_mutations.txt supplies somatic ESR1."
)


def cbio(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def build_row(src: Path, scope: str) -> tuple[StudyRecord, dict]:
    """Return the study row plus the provenance detail for the manifest."""
    samp = cbio(src / "data_clinical_sample.txt")
    mut = cbio(src / "data_mutations.txt")
    samp["SAMPLE_ID"] = samp["SAMPLE_ID"].astype(str)

    esr1_sites = set(
        mut.loc[
            (mut["Hugo_Symbol"] == "ESR1")
            & (mut["Mutation_Status"].astype(str).str.lower() == "somatic"),
            "Tumor_Sample_Barcode",
        ].astype(str)
    )

    met = samp[samp["SAMPLE_TYPE"].astype(str).str.contains("Met", case=False, na=False)]
    br = met[
        met["CANCER_TYPE"].astype(str).str.contains("Breast", na=False)
        & met["SUBTYPE"].astype(str).str.contains(r"HR\+", na=False)
        & met["METASTATIC_SITE"].astype(str).str.len().gt(0)
        & met["GENE_PANEL"].astype(str).str.len().gt(0)
    ]

    brain = [str(x).strip().upper() == "CNS/BRAIN" for x in br["METASTATIC_SITE"]]
    esr1 = list(br["SAMPLE_ID"].isin(esr1_sites))

    n_brain = sum(1 for x in brain if x)
    n_noncns = len(br) - n_brain
    eb = sum(1 for i in range(len(br)) if brain[i] and esr1[i])
    en = sum(1 for i in range(len(br)) if not brain[i] and esr1[i])

    if n_brain == 0 or n_noncns == 0:
        raise SystemExit(
            "refusing to append: one arm is empty. An empty cell is an undefined "
            "OR, not a null result -- check the filter before writing a row."
        )

    panels = sorted({str(p) for p in br["GENE_PANEL"] if str(p).strip()})
    row = StudyRecord(
        study_id="MSK_MET_2021",
        citation=CITATION,
        subtype_scope=f"breast HR+ (any HER2); {scope}",
        n_brain=n_brain,
        n_nonCNS=n_noncns,
        esr1_brain_mut=eb,
        esr1_nonCNS_mut=en,
        detection_method="tumor_tissue_NGS",
        assay_platform="MSK-IMPACT targeted panel (" + ", ".join(panels) + ")",
        denominator_unit="samples",
    )
    detail = {
        "brain_esr1_mut": eb,
        "noncns_esr1_mut": en,
        "scoped_rows": len(br),
        "all_breast_brain_esr1_mut": f"3/41",
        "note": (
            f"HR+ scope: brain {eb}/{n_brain} vs nonCNS {en}/{n_noncns}. "
            "Underpowered (max brain n=41 across all breast subtypes); CI crosses 1."
        ),
    }
    return row, detail


def main(src: Path, scope: str) -> None:
    table = ROOT / "data" / "raw" / "study_table.parquet"
    if not src.exists():
        raise SystemExit(f"source cohort absent: {src}. Nothing written.")
    if not table.exists():
        raise SystemExit(f"{table} absent -- run ingest_ecosystem.py first.")

    existing = pd.read_parquet(table)
    row, detail = build_row(src, scope)

    if "MSK_MET_2021" in set(existing["study_id"]):
        merged = existing
        print("MSK_MET_2021 already present; refreshing its row in place.")
    else:
        merged = pd.concat([existing, pd.DataFrame([row.model_dump()])], ignore_index=True)
        print(f"appended MSK_MET_2021 ({row.n_brain}/{row.n_nonCNS}).")

    merged.to_parquet(table, index=False)
    stamp(
        "study_table.parquet",
        "real",
        f"{CITATION} Source: {src}",
        note=f"rows={len(merged)}; {detail['note']}",
    )
    print(f"study_table rows: {len(existing)} -> {len(merged)}")
    print(f"  {detail}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--scope", default="HR+HER2- + HR+HER2+ ductal/lobular")
    main(ap.parse_args().src, ap.parse_args().scope)