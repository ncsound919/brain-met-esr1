# data/

Five tables, all Parquet (plus one footprint CSV), all validated against
`schema.py` before any module runs.

## detection_history.parquet

One row per patient, compartment, timepoint, and method.
See `DetectionHistory` in `schema.py`.

## study_table.parquet

One row per included study (see `StudyRecord`). Owned by
`ingest_ecosystem.py`: row 1 is the real Giannoudis2024 cohort ingested
from the ecosystem strata CSV. Rows append as the 10-study extraction
lands (study list currently unseen — abstract only).

## meta_anchor.json

External meta-analytic anchor (Annals of Oncology abstract 106P: OR 0.35,
10 studies, abstract-only provenance). Owned by `ingest_ecosystem.py`;
M1 reads it, nothing else may.

## first_progression.parquet (M5)

One row per patient: time to first progression site or censoring.
See `FirstProgression`. Brain-first vs non-CNS-first are competing events.

## expression_samples.parquet + expression.parquet (M6)

Sample metadata (`ExpressionSample`) plus a tidy gene-expression matrix
(`ExpressionRow`: sample_id, gene, log_norm_count). Full-transcriptome
scale in production; 200 genes synthetic.

## footprint_estrogen.csv (M6)

PROGENy-style Estrogen footprint weights (gene, weight, p_value).
Production: vendor once via `decoupler.get_progeny(organism='human',
top=500)` and commit — never re-download per run (reproducibility).

Placeholders: `make_synthetic.py` generates schema-valid synthetic tables
for everything ingest does not yet cover, so the pipeline and tests run.
Numbers from synthetic data are plumbing checks, not findings.

## What "ambiguous" means

The loader (`validate.py`) fails if:

- a required column is missing;
- a literal-typed field has a value outside its allowed set;
- `esr1_mut > n`;
- a ctDNA row lacks `ctdna_fraction`;
- the same `(patient_id, compartment, timepoint_days, detection_method)`
  appears twice with different `detected` values.
