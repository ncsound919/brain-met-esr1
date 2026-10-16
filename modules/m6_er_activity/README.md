# Module 6 — ER-pathway activity (functional readout for ESR1 status)

Mutation status alone doesn't say whether ER signaling still runs.
M6 scores PROGENy Estrogen-footprint activity per sample (MLM t-statistic)
and compares ESR1-mut vs wt × compartment.

## Math

Per sample: `β = (W'W)⁻¹W'y`, `t = β/se`. Positive = active, negative =
inactive. Hand-rolled in `er_activity.py` (tested); `decoupler` 2.x
(scverse, Saez lab; `dc.mt.mlm`, verified against 2.2.0) cross-checks when
installed. Footprint weights are **vendored** (`footprint_estrogen.csv`,
via `dc.op.progeny(organism='human', top=500)` once) — never
re-downloaded per run. VIPER excluded (non-commercial license); mlm/wmean
logic only.

## Pre-registration

Mut+brain LOW activity = decoupled signaling (mechanistic bridge for L1).
Mut HIGH activity = ER output intact (challenges simple L1; look to
L2–L4, M5). Both directions stated before running.

## Data

`expression_samples.parquet` + `expression.parquet` (tidy matrix).
Real full-transcriptome RNA replaces the 200-gene synthetic placeholder.

## Dependencies

`pandas`, `numpy`. Optional: `decoupler`, `anndata` (cross-check only).
