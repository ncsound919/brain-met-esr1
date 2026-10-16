# Module 5 — time to first progression as competing risks

Brain-first vs non-CNS-first progression are mutually exclusive first
events. Standard KM treating the other site as censored overstates
absolute risk; the cause-specific cumulative incidence function (CIF)
is the correct target.

## Math

Aalen-Johansen: `CIF_k(t) = Σ S(t-) · d_kj / n_j` over distinct event
times, `S` = overall survival. `competing_risks.py` implements this by
hand (tested, no new deps) and layers **comprisk** (PyPI, sklearn-style,
July 2026 preprint) for Gray's K-sample test, Fine-Gray regression and
`CompetingRiskForest` + TreeSHAP when installed.

## Pre-registration

Gray's test, brain-first CIF, ESR1-mut vs wt. Higher wt CIF with
overlapping CIs challenges L1; identical CIFs after Fine-Gray adjustment
for detection method supports L3 over L1.

## Data

`first_progression.parquet` (one row per patient; see `FirstProgression`).
Real longitudinal follow-up replaces the synthetic placeholder.

## Dependencies

`pandas`, `numpy`. Optional: `comprisk` (preprint-stage — pin version,
expect API churn; stable fallback: scikit-survival ≥ 0.24
`cumulative_incidence_competing_risks`, GPL-3, nonparametric only).
