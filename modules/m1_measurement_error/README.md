# Module 1 — measurement error

Dual-anchored to two depletion estimates (see register for provenance):

- **pooled measured**: OR 0.263, Giannoudis2024 single cohort
  (data-driven; recomputed every run from `study_table.parquet`);
- **external meta**: OR 0.35 (0.24–0.52), Annals of Oncology abstract 106P,
  10 studies, abstract ONLY (methods unseen).

Attenuation identity per anchor: `λ · log(OR_true) = log(OR_obs)`,
so `OR_true = exp(log(OR_obs) / λ)`.
Reference points (measured): pooled λ = 0.6 → ≈ 0.11, λ = 0.8 → ≈ 0.19;
meta λ = 0.6 → ≈ 0.17, λ = 0.8 → ≈ 0.27.

## When to use which script

- `simex_correction.R` — fast point correction + sensitivity grid over λ.
  Build this first.
- `calibration_pymc.py` — posterior on λ and OR_true. Build only if you
  can afford sampling and want λ uncertainty propagated.

Do not build both in parallel; SIMEX first, PyMC only if needed.

## Design note

SIMEX assumes **nondifferential** measurement error. If brain biopsies
have systematically higher tumor purity than plasma ctDNA, the error is
differential and the correction **undercorrects**. Read each λ-grid output
as a sensitivity analysis, not a single corrected estimate. Report the
full grid plus interval width at each λ, per anchor.

## R dependencies

`metafor`, `arrow`, `dplyr`, `simex`, `readr`, `jsonlite`, `mecor` (v1.0.0, GPL-3).

## Extension: mecor regression calibration

`mecor_calibration.R` corrects an **error-prone continuous moderator**
(e.g. study-mean tumor purity) in the `log_or ~ moderator` meta-regression,
via regression calibration with bootstrap CIs — plus differential-error
methods for the brain-purity-vs-plasma worry. Needs replicate measures
(`mod_value`, `mod_rep1`, `mod_rep2`); without them it writes a collection
template and a reliability sensitivity grid. Complements, not replaces,
the SIMEX grid on the pooled OR.
