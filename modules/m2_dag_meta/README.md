# Module 2 — DAG + meta-regression

## What the DAG encodes

Detection method (M) and tumor fraction (C) both affect measured ESR1 (E)
and site composition (S). In the base graph they are **confounders**,
not colliders. The "ctDNA studies open a spurious path" claim requires a
second graph where study inclusion (I) is a descendant of both E and S
(selection collider). `dag.py` encodes **both** graphs so the adjustment
sets can differ. The point is to make identification assumptions explicit,
not to assert them.

## Pre-registration requirement

Before fitting `meta_regression.R`, write down which strata you predict
will flip and in which direction. The meta-regression can show method is
a significant moderator (Q_M); it cannot show the DAG is correct.
The DAG gives the interpretation; the regression gives the evidence.

## Dependencies

- `dag.py`: `pydagitty` 0.1.x, DAG/ADMG only (PAG behavior experimental).
- `meta_regression.R`: `metafor`, `arrow`, `dplyr`.
- `baggr_meta.R`: `baggr` v0.8.2 (GPL-3) + `rstan` ≥ 2.26 C++ toolchain.

## Extension: baggr partial pooling

One cohort plus one abstract-only meta is no basis for REML moderation
— today the table holds Giannoudis2024 alone (single method), so
`baggr_meta.R` and `meta_regression.R` run pooled-only by construction
and go live when study 2 lands. `baggr_meta.R` fits partial-pooling
(plus full-pooling and method-moderated) Stan models, adjudicates pooling
via `baggr_compare`, and offers opt-in leave-one-study-out CV
(`BAGGR_LOO=1`; `BAGGR_ITER=500` for smoke runs only). Upgrade path to a
binomial-likelihood model on raw 2x2 counts (no continuity correction) is
sketched in-script.
