# brain-met-esr1

ESR1-mutation organotropism in breast-cancer brain metastasis: is the
depletion real, or measurement, confounding, detection, or dormancy?

Observed anchors (both depletion; magnitude differs — see register):
measured pooled OR **0.263** (53/822 brain vs 3220/15516 non-CNS,
Fisher p ≈ 7e-29; Giannoudis 2024, Foundation Medicine, `samples`
denominator) plus external meta OR **0.35 (0.24–0.52), 10 studies**
(Annals of Oncology abstract 106P — abstract ONLY, methods unseen).
Naive pooled contrast 6.4% vs 20.7% is consistent with L1–L4
simultaneously — do not read it as evidence for L1.

## Ecosystem wiring

`data/ingest_ecosystem.py` reads the single source of truth in
`SCIENCE/Oncology Ecosystem` (strata CSV + `ESR1_META` ledger anchor)
into `study_table.parquet` + `meta_anchor.json`
(`--eco-root` / `$ECO_ROOT` override the default relative path).
L1/L2 verdicts are imported from the ecosystem register, not re-derived.

## Structure

Six modules, one data spine (`data/schema.py`), Snakemake pipeline.
Each module's README states what it can and cannot show.

| Step | Module | Dependency | Why |
|------|--------|-----------|-----|
| 1 | Schema + validation + ingest | Ecosystem repo present | Real rows before anything reads them |
| 2 | M1 SIMEX (dual anchor) | Study table + meta anchor | Cheapest, data-driven pooled + external 0.35 |
| 3 | M3 occupancy | Detection histories | Same schema, tests nondetection claim |
| 4 | M2 DAG + meta-regression | Study table | Pooled-only until study 2 lands |
| 5 | M4 dormancy + OncoSimulR | None / Bioconductor | Independent, most speculative |
| 6 | M5 competing risks | Longitudinal follow-up | Timing-domain L1 test |
| 7 | M6 ER activity | Tumor RNA | Functional L1 mechanism |

## Dependencies that must resolve before code runs

| # | Dependency | Status | Blocks |
|---|-----------|--------|--------|
| 1 | Remaining 9 study slots (study list unseen — abstract only) | **OPEN** (ingest appends rows as extraction lands) | M2 moderator, M1 multi-study pooling |
| 2 | Patient-level detection histories with `ctdna_fraction`, `limit_of_detection`, repeat occasions | **OPEN** | M3 identifiability |
| 3 | R stack: `metafor`, `arrow`, `dplyr`, `tidyr`, `simex`, `readr`, `spOccupancy` | install | M1(R), M2(R), M3(R) |
| 4 | Python stack: `pydantic>=2`, `pandas`, `pyarrow`, `numpy`, `pymc`, `arviz`, `gillespy2`, `pydagitty==0.1.x` | `pip install -e .` | schema/M1(PyMC)/M3(PyMC)/M4/M2-DAG |
| 5 | Prior endocrine-therapy / subtype strata for L2, L4 | **OPEN** (supplement may lack) | L2/L4 falsification |
| 6 | Assay sensitivity validation priors for p | **SHELF workshop** (`m3_shelf_priors`; placeholders until run) | M3 single-survey runs |
| 7 | R extensions: `mecor` 1.0.0, `baggr` 0.8.2 + `rstan` ≥ 2.26 toolchain, `SHELF` 1.13.0 | install (baggr = minutes of HMC) | M1-moderator, M2-pooling, M3-priors |
| 8 | Replicate purity/fraction measures per study (`mod_rep1/2`) | **OPEN** (template via `m1_mecor`) | mecor calibration |
| 9 | Longitudinal time-to-first-progression by site + ESR1 | **OPEN** (synthetic placeholder) | M5 competing risks |
| 10 | Tumor RNA (bulk or sc) + vendored PROGENy weights | **OPEN** (200-gene synthetic placeholder) | M6 ER activity |
| 11 | Bioconductor `OncoSimulR` 4.14.0 | install (no new data needed) | M4 selection challenger |
| 12 | Python `comprisk`, `decoupler` + `anndata` | install (both have tested hand-rolled fallbacks) | M5 inference, M6 cross-check |

Until #1–#2 land: `python data/make_synthetic.py` writes schema-valid
placeholders so `snakemake --cores 1` and `pytest` run end to end.
Numbers from synthetic data are plumbing checks, not findings.

## Quickstart

```bash
uv sync --extra all                    # or: pip install -e ".[all]"
uv run python data/ingest_ecosystem.py # real study row + meta anchor
uv run python data/make_synthetic.py   # synthetic tables (not study_table)
uv run python data/validate.py
uv run pytest -q
ruff check . && ruff format --check .  # lint gate (also pre-commit)
snakemake --cores 1          # needs R stack for R rules; Python rules run regardless
uv run python modules/m4_dormancy/simulate.py --n-trajectories 10 --t-end 50
marimo edit notebooks/explore_results.py  # needs --extra explore
```

## R toolchain (measured 2026-10-10, this machine)

- R 4.6.1 (`C:\Program Files\R\R-4.6.1`, pre-existing via winget, not on PATH).
- User library `C:\Users\User\R\win-library\4.6` — every `Rscript` call needs
  `$env:R_LIBS` pointed there (system library is not writable); Snakefile
  R rules assume `Rscript` on PATH with that lib.
- Versions: metafor 5.2.1, arrow 25.0.1, dplyr 1.2.1, tidyr 1.3.2,
  readr 2.2.0, jsonlite 2.0.0, simex 1.8, mecor 1.0.0, SHELF 1.13.0,
  spOccupancy 0.8.1, OncoSimulR 4.14.0, baggr 0.8.2, rstan 2.32.7.
- All R scripts execute: simex dual grid, pooled-only meta-regression,
  mecor template mode, SHELF Beta fits, occupancy PGOcc (ψ̂≈0.83,
  p̄≈0.22 on synthetic single-survey data), oncosimulr both regimes
  (baseline entropy 0.90/ESR1 10% vs therapy 0.61/85%, N=20, toy params).
  baggr parks pooled-only until study 2 (Stan path untested — no HMC run yet).

## Python tooling

| Tool | What | Where |
|------|------|-------|
| `uv` | installs + runs (`.python-version` pins 3.12) | quickstart |
| `ruff` | lint + format gate (`pyproject` config, pre-commit) | `ruff check .` |
| `hypothesis` | property tests: schema invariants hold for all inputs | `tests/test_schema_properties.py` (lean profile, 20 examples) |
| `nutpie` | Rust NUTS for the PyMC models (`--sampler nutpie`, `[fast]` extra) | `calibration_pymc.py`, `occupancy_pymc.py` |
| `rich` | progress bar + results table in M4 grid runs | `simulate.py` (`--quiet` disables) |
| `marimo` | reactive results notebook (pure `.py`, git-friendly) | `notebooks/explore_results.py` |

## Verification checklist

- Dual anchors: pooled measured OR 0.263 (Fisher p ≈ 7e-29) drives the
  `giannoudis2024_pooled` grid (λ = 0.6 → ≈ 0.11; λ = 0.8 → ≈ 0.19);
  external 0.35 (0.24–0.52, 10 studies, abstract-only) drives the
  `annals_106P_meta` grid (λ = 0.6 → ≈ 0.17; λ = 0.8 → ≈ 0.27).
  Correction: earlier drafts said "14 studies" — ledger says 10.
- `pydagitty` 0.1.x: DAG/ADMG only, PAG experimental.
- `spOccupancy::PGOcc`: priors default mean 0, var 2.72; separate
  occupancy/detection priors.
- Tellurium: `r.integrator = 'gillespie'`; GillesPy2: SSA/tau-leap/ODE, GPL-3.
- `metafor::rma(mods = ~ factor(detection_method))`, Q_M omnibus.
- `simex`: logistic / proportional-odds logistic, GPL-2+.
- `mecor` 1.0.0 (GPL-3): regression calibration / efficient RC / ML for
  covariate error, MoM for outcome error, in linear models — applied here
  to the study-level `log_or ~ moderator` regression, not the pooled OR.
- `baggr` 0.8.2 (GPL-3, Stan): partial-pooling meta-analysis +
  meta-regression + `baggr_compare`; binary-data model avoids the
  continuity correction.
- `SHELF` 1.13.0 (GPL-2|GPL-3): scriptable `fitdist` Beta elicitation;
  `elicitHeterogen` for M2's tau prior (interactive).
- Python: `ruff` gate clean; `hypothesis` invariants (ctDNA rule,
  count bounds, attenuation monotonicity) in lean profile; `nutpie`
  via `--sampler` (continuous models only — holds for M1/M3);
  `marimo` notebook is `py_compile`-clean (runtime needs `[explore]`).
- DuckDB + Pydantic is the Data Validation Engine pattern (Pydantic for
  schema here; DuckDB optional later).
