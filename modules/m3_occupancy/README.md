# Module 3 — occupancy

Formalizes "ESR1 not detected in brain" as **nondetection, not absence**.

Empirical precedent (ecosystem ledger `ESR1_CTDNA_BRAIN`): ESR1 Y537S found
in a resected brain metastasis but **not** in matched plasma ctDNA (Curr
Oncol 2026 case report). Single case — but exactly the `p_brain < p_plasma`
mechanism this module quantifies. Real-world survey cadence anchor:
Guardant360 CDx ESR1 surveillance every 2–3 months (ledger `S6_CTDNA_INTERVAL`;
labeling recommends q3mo) — the repeat-occasion design M3 needs already
exists in trial practice.

If `p_brain < p_nonCNS`, naive brain prevalence underestimates true
occupancy ψ. The model estimates ψ (true occupancy) and p (detection)
separately. Output of interest: ψ_brain vs ψ_nonCNS, with p_brain as the
nuisance parameter quantifying naive-estimate bias.

## Identifiability warning

Single-survey-per-patient data cannot separate ψ from p without
informative priors. Genuine occupancy needs **repeated surveys**
(serial scans / serial draws). With one row per patient×compartment,
use the PyMC marginal model (`occupancy_pymc.py`) with assay-sensitivity
priors and read the result as prior-sensitive, or collect repeat
occasions before running `occupancy.R` (`PGOcc`).

## Novelty check

No published paper found applying occupancy modeling to metastatic
organotropism with detection–nondetection data (suggestive, not
conclusive). Run a focused PubMed + bioRxiv search before claiming
novelty in a manuscript.

## Dependencies

- `occupancy.R`: `spOccupancy` (`PGOcc`; priors default mean 0, var 2.72),
  `arrow`, `dplyr`, `tidyr`.
- `occupancy_pymc.py`: `pymc`, `arviz`, `pandas`, `numpy`.
- `prior_elicitation_shelf.R`: `SHELF` v1.13.0 (GPL-2|GPL-3).

## Extension: SHELF prior elicitation

`prior_elicitation_shelf.R` replaces the asserted `Beta(7,3)` detection
priors with expert-elicited `Beta(a,b)` per compartment
(`SHELF::fitdist`, scriptable; Shiny roulette apps reserved for the
workshop). It ships with PLACEHOLDER judgements and refuses to be trusted
until the workshop replaces them — see the in-script protocol. Output
`results/m3/elicited_beta_priors.csv` feeds
`occupancy_pymc.py --beta-priors`. Same workshop can elicit M2's
heterogeneity prior via `SHELF::elicitHeterogen`.
