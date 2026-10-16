# Module 4 — dormancy (bet-hedging)

Two-state stochastic model: dormant (D) ↔ active (A).
Pre-registered prediction: pause-targeting the active state **increases**
diversification risk.

Real pause biology (ecosystem ledger `MXD4_*`, Garcia-Gomez et al. Cancer
Cell 2026): an MXD4-dependent proliferative pause enables brain
micrometastasis survival, marked SOX2/p27/p21/BCL-2. The D-state now has
candidate markers — and a comparator result this model must address:
BCL-2/ATR/HSP90 inhibition **reduced** pause-state cells and recurrence,
i.e. killing the pause compartment helped. The bet-hedging prediction
(diversification risk rises under pause-targeting) and the MXD4 direction
are in tension; the OncoSimulR challenger plus marker-stratified data
adjudicate. First measured run (mu=1e-4, T=150, N=20, toy params):
baseline entropy 0.90 with ESR1 in 10% of tumors vs therapy entropy 0.61
with ESR1 in 85% — directional selection, not diversification. TRAP-seq GEO data (accessions unseen) is a standing M6 lead.

## Falsification condition

Flat or reversed effect of pause-targeting across the plausible parameter
range → bet-hedging hypothesis not supported by this model. Write the
prediction and this condition down **before** running `simulate.py`.

## Most speculative module

Independent of M1–M3 (no data dependency). Build last.

## Parameterization note

Default rates (proliferation 0.1 > death 0.05) give net growth of A, so
pause = 0 trajectories explode — that is a toy-parameter property, not a
finding. Calibrate rates to a plausible regime before reading direction.

## Challenger: OncoSimulR selection model

`oncosimulr_selection.R` (Bioconductor 3.23, v4.14.0, GPL-3) asks whether
selection alone reproduces diversification, without phenotypic switching.
Two regimes (ESR1 neutral vs endocrine-selected, PIK3CA epistasis),
forward-time BNB simulation, Shannon entropy of whole-tumor-sampled
genotypes as the diversification metric. Sweep (entropy collapse + ESR1
fixation) under therapy favours bet-hedging; maintained entropy under
selection makes bet-hedging unnecessary. Death-rate fitness effects and
adaptive-therapy schedules are noted extensions, not implemented.

## Dependencies

`gillespy2` (SSA / tau-leaping / ODE solvers; GPL-3, verified against 1.8.3).
Grid runs use `NumPySSASolver` (pure Python, no toolchain);
`SSACSolver` is the faster drop-in where g++ exists. `simulate.py` keeps a
pure-NumPy tau-leap fallback so smoke runs work without gillespy2 at all.
Tellurium alternative: `r.integrator = 'gillespie'` with parameters
changed between segments. `oncosimulr_selection.R` needs Bioconductor
`OncoSimulR` + Rcpp toolchain.
