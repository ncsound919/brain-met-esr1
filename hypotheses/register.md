# hypotheses/register.md

Synced with the ecosystem register (`SCIENCE/Oncology Ecosystem/brain-met-esr1`,
updated 2026-10-10, Giannoudis 2024 supplement). Statuses below are **measured**
where marked; the rest remain untested.

## Anchors (both depletion; magnitude differs)

- Pooled measured: OR **0.263** (53/822 brain vs 3220/15516 non-CNS),
  Fisher p ≈ 7e-29. Single Foundation Medicine cohort, `samples` denominator.
- External meta: OR **0.35 (0.24–0.52), 10 studies**, Annals of Oncology
  abstract 106P, search to 2025-06-05. Abstract ONLY — study list, detection
  methods, heterogeneity all unseen (ledger `ESR1_META`).

## L1 — True site effect → SUPPORTED within ER+/HER2−

16.98% brain vs 32.49% other-met, p = 1e-5 (paper's own stratified table).
Limit: stratum counts unreleased, independent Fisher/CI blocked.

## L2 — Subtype composition → PARTIAL (explains magnitude, not direction)

Gap attenuates pooled OR 0.26 → stratified OR ≈ 0.42 but persists significant.
Composition explains part of the magnitude, not the effect.

## L3 — Detection-method artefact → UNTESTED (cross-study question)

Giannoudis2024 is single-method (Foundation Medicine), so L3 is tested by
the cross-study meta-regression (M2), not within this cohort.
**Would be falsified by:** meta-regression showing no moderation by method,
or a within-method comparison showing the gap persists.

## L4 — Treatment-history confounding → UNTESTED

Supplement carries no endocrine-therapy history.
**Would be falsified by:** adjustment for prior endocrine therapy
attenuating the gap toward zero.

## Resistance warning (kept, updated)

The pooled 6.4% vs 20.75% contrast was consistent with L1 and L2
simultaneously; the stratified table distinguishes them and favors L1
with L2 attenuation. Do not cite the pooled contrast alone.

## M5 — Competing-risks timing

**Would be falsified by:** Gray's test on brain-first CIF (mut vs wt)
with p ≥ 0.05 AND Fine-Gray ESR1 coefficient ≈ 0 after adjusting
for detection method → no timing difference, L1 unsupported.
Supports L3 over L1 if method-adjusted CIFs coincide.

## M6 — ER-pathway decoupling

**Would be falsified (as L1 mechanism) by:** ESR1-mut brain mets showing
Estrogen activity as high as wt → mutation does not blunt ER output.
Supports L1 mechanism if mut+brain activity is specifically low.
Needs real RNA; synthetic placeholder encodes the high-activity case
by construction and proves nothing.

## Status

| Lens | Status | Required data |
|------|--------|---------------|
| L1 | **supported** within ER+/HER2− (paper table, counts unreleased) | stratum counts for independent CI |
| L2 | **partial** (0.26 → 0.42) | — |
| L3 | untested | 9 remaining studies (list unseen) |
| L4 | untested | treatment-history variables |
| M5 | untested | longitudinal time-to-progression |
| M6 | untested | tumor RNA |
