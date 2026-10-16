# Elicit Beta priors for compartment detection probabilities via SHELF.
#
# Resolves the OPEN dependency "assay sensitivity validation priors for p":
# instead of asserting Beta(7,3), fit Beta(a,b) to expert probability
# judgements (pathologist / assay scientist) with SHELF::fitdist, which is
# non-interactive and scriptable (Shiny roulette apps are for workshops).
#
# WORKSHOP PROTOCOL (before trusting any output of this script):
#  1. Replace PLACEHOLDER judgements below with elicited 5/50/95 percentiles
#     per compartment from >=1 assay expert (SHELF protocol: individual ->
#     discussion -> RIO consensus; see shelf.sites.sheffield.ac.uk).
#  2. Record expert ids + date in elicitation_log.csv (same folder).
#  3. Re-run; downstream occupancy_pymc.py --beta-priors consumes the CSV.
#
# Also see SHELF::elicitHeterogen for a genuine prior on M2's between-study
# heterogeneity tau (interactive Shiny; run in the same workshop).
#
# Input:  none (expert judgements live in-script until the workshop)
# Output: results/m3/elicited_beta_priors.csv (compartment,alpha,beta,source)
#
# Ref: SHELF v1.13.0 (Oakley, GPL-2|GPL-3).

suppressPackageStartupMessages({
  library(SHELF)
})

# PLACEHOLDER judgements: P(p <= lo)=0.05, median, P(p <= hi)=0.95.
# p_brain lower (biopsy/ctDNA yield worry); p_nonCNS higher. REPLACE.
judgements <- list(
  brain   = list(vals = c(0.30, 0.55, 0.80), source = "PLACEHOLDER"),
  liver   = list(vals = c(0.50, 0.70, 0.85), source = "PLACEHOLDER"),
  bone    = list(vals = c(0.50, 0.70, 0.85), source = "PLACEHOLDER"),
  lung    = list(vals = c(0.50, 0.70, 0.85), source = "PLACEHOLDER"),
  lymph_node = list(vals = c(0.50, 0.70, 0.85), source = "PLACEHOLDER"),
  other   = list(vals = c(0.40, 0.60, 0.80), source = "PLACEHOLDER")
)
probs <- c(0.05, 0.50, 0.95)

rows <- lapply(names(judgements), function(comp) {
  fit <- fitdist(vals = judgements[[comp]]$vals, probs = probs,
                 lower = 0, upper = 1)
  # fit$Beta holds shape1/shape2 per fitted expert row; single judgement
  # set -> first row.
  data.frame(
    compartment = comp,
    alpha = fit$Beta[1, "shape1"],
    beta = fit$Beta[1, "shape2"],
    source = judgements[[comp]]$source,
    stringsAsFactors = FALSE
  )
})
priors <- do.call(rbind, rows)
print(priors)

dir.create("results/m3", showWarnings = FALSE, recursive = TRUE)
write.csv(priors, "results/m3/elicited_beta_priors.csv", row.names = FALSE)

if (any(priors$source == "PLACEHOLDER")) {
  message("WARNING: priors contain PLACEHOLDER judgements -- run the workshop.")
}
