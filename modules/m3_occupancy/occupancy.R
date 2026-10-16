# Bayesian occupancy model with compartment-specific detection.
# Compartment = site, serial scan = survey, detection = per-survey p.
#
# REQUIRES repeated occasions per patient x compartment. If the input has
# one row per patient x compartment, the model is not identified without
# informative priors -- see README.md and prefer occupancy_pymc.py there.
#
# API verified against spOccupancy::PGOcc docs (rdrr, Oct 2026):
#   PGOcc(occ.formula, det.formula, data, inits, priors, n.samples,
#         n.burn, n.thin, n.chains, ...).
# Return object carries beta/alpha/z/psi.samples (coda); per-site detection
# probabilities are NOT stored -- use fitted(). Priors default mean 0,
# variance 2.72. det.covs elements must be NUMERIC (site vector length J or
# J x occasions matrix), so compartments are integer-coded below.
#
# Input:  data/raw/detection_history.parquet
# Output: results/m3/occupancy_summary.csv

suppressPackageStartupMessages({
  library(spOccupancy)
  library(arrow)
  library(dplyr)
  library(tidyr)
})

det <- read_parquet("data/raw/detection_history.parquet")

# Wide detection-history matrix: rows = patients, cols = occasions.
# Expects an `occasion` column when repeated surveys exist; otherwise each
# patient x compartment contributes a single occasion and PGOcc priors
# carry the identification (flagged, not hidden).
if (!"occasion" %in% names(det)) {
  det$occasion <- 1L
  message("No occasion column: single-survey data; priors carry identification.")
}

wide <- det |>
  select(patient_id, compartment, occasion, detected) |>
  mutate(occasion = paste0(compartment, "_", occasion)) |>
  # id_cols MUST be patient_id only: compartment varies within patient and
  # would otherwise become an id column (180 patient x compartment rows of
  # mostly zeros) and coerce y to character. Occasion names embed the
  # compartment already ("brain_1"), so no information is lost.
  pivot_wider(
    id_cols = patient_id, names_from = occasion,
    values_from = detected, values_fill = 0
  )

y <- as.matrix(wide[, -1])
comp_levels <- sort(unique(det$compartment))
occasion_comp <- sub("_[0-9]+$", "", colnames(y))
comp_codes <- matrix(
  match(occasion_comp, comp_levels),
  nrow = nrow(y), ncol = ncol(y), byrow = TRUE
)
message("Compartment codes (1-based, from match()): ", paste0(seq_along(comp_levels),
  "=", comp_levels,
  collapse = ", "
))

occ_covs <- data.frame(intercept = rep(1, nrow(y)))
det_covs <- list(compartment = comp_codes)

# psi ~ 1 (intercept; extend with covariates), p ~ compartment (occasion effect).
# Informative priors on p should come from assay sensitivity validation
# (see prior_elicitation_shelf.R) -- pass via priors$alpha.normal.
data_list <- list(y = y, occ.covs = occ_covs, det.covs = det_covs)
inits_list <- list(z = apply(y, 1, max))
prior_list <- list(
  beta.normal = list(mean = 0, var = 2.72),
  alpha.normal = list(mean = 0, var = 2.72)
)

out <- PGOcc(
  occ.formula = ~ 1,
  det.formula = ~ compartment,
  data = data_list,
  inits = inits_list,
  priors = prior_list,
  n.samples = 5000,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 1,
  verbose = FALSE
)
print(summary(out))

psi_mean <- mean(colMeans(as.matrix(out$psi.samples)))
# fitted() returns a list with p.samples = [kept draws x sites x occasions].
p_hat <- tryCatch(
  fitted(out)$p.samples,
  error = function(e) {
    message("fitted() unavailable: ", conditionMessage(e))
    NULL
  }
)
p_mean <- if (!is.null(p_hat)) mean(p_hat, na.rm = TRUE) else NA_real_

dir.create("results/m3", showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(stat = c("psi_mean", "p_mean"), value = c(psi_mean, p_mean)),
  "results/m3/occupancy_summary.csv",
  row.names = FALSE
)
