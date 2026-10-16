# mecor regression calibration for an error-prone continuous moderator.
#
# Scope (honest): mecor corrects measurement error in a continuous covariate
# (or outcome) in LINEAR models with continuous outcomes. Here the linear
# model is the study-level meta-regression log_or ~ moderator, where the
# moderator (e.g. study-mean tumor purity / ctDNA fraction) is error-prone.
# It does NOT replace the M1 SIMEX grid on the pooled OR; it complements it
# by correcting the moderator slope, with bootstrap CIs and support for
# differential error (brain purity > plasma purity -- our key worry).
#
# Input:  data/raw/study_table.parquet, optionally with columns
#         mod_value, mod_rep1, mod_rep2 (replicate measures of the moderator)
# Output: results/m1/mecor_report.txt (+ mecor_template.csv if columns absent)
#
# Ref: Nab et al. 2021, mecor v1.0.0 (GPL-3).

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
})

study_table <- read_parquet("data/raw/study_table.parquet") |>
  mutate(
    a = esr1_brain_mut + 0.5,
    b = n_brain - esr1_brain_mut + 0.5,
    c = esr1_nonCNS_mut + 0.5,
    d = n_nonCNS - esr1_nonCNS_mut + 0.5,
    log_or = log((a * d) / (b * c)),
    se_log_or = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  )

dir.create("results/m1", showWarnings = FALSE, recursive = TRUE)
report <- "results/m1/mecor_report.txt"

has_replicates <- all(c("mod_value", "mod_rep1", "mod_rep2") %in% names(study_table))

# The replicate encoding below follows the MeasError(substitute, reference,
# replicate, differential) formals but is UNVERIFIED against an installed
# mecor (no R here): confirm with ?MeasError on first run. Any failure
# degrades to template mode -- the pipeline rule must still write `report`.
fit <- NULL
if (has_replicates) {
  library(mecor)
  # Regression calibration with internal replicates; B=999 bootstrap CIs.
  fit <- tryCatch(
    mecor(
      formula = log_or ~ MeasError(
        substitute = mod_value, replicate = cbind(mod_rep1, mod_rep2)
      ) + factor(detection_method),
      data = study_table,
      method = "standard",
      B = 999
    ),
    error = function(e) {
      message("mecor fit failed (", conditionMessage(e), "); falling back to template.")
      NULL
    }
  )
}

if (!is.null(fit)) {
  sink(report)
  print(summary(fit))
  cat("\nNote: if error is differential by site (brain purity > plasma),\n")
  cat("prefer mecor's differential-error methods over standard calibration.\n")
  sink()
} else {
  # No replicates yet: write the collection template + a reliability
  # sensitivity grid for the moderator slope (attenuation identity).
  write_csv(
    tibble(study_id = study_table$study_id,
           mod_value = NA_real_, mod_rep1 = NA_real_, mod_rep2 = NA_real_),
    "results/m1/mecor_template.csv"
  )
  # Naive slope of log_or on nothing-observable: show how a hypothetical
  # moderator coefficient beta_hat attenuates under reliability lambda.
  lambda_grid <- c(0.6, 0.7, 0.8, 0.9)
  sink(report)
  if (has_replicates) {
    cat("Replicate columns present but the mecor fit failed; template written to\n")
  } else {
    cat("No mod_value/mod_rep1/mod_rep2 columns: template written to\n")
  }
  cat("results/m1/mecor_template.csv. Collect replicate purity/fraction\n")
  cat("measures per study, then re-run for regression calibration.\n\n")
  cat("Sensitivity: corrected slope = naive slope / lambda:\n")
  print(tibble(lambda = lambda_grid, correction_factor = 1 / lambda_grid))
  sink()
}
