# Meta-regression with detection method as moderator.
# metafor::rma() with mods = ~ factor(detection_method).
# Input:  data/raw/study_table.parquet
# Output: results/m2/meta_regression.txt

suppressPackageStartupMessages({
  library(metafor)
  library(arrow)
  library(dplyr)
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

if (nrow(study_table) >= 3 && n_distinct(study_table$detection_method) >= 2) {
  res <- rma(
    yi = log_or,
    sei = se_log_or,
    mods = ~ factor(detection_method),
    data = study_table,
    method = "REML"
  )
  mod_note <- "moderator fitted"
} else {
  # Single-method/single-study table (today: Giannoudis2024 alone):
  # moderator unidentifiable. Pooled-only; L3 awaits the cross-study table.
  message("Moderator unidentifiable: pooled-only fit.")
  res <- rma(yi = log_or, sei = se_log_or, data = study_table, method = "REML")
  mod_note <- "POOLED ONLY -- moderator unidentifiable (single study/method)"
}

res_no_mod <- rma(yi = log_or, sei = se_log_or, data = study_table, method = "REML")

dir.create("results/m2", showWarnings = FALSE, recursive = TRUE)
sink("results/m2/meta_regression.txt")
cat(mod_note, "\n")
print(summary(res))
cat("I2 without moderator:", res_no_mod$I2, "\n")
cat("I2 with moderator:   ", res$I2, "\n")
cat("R2 explained by method:",
    (res_no_mod$tau2 - res$tau2) / res_no_mod$tau2, "\n")
sink()

# Omnibus test of moderators: Q_M statistic and p-value are in the summary.
# Significant Q_M -> detection method explains heterogeneity.
# Non-significant -> moderator not supported; the gap persists within method.
