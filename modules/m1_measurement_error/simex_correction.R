# Module 1: SIMEX measurement-error correction, dual-anchored.
# Input:  data/raw/study_table.parquet + data/raw/meta_anchor.json
# Output: results/m1/corrected_or.csv + console meta-analysis summary.
#
# Anchors (both depletion; magnitude differs):
#   pooled: data-driven pooled log OR (today: Giannoudis2024 single cohort,
#           measured OR 0.263, Fisher p ~ 7e-29). Recomputed every run.
#   meta:   external Annals of Oncology abstract 106P anchor, OR 0.35
#           (0.24-0.52), 10 studies, abstract ONLY (methods unseen).
# Attenuation identity per anchor: log(OR_obs) = lambda * log(OR_true).
# Sensitivity grid, NOT a single corrected estimate -- see README.md.

suppressPackageStartupMessages({
  library(metafor)
  library(arrow)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

study_table <- read_parquet("data/raw/study_table.parquet")

# Effect sizes: log OR per study (0.5 continuity correction).
study_table <- study_table |>
  mutate(
    a = esr1_brain_mut + 0.5,
    b = n_brain - esr1_brain_mut + 0.5,
    c = esr1_nonCNS_mut + 0.5,
    d = n_nonCNS - esr1_nonCNS_mut + 0.5,
    log_or = log((a * d) / (b * c)),
    se_log_or = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  )

# Pooled estimate (intercept-only; == the single-study estimate while n == 1).
pooled <- rma(yi = log_or, sei = se_log_or, data = study_table, method = "REML")
print(summary(pooled))

# Moderator model only when identifiable (else L3 awaits the cross-study table).
if (nrow(study_table) >= 3 && n_distinct(study_table$detection_method) >= 2) {
  res <- rma(
    yi = log_or, sei = se_log_or,
    mods = ~ factor(detection_method),
    data = study_table, method = "REML"
  )
  print(summary(res))
} else {
  message("Single-method/single-study table: moderator unidentifiable; pooled only.")
}

anchor <- fromJSON("data/raw/meta_anchor.json")
lambda_grid <- c(0.6, 0.7, 0.8, 0.9)

grid_for <- function(name, or_obs, lo, hi) {
  tibble(
    anchor = name,
    lambda = lambda_grid,
    or_true = exp(log(or_obs) / lambda_grid),
    or_lo = exp(log(lo) / lambda_grid),
    or_hi = exp(log(hi) / lambda_grid)
  )
}

corrected <- bind_rows(
  grid_for(
    "giannoudis2024_pooled",
    exp(coef(pooled)),
    exp(coef(pooled) - 1.96 * pooled$se),
    exp(coef(pooled) + 1.96 * pooled$se)
  ),
  grid_for(
    "annals_106P_meta",
    anchor$annals_106P$or, anchor$annals_106P$lo, anchor$annals_106P$hi
  )
)
corrected$interval_width <- corrected$or_hi - corrected$or_lo
print(corrected)

dir.create("results/m1", showWarnings = FALSE, recursive = TRUE)
write_csv(corrected, "results/m1/corrected_or.csv")
