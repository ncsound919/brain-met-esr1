# Bayesian meta-analysis with partial pooling (baggr/Stan).
#
# Why: REML moderation needs a cross-study table; today the table holds
# Giannoudis2024 alone (single method), so this rule parks pooled-only
# until study 2 lands. Partial pooling then shrinks extreme small-study
# log-ORs toward the hypermean, priors on heterogeneity are explicit, and
# baggr_compare + LOO-CV adjudicate none/partial/full pooling instead of
# asserting random effects.
#
# Input:  data/raw/study_table.parquet
# Output: results/m2/baggr_summary.txt
#
# Cost: rstan C++ toolchain + minutes of HMC. For smoke runs, set
# BAGGR_ITER=500 (illustration only, not for inference).
#
# Upgrade path (no continuity correction): baggr's binary-data model
# ("Meta-analysis of binary data" vignette) fits the raw 2x2 counts with a
# binomial likelihood instead of the normal approx on log_or. Kept as a
# commented block below until the toolchain is confirmed on this machine.
#
# Ref: baggr v0.8.2 (Wiecek/Meager, GPL-3); needs rstan >= 2.26.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(baggr)
})

study_table <- read_parquet("data/raw/study_table.parquet") |>
  mutate(
    a = esr1_brain_mut + 0.5,
    b = n_brain - esr1_brain_mut + 0.5,
    c = esr1_nonCNS_mut + 0.5,
    d = n_nonCNS - esr1_nonCNS_mut + 0.5,
    tau = log((a * d) / (b * c)),
    se = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  )

# baggr rubin model expects columns tau (effect), se, group.
ma_data <- study_table |>
  transmute(tau = tau, se = se, group = study_id,
            method = factor(detection_method))

dir.create("results/m2", showWarnings = FALSE, recursive = TRUE)
if (nrow(ma_data) < 2) {
  # Partial pooling over a single group is undefined. Pooled-only M1 covers
  # today; this rule goes live when study 2 lands.
  writeLines(
    "Single-study table: baggr skipped (partial pooling needs >= 2 groups).",
    "results/m2/baggr_summary.txt"
  )
  quit(save = "no")
}

iter <- as.integer(Sys.getenv("BAGGR_ITER", unset = "2000"))

fit_partial <- baggr(ma_data, pooling = "partial", iter = iter)
fit_full <- baggr(ma_data, pooling = "full", iter = iter)
fit_metareg <- baggr(ma_data, pooling = "partial",
                     covariates = "method", iter = iter)

dir.create("results/m2", showWarnings = FALSE, recursive = TRUE)
sink("results/m2/baggr_summary.txt")
cat("=== partial pooling (random-effects analog) ===\n")
print(fit_partial)
cat("\n=== pooling comparison (none vs partial vs full) ===\n")
print(baggr_compare("partial" = fit_partial, "full" = fit_full))
cat("\n=== meta-regression on detection method ===\n")
print(fit_metareg)
if (identical(Sys.getenv("BAGGR_LOO", unset = "0"), "1")) {
  cat("\n=== leave-one-study-out CV (K=14 Stan refits; slow) ===\n")
  print(loocv(fit_metareg))
} else {
  cat("\n(LOO-CV skipped; set BAGGR_LOO=1 to run K=14 refits.)\n")
}
sink()

# --- Binary-data upgrade (raw counts, binomial likelihood) ---
# See vignette("Meta-analysis of binary data", package = "baggr").
# df_bin <- study_table |>
#   transmute(group = study_id,
#             a = esr1_brain_mut, n1 = n_brain,       # brain events/trials
#             c = esr1_nonCNS_mut, n2 = n_nonCNS)     # non-CNS events/trials
# fit_bin <- baggr(prepare_ma(df_bin, effect = "logOR"), model = "logit",
#                  pooling = "partial", iter = iter)
