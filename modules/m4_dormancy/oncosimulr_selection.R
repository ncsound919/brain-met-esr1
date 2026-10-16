# Genetic-selection challenger to the M4 bet-hedging model (OncoSimulR).
#
# Question: can SELECTION alone (endocrine therapy favouring ESR1-mutant
# clones) reproduce clonal-diversification patterns, or is phenotypic
# switching (D<->A bet-hedging) needed? Two regimes, same mutation process:
#   baseline: ESR1 neutral (birth rate == WT)
#   therapy:  ESR1 advantageous (endocrine selection) + PIK3CA epistasis
# Metric: Shannon entropy of sampled genotypes + ESR1-clone frequency.
# If therapy-regime entropy matches/exceeds the bet-hedging model's
# diversification, selection suffices and bet-hedging is not required.
# If therapy causes a selective sweep (entropy collapse + ESR1 fixation),
# maintained diversity in real data favours bet-hedging.
#
# API note: written against OncoSimulR 4.x (Bioconductor 3.23, GPL-3).
# Confirm arg names in the installed vignette: browseVignettes("OncoSimulR").
# Forward-time BNB algorithm (Mather et al. 2012); fitness = birth rates here
# (death-modifier extension noted in README, not implemented).
#
# Input:  none (pure simulation; seeds for replicates)
# Output: results/m4/oncosimulr_diversity.csv

suppressPackageStartupMessages({
  library(OncoSimulR)
})

N_REP <- as.integer(Sys.getenv("ONCOSIMULR_REP", unset = "20"))
FINAL_T <- 150

run_regime <- function(regime, seed_start) {
  # mu=1e-4 x initSize 500 x T=150 ~= 7.5 expected mutants/run: regimes must
  # actually generate diversity, else the comparison is vacuous (REP=2 smoke
  # at mu=1e-5/T=100 returned monoclonal runs).
  mu <- c("ESR1" = 1e-4, "PIK3CA" = 1e-4)
  if (regime == "baseline") {
    fe <- allFitnessEffects(genotFitness = data.frame(
      Genotype = c("WT", "ESR1", "PIK3CA", "ESR1, PIK3CA"),
      Birth = c(1.0, 1.0, 1.1, 1.15)
    ))
  } else { # therapy: endocrine selection favours ESR1; epistasis with PIK3CA
    fe <- allFitnessEffects(genotFitness = data.frame(
      Genotype = c("WT", "ESR1", "PIK3CA", "ESR1, PIK3CA"),
      Birth = c(1.0, 1.3, 1.1, 1.6)
    ))
  }
  out <- lapply(seq_len(N_REP), function(r) {
    sim <- oncoSimulPop(
      Nindiv = 1, fp = fe, mu = mu, initSize = 500,
      finalTime = FINAL_T, onlyCancer = FALSE, keepPhylog = TRUE,
      seed = seed_start + r, model = "Exp", mc.cores = 1
    )
    samp <- suppressWarnings(samplePop(sim, timeSample = "last", typeSample = "wholeTumor"))
    # One row per tumor (Nindiv=1): the per-tumor genotype profile.
    # Diversity is measured ACROSS replicate tumors per regime (below),
    # not within a single whole-tumor sample (which is degenerate by
    # construction). Extinct runs are flagged, not silently NA'd.
    mat <- as.matrix(samp)
    extinct <- all(is.na(mat))
    geno <- if (extinct) "EXTINCT" else paste((mat > 0) * 1, collapse = "")
    ecol <- grep("^ESR1$", colnames(mat), value = TRUE)
    esr1 <- if (extinct || length(ecol) == 0) NA_real_ else as.numeric(mat[, ecol] > 0)
    data.frame(regime = regime, replicate = r, genotype = geno,
               esr1_present = esr1, extinct = extinct)
  })
  do.call(rbind, out)
}

rep_df <- rbind(run_regime("baseline", 1000), run_regime("therapy", 2000))

# Regime-level diversification: Shannon entropy over non-extinct sampled
# genotypes + ESR1 presence fraction + extinction fraction.
summarize_regime <- function(d) {
  alive <- d[!d$extinct, ]
  freq <- table(alive$genotype) / nrow(alive)
  data.frame(
    regime = d$regime[1],
    n_tumors = nrow(d),
    extinction_frac = mean(d$extinct),
    n_distinct_genotypes = length(freq),
    shannon_entropy = -sum(freq * log(freq)),
    esr1_freq = mean(alive$esr1_present, na.rm = TRUE)
  )
}
res <- do.call(rbind, lapply(split(rep_df, rep_df$regime), summarize_regime))
print(res)

dir.create("results/m4", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/m4/oncosimulr_diversity.csv", row.names = FALSE)
write.csv(rep_df, "results/m4/oncosimulr_replicates.csv", row.names = FALSE)
