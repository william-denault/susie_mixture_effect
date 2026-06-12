# =============================================================================
#  SER-Mix simulation study (O2 of CHIMERA)
#
#  Demonstrates the core claims of the proposal on simulated data:
#   (S1) Coding recovery: SER-Mix labels a single causal variant's mode of
#        action (additive / dominant / recessive) correctly.
#   (S2) Heterogeneous regulation: when two co-localised causal variants act
#        through DIFFERENT codings, standard additive SuSiE either loses power
#        or inflates credible sets, while SuSiE-Mix recovers both signals with
#        the correct coding.
#   (S3) Null calibration: under pure additivity SER-Mix does NOT spuriously
#        prefer a non-additive coding (false-positive control).
#
#  Requires: susieR (baseline) and the MixSER package (this project).
# =============================================================================

## ---- 0. setup ---------------------------------------------------------------
suppressMessages({
  library(susieR)
  # Load MixSER from source (adjust path if installed):
  if (requireNamespace("MixSER", quietly = TRUE)) {
    library(MixSER)
  } else {
    pkg <- "C:/Document/Serieux/Travail/Package/MixSER"
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(pkg)
    } else {
      for (f in list.files(file.path(pkg, "R"), full.names = TRUE)) source(f)
    }
  }
})
set.seed(1)

## ---- helper: simulate a genotype block with realistic LD --------------------
sim_genotypes <- function(n, p, maf = 0.4, rho = 8) {
  # Two haplotypes per individual; smooth LD via an exponential-decay latent.
  Sig <- exp(-abs(outer(seq_len(p), seq_len(p), "-")) / rho)
  Lc  <- chol(Sig + 1e-6 * diag(p))
  draw_hap <- function() {
    z   <- matrix(rnorm(n * p), n, p) %*% Lc
    thr <- apply(z, 2, quantile, probs = 1 - maf)
    sweep(z, 2, thr, ">") * 1
  }
  draw_hap() + draw_hap()          # n x p dosage in {0,1,2}
}

## =============================================================================
## S1 + S3.  Single causal variant: coding recovery and null calibration
## =============================================================================
single_variant_experiment <- function(n = 600, p = 50, n_rep = 100,
                                       beta = 1.3, maf = 0.4) {
  models <- c("additive", "dominant", "recessive")
  out <- data.frame()
  for (truth in models) {
    correct <- 0L
    for (r in seq_len(n_rep)) {
      G   <- sim_genotypes(n, p, maf)
      cod <- recode_genotypes(G)
      j   <- sample(which(cod$identifiable), 1)
      xj  <- switch(truth,
                    additive  = cod$additive[, j],
                    dominant  = cod$dominant[, j],
                    recessive = cod$recessive[, j])
      y   <- beta * as.numeric(scale(xj)) + rnorm(n)

      s <- ser_mix(y, cod[c("additive", "dominant", "recessive")],
                   residual_variance = var(y))
      top  <- which.max(s$pip)
      best <- names(which.max(s$coding_post[top, ]))
      correct <- correct + (top == j && best == truth)
    }
    out <- rbind(out, data.frame(true_model = truth,
                                 coding_accuracy = correct / n_rep))
  }
  out
}

cat("\n=== S1/S3: coding recovery (single causal variant) ===\n")
s1 <- single_variant_experiment()
print(s1)

## =============================================================================
## S2.  Heterogeneous regulation: additive vs SuSiE-Mix credible sets
## =============================================================================
het_regulation_experiment <- function(n = 600, p = 50, n_rep = 100,
                                       beta_add = 0.7, beta_rec = 1.1,
                                       maf = 0.42, L = 8) {
  res <- data.frame()
  for (r in seq_len(n_rep)) {
    G   <- sim_genotypes(n, p, maf)
    cod <- recode_genotypes(G)
    ids <- which(cod$identifiable)
    ca  <- ids[round(length(ids) * 0.25)]   # additive causal
    cb  <- ids[round(length(ids) * 0.75)]   # recessive causal
    y   <- beta_add * as.numeric(scale(cod$additive[, ca])) +
           beta_rec * as.numeric(scale(cod$recessive[, cb])) + rnorm(n)

    ## baseline 1: standard additive SuSiE on additive coding
    fit_add <- susie(cod$additive, y, L = L, verbose = FALSE)
    n_cs_add <- length(fit_add$sets$cs)

    ## baseline 2: naive "stacked" additive SuSiE on [A | D | R]
    Xstk <- cbind(cod$additive, cod$dominant, cod$recessive)
    fit_stk <- susie(Xstk, y, L = L, verbose = FALSE)
    n_cs_stk <- length(fit_stk$sets$cs)

    ## SuSiE-Mix
    fit_mix <- susie_mix(G, y, L = L, codings = cod[c("additive","dominant","recessive")])
    n_cs_mix <- length(fit_mix$sets)
    tops <- vapply(fit_mix$sets, `[[`, integer(1), "top_variant")
    cods <- vapply(fit_mix$sets, `[[`, character(1), "best_coding")
    got_add <- any(tops == ca & cods == "additive")
    got_rec <- any(tops == cb & cods == "recessive")

    res <- rbind(res, data.frame(
      rep = r,
      n_cs_additive = n_cs_add,
      n_cs_stacked  = n_cs_stk,
      n_cs_mix      = n_cs_mix,
      mix_both_correct = got_add && got_rec))
  }
  res
}

cat("\n=== S2: heterogeneous regulation (truth = 2 causal variants) ===\n")
s2 <- het_regulation_experiment(n_rep = 100)
summ <- data.frame(
  method = c("additive SuSiE", "stacked SuSiE", "SuSiE-Mix"),
  mean_n_cs = c(mean(s2$n_cs_additive), mean(s2$n_cs_stacked), mean(s2$n_cs_mix)),
  pct_exactly_2 = c(mean(s2$n_cs_additive == 2),
                    mean(s2$n_cs_stacked == 2),
                    mean(s2$n_cs_mix == 2)) * 100,
  pct_inflated_gt2 = c(mean(s2$n_cs_additive > 2),
                       mean(s2$n_cs_stacked > 2),
                       mean(s2$n_cs_mix > 2)) * 100)
print(summ, row.names = FALSE)
cat(sprintf("\nSuSiE-Mix recovered BOTH causal codings in %.0f%% of reps\n",
            100 * mean(s2$mix_both_correct)))

## ---- optional plots ---------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  cs_long <- data.frame(
    n_cs = c(s2$n_cs_additive, s2$n_cs_stacked, s2$n_cs_mix),
    method = rep(c("additive SuSiE", "stacked SuSiE", "SuSiE-Mix"),
                 each = nrow(s2)))
  p1 <- ggplot(cs_long, aes(factor(n_cs), fill = method)) +
    geom_bar(position = "dodge") +
    geom_vline(xintercept = 2, linetype = "dashed") +
    labs(x = "# credible sets (truth = 2)", y = "reps",
         title = "Credible-set count under heterogeneous regulation") +
    theme_bw()
  p2 <- ggplot(s1, aes(true_model, coding_accuracy, fill = true_model)) +
    geom_col() + ylim(0, 1) +
    labs(x = "true coding", y = "SER-Mix coding accuracy",
         title = "Coding recovery (single causal variant)") +
    theme_bw() + theme(legend.position = "none")
  ggsave("sermix_sim_credible_sets.png", p1, width = 7, height = 4, dpi = 130)
  ggsave("sermix_sim_coding_recovery.png", p2, width = 6, height = 4, dpi = 130)
  cat("\nSaved plots: sermix_sim_credible_sets.png, sermix_sim_coding_recovery.png\n")
}

cat("\nDone.\n")
