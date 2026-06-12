library(susieR)
library(ggplot2)
library(gridExtra)

# ─────────────────────────────────────────────────────────────────────────────
# 1. DATA SETUP
# ─────────────────────────────────────────────────────────────────────────────
data(N3finemapping)
X_raw <- N3finemapping$X
n     <- nrow(X_raw)
p     <- ncol(X_raw)

# Recode each SNP into additive (0/1/2), dominant (0/1/1), recessive (0/0/1)
# Most-frequent genotype coded as 0
recode_snp <- function(x) {
  sorted_vals <- names(sort(table(x), decreasing = TRUE))
  ad          <- match(x, sorted_vals) - 1
  list(ad  = ad,
       dom = as.integer(ad >= 1),
       rec = as.integer(ad == 2))
}

AD <- DOM <- REC <- matrix(0, n, p)
for (j in 1:p) {
  rc       <- recode_snp(X_raw[, j])
  AD[, j]  <- rc$ad
  DOM[, j] <- rc$dom
  REC[, j] <- rc$rec
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. SNP SUITABILITY FILTER
#    Mixture model requires all three genotype classes to be populated.
#    Low-MAF SNPs have AD ≈ DOM, making components unidentifiable.
# ─────────────────────────────────────────────────────────────────────────────
is_suitable_snp <- function(j, min_per_class = 20)
  length(table(AD[, j])) == 3 && min(table(AD[, j])) >= min_per_class

suitable_snps <- which(sapply(1:p, is_suitable_snp))
cat("SNPs suitable for mixture testing:", length(suitable_snps), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. CORE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

# ── Profile log-likelihoods (for BIC-based model comparison) ─────────────────
profile_loglik_null <- function(y) {
  y <- y - mean(y); n <- length(y)
  -n/2 * log(2*pi) - n/2 * log(sum(y^2)/n) - n/2
}

profile_loglik_single <- function(x, y) {
  x <- x - mean(x); y <- y - mean(y)
  if (var(x) == 0) return(-Inf)
  n    <- length(x)
  beta <- sum(x * y) / sum(x^2)
  rss  <- sum((y - beta * x)^2)
  -n/2 * log(2*pi) - n/2 * log(rss/n) - n/2
}

# ── 3-component mixture EM ────────────────────────────────────────────────────
#
#   z_i ~ Categorical(pi_1, pi_2, pi_3)
#   Component 1 (additive):  y_i = beta_1 * x_i^AD  + eps_i
#   Component 2 (dominant):  y_i = beta_2 * x_i^DOM + eps_i
#   Component 3 (recessive): y_i = beta_3 * x_i^REC + eps_i
#   Shared sigma^2 across components.
#
#   Separate beta per component — identifiable when betas differ.
#   Returns the observed-data log-likelihood used as profile loglik for BIC.

run_mixture3_em <- function(ad_j, dom_j, rec_j, y,
                            q_init   = NULL,
                            max_iter = 300,
                            tol      = 1e-8) {
  y_c   <- y     - mean(y)
  ad_c  <- ad_j  - mean(ad_j)
  dom_c <- dom_j - mean(dom_j)
  rec_c <- rec_j - mean(rec_j)

  if (var(ad_c) == 0 | var(dom_c) == 0) return(NULL)

  n <- length(y_c)

  # initialise responsibilities (n x 3)
  Q <- if (is.null(q_init)) {
    matrix(1/3, n, 3)
  } else {
    pmax(q_init, 1e-8)
  }
  Q <- Q / rowSums(Q)

  pi <- colMeans(Q)
  b1 <- sum(Q[,1] * ad_c  * y_c) / (sum(Q[,1] * ad_c^2)  + 1e-8)
  b2 <- sum(Q[,2] * dom_c * y_c) / (sum(Q[,2] * dom_c^2) + 1e-8)
  b3 <- sum(Q[,3] * rec_c * y_c) / (sum(Q[,3] * rec_c^2) + 1e-8)
  s2 <- var(y_c)
  ll_old <- -Inf

  for (iter in 1:max_iter) {

    # ── E-step ──────────────────────────────────────────────────────────────
    lp1 <- log(pi[1]) - (y_c - b1 * ad_c)^2  / (2 * s2)
    lp2 <- log(pi[2]) - (y_c - b2 * dom_c)^2 / (2 * s2)
    lp3 <- log(pi[3]) - (y_c - b3 * rec_c)^2 / (2 * s2)

    mx  <- pmax(lp1, lp2, lp3)
    lse <- mx + log(exp(lp1-mx) + exp(lp2-mx) + exp(lp3-mx))
    Q   <- cbind(exp(lp1-lse), exp(lp2-lse), exp(lp3-lse))

    # ── M-step ──────────────────────────────────────────────────────────────
    pi <- pmax(colMeans(Q), 1e-8); pi <- pi / sum(pi)
    b1 <- sum(Q[,1] * ad_c  * y_c) / (sum(Q[,1] * ad_c^2)  + 1e-8)
    b2 <- sum(Q[,2] * dom_c * y_c) / (sum(Q[,2] * dom_c^2) + 1e-8)
    b3 <- sum(Q[,3] * rec_c * y_c) / (sum(Q[,3] * rec_c^2) + 1e-8)
    s2 <- max(mean(
      Q[,1] * (y_c - b1*ad_c)^2 +
        Q[,2] * (y_c - b2*dom_c)^2 +
        Q[,3] * (y_c - b3*rec_c)^2
    ), 1e-8)

    ll <- sum(lse)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }

  list(Q = Q, pi = pi, beta = c(b1, b2, b3), sigma2 = s2, ll = ll)
}

# ── Multi-start EM ────────────────────────────────────────────────────────────
run_mixture3_multistart <- function(ad_j, dom_j, rec_j, y,
                                    n_starts = 20, ...) {
  best_ll <- -Inf; best_res <- NULL
  for (s in 1:n_starts) {
    r   <- matrix(runif(length(y) * 3), length(y), 3)
    q0  <- r / rowSums(r)
    res <- tryCatch(
      run_mixture3_em(ad_j, dom_j, rec_j, y, q_init = q0, ...),
      error = function(e) NULL
    )
    if (!is.null(res) && res$ll > best_ll) {
      best_ll <- res$ll; best_res <- res
    }
  }
  best_res
}

# ── Model assessment for one SNP given a (partial) residual y ─────────────────
#
#   BIC parameter counts (vs null with sigma^2 only):
#     Pure models:  1 extra param (beta)       → penalty = log(n)/2
#     Mixture K=3:  5 extra params             → penalty = 5*log(n)/2
#       (beta_1, beta_2, beta_3, pi_1, pi_2; pi_3 = 1 - pi_1 - pi_2)
#
#   IMPORTANT: y should be the partial residual from SuSiE
#   (i.e. y minus all other single effects), not raw phenotype.
#   Using raw y conflates multiple causal SNPs and the EM will partition
#   individuals based on all effects rather than this SNP alone.

assess_genetic_model <- function(snp_idx, AD, DOM, REC, y,
                                 n_starts  = 20,
                                 min_group = 10) {

  if (!(snp_idx %in% suitable_snps))
    warning(sprintf(
      "SNP %d may have low MAF — mixture components may not be identifiable",
      snp_idx))

  n_obs  <- length(y)
  L_null <- profile_loglik_null(y)
  pen1   <- log(n_obs) / 2
  pen5   <- 5 * log(n_obs) / 2

  # pure model BIC scores
  score_ad  <- profile_loglik_single(AD[,  snp_idx], y) - L_null - pen1
  score_dom <- profile_loglik_single(DOM[, snp_idx], y) - L_null - pen1
  score_rec <- profile_loglik_single(REC[, snp_idx], y) - L_null - pen1

  # mixture BIC score
  em <- run_mixture3_multistart(
    AD[, snp_idx], DOM[, snp_idx], REC[, snp_idx], y,
    n_starts = n_starts
  )

  if (is.null(em)) {
    score_mix <- -Inf; pi_est <- c(NA, NA, NA)
  } else {
    group_sizes <- colSums(em$Q > 1/3)
    score_mix   <- if (any(group_sizes < min_group)) -Inf else
      em$ll - L_null - pen5
    pi_est <- em$pi
  }

  log_bfs <- c(Additive  = score_ad,
               Dominant  = score_dom,
               Recessive = score_rec,
               Mixture   = score_mix)

  s          <- log_bfs - max(log_bfs[is.finite(log_bfs)])
  post_probs <- exp(s) / sum(exp(s))

  data.frame(
    SNP          = snp_idx,
    Model        = names(log_bfs),
    log_BF       = round(log_bfs, 2),
    Post_Prob    = round(post_probs, 4),
    pi_additive  = round(pi_est[1], 3),
    pi_dominant  = round(pi_est[2], 3),
    pi_recessive = round(pi_est[3], 3),
    Is_Mixture   = which.max(post_probs) == 4,
    row.names    = NULL
  )
}

# ── Partial residual from SuSiE fit ──────────────────────────────────────────
#
#   For effect l, the partial residual is:
#     r_l = y - sum_{l' != l} X * E[b_{l'}]
#         = y - X * (sum_all E[b_l']) + X * E[b_l]
#
#   This isolates the signal from SNP l, removing contributions from
#   all other single effects — exactly what SuSiE's IBSS computes internally.

get_partial_residual <- function(res_susie, X_stacked, y, l_idx) {
  # posterior mean for each single effect l: alpha_l * mu_l (p-vector)
  fitted_all  <- X_stacked %*% colSums(res_susie$alpha * res_susie$mu)
  fitted_this <- X_stacked %*% (res_susie$alpha[l_idx, ] * res_susie$mu[l_idx, ])
  as.numeric(y - fitted_all + fitted_this)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. SIMULATE PHENOTYPE WITH KNOWN ARCHITECTURE
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)

causal_snps <- suitable_snps[c(5, 15, 25, 35)]
true_models <- c("additive", "dominant", "recessive", "mixture")

cat("Causal SNPs:", causal_snps, "\n")
cat("Genotype counts for each causal SNP:\n")
for (k in seq_along(causal_snps))
  cat(sprintf("  SNP %d (%s): ", causal_snps[k], true_models[k]),
      table(AD[, causal_snps[k]]), "\n")

y <- rnorm(n, sd = 0.5)
y <- y + 2.0 * scale(AD[,  causal_snps[1]])   # purely additive
y <- y + 2.0 * scale(DOM[, causal_snps[2]])   # purely dominant
y <- y + 2.0 * scale(REC[, causal_snps[3]])   # purely recessive

# mixture SNP: pi = (0.5 additive, 0.3 dominant, 0.2 recessive)
# same-sign effects of different magnitude — biologically realistic
z_mix <- sample(1:3, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
y <- y +
  (z_mix == 1) * 1.0 * AD[,  causal_snps[4]] +
  (z_mix == 2) * 3.0 * DOM[, causal_snps[4]] +
  (z_mix == 3) * 2.0 * REC[, causal_snps[4]]
y <- as.numeric(scale(y))

# ─────────────────────────────────────────────────────────────────────────────
# 5. FINE-MAP WITH SUSIE ON STACKED MATRIX
# ─────────────────────────────────────────────────────────────────────────────
X_stacked <- cbind(AD, DOM, REC)
colnames(X_stacked) <- c(paste0("AD_",  1:p),
                         paste0("DOM_", 1:p),
                         paste0("REC_", 1:p))

res_susie <- susie(X_stacked, y, L = 10, verbose = FALSE)
cat("\nCredible sets found:\n")
print(res_susie$sets)

# ─────────────────────────────────────────────────────────────────────────────
# 6. MODEL ASSESSMENT USING PARTIAL RESIDUALS
#
#   For each CS from SuSiE:
#     1. Compute the partial residual (removes all other effects)
#     2. Identify the top SNP in the CS
#     3. Run model assessment on the partial residual
#
#   This ensures the EM sees only the signal from that SNP,
#   not confounded variation from other causal SNPs.
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Model assessment per credible set ===\n")

cs_results <- do.call(rbind, lapply(seq_along(res_susie$sets$cs), function(k) {

  l_idx   <- res_susie$sets$cs_index[k]      # which L in the SuSiE model
  cs_cols <- res_susie$sets$cs[[k]]           # columns in X_stacked

  # top SNP = column with highest PIP within this CS
  pips_cs  <- res_susie$pip[cs_cols]
  top_col  <- cs_cols[which.max(pips_cs)]

  # decode stacked column back to original SNP index and encoding
  if      (top_col <= p)      { snp_idx <- top_col;       encoding <- "AD"  }
  else if (top_col <= 2*p)    { snp_idx <- top_col - p;   encoding <- "DOM" }
  else                        { snp_idx <- top_col - 2*p; encoding <- "REC" }

  # partial residual isolates this SNP's contribution
  y_resid <- get_partial_residual(res_susie, X_stacked, y, l_idx)

  df              <- assess_genetic_model(snp_idx, AD, DOM, REC, y_resid,
                                          n_starts = 20)
  df$CS           <- k
  df$L_idx        <- l_idx
  df$Top_col      <- top_col
  df$SNP_encoding <- encoding
  df$Coverage     <- round(res_susie$sets$coverage[k], 3)
  df
}))

print(cs_results)

# ─────────────────────────────────────────────────────────────────────────────
# 7. SIMULATION BENCHMARK
#    For each of the 4 architectures, repeat 100 times using partial residuals
#    from a new SuSiE fit each time.
# ─────────────────────────────────────────────────────────────────────────────
run_benchmark <- function(snp_idx, true_model,
                          other_snps,            # other causal SNPs for background
                          beta_ad  = 1.0,
                          beta_dom = 2.0,
                          beta_rec = 1.5,
                          pi_mix   = c(0.5, 0.3, 0.2),
                          n_sim    = 100,
                          n_starts = 10) {

  expected <- tools::toTitleCase(true_model)
  correct  <- logical(n_sim)

  for (s in 1:n_sim) {
    noise <- rnorm(n, sd = 0.5)

    # simulate background effects from other SNPs (additive)
    y_sim <- noise
    for (j in other_snps)
      y_sim <- y_sim + scale(AD[, j])

    # add effect of focal SNP
    focal_effect <- switch(true_model,
                           additive  = beta_ad  * scale(AD[,  snp_idx]),
                           dominant  = beta_dom * scale(DOM[, snp_idx]),
                           recessive = beta_rec * scale(REC[, snp_idx]),
                           mixture   = {
                             z <- sample(1:3, n, replace = TRUE, prob = pi_mix)
                             (z==1)*beta_ad *AD[, snp_idx] +
                               (z==2)*beta_dom*DOM[,snp_idx] +
                               (z==3)*beta_rec*REC[,snp_idx]
                           }
    )
    y_sim <- as.numeric(scale(y_sim + focal_effect))

    # fit SuSiE on stacked matrix
    res_s <- tryCatch(
      susie(X_stacked, y_sim, L = 10, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(res_s) || length(res_s$sets$cs) == 0) next

    # find the CS containing our focal SNP (in any encoding)
    focal_cols <- c(snp_idx, snp_idx + p, snp_idx + 2*p)
    l_found    <- NULL
    for (k in seq_along(res_s$sets$cs)) {
      if (any(focal_cols %in% res_s$sets$cs[[k]])) {
        l_found <- res_s$sets$cs_index[k]; break
      }
    }
    if (is.null(l_found)) next

    y_resid    <- get_partial_residual(res_s, X_stacked, y_sim, l_found)
    res_assess <- assess_genetic_model(snp_idx, AD, DOM, REC, y_resid,
                                       n_starts = n_starts)
    best       <- res_assess$Model[which.max(res_assess$Post_Prob)]
    correct[s] <- (best == expected)
  }

  data.frame(
    True_Model = true_model,
    SNP        = snp_idx,
    Accuracy   = mean(correct),
    SE         = sqrt(mean(correct) * (1 - mean(correct)) / n_sim)
  )
}

set.seed(123)
cat("\n=== Simulation benchmark (100 reps per model) ===\n")
sim_results <- do.call(rbind, lapply(seq_along(causal_snps), function(k) {
  cat("  Benchmarking:", true_models[k], "\n")
  run_benchmark(
    snp_idx    = causal_snps[k],
    true_model = true_models[k],
    other_snps = causal_snps[-k],   # background effects from other causal SNPs
    n_sim      = 100,
    n_starts   = 10
  )
}))
print(sim_results)

# ─────────────────────────────────────────────────────────────────────────────
# 8. VISUALISE
# ─────────────────────────────────────────────────────────────────────────────
model_colors <- c(Additive  = "#1A85FF",
                  Dominant  = "#D41159",
                  Recessive = "darkgreen",
                  Mixture   = "orange")

# ── Panel A: posterior model probabilities per CS ─────────────────────────────
# add true model label by matching CS SNP to known causal SNPs
cs_top_snps <- cs_results[!duplicated(cs_results$CS), c("CS", "SNP")]
cs_true <- sapply(cs_top_snps$SNP, function(s) {
  idx <- which(causal_snps == s)
  if (length(idx) == 0) "unknown" else true_models[idx]
})
cs_label_map <- setNames(
  paste0("CS ", cs_top_snps$CS, " SNP ", cs_top_snps$SNP,
         "\n(true: ", cs_true, ")"),
  cs_top_snps$CS
)
cs_results$CS_label <- cs_label_map[as.character(cs_results$CS)]

pA <- ggplot(cs_results, aes(x = Model, y = Post_Prob, fill = Model)) +
  geom_col() +
  facet_wrap(~CS_label, nrow = 1) +
  scale_fill_manual(values = model_colors) +
  labs(y     = "Posterior model probability",
       x     = NULL,
       title = "A: Model assessment per credible set (partial residuals)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.text      = element_text(size = 8),
        axis.text.x     = element_text(angle = 30, hjust = 1))

# ── Panel B: estimated mixture proportions for mixture CS ────────────────────
mix_cs   <- cs_results[cs_results$CS_label == cs_label_map[
  which(cs_top_snps$SNP %in% causal_snps[true_models == "mixture"])
], ]
mix_row  <- mix_cs[1, ]

pi_df <- data.frame(
  Component = c("Additive", "Dominant", "Recessive"),
  Estimated = c(mix_row$pi_additive, mix_row$pi_dominant, mix_row$pi_recessive),
  True      = c(0.5, 0.3, 0.2)
)
pi_long <- reshape(pi_df, varying = c("Estimated", "True"),
                   v.names = "pi", timevar = "Type",
                   times = c("Estimated", "True"), direction = "long")

pB <- ggplot(pi_long, aes(x = Component, y = pi,
                          fill = Component, alpha = Type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(Additive  = "#1A85FF",
                               Dominant  = "#D41159",
                               Recessive = "darkgreen")) +
  scale_alpha_manual(values = c(Estimated = 1.0, True = 0.4)) +
  labs(y     = "Mixture proportion (pi)",
       x     = NULL,
       alpha = NULL,
       title = "B: Estimated vs true mixture proportions (mixture SNP)") +
  theme_bw(base_size = 12)

# ── Panel C: classification accuracy ─────────────────────────────────────────
sim_results$True_Model <- factor(sim_results$True_Model, levels = true_models)

pC <- ggplot(sim_results,
             aes(x      = True_Model,
                 y      = Accuracy,
                 ymin   = pmax(0, Accuracy - 1.96 * SE),
                 ymax   = pmin(1, Accuracy + 1.96 * SE),
                 colour = True_Model)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(additive  = "#1A85FF",
                                 dominant  = "#D41159",
                                 recessive = "darkgreen",
                                 mixture   = "orange")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(y     = "Proportion correct model selected",
       x     = "True genetic model",
       title = "C: Classification accuracy (100 simulations per model)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

grid.arrange(pA, pB, pC, nrow = 3)
