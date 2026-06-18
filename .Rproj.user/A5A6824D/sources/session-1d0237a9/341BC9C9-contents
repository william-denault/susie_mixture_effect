
#set.seed(1)
n  <- 400

## A genuine per-individual mixture at ONE SNP: a fraction pi_R of individuals
## respond recessively, the rest additively.
g  <- rbinom(n, 2, 0.45); while (length(unique(g)) < 3) g <- rbinom(n, 2, 0.45)
rc <- recode_snp(g)
pi_R <- 0.5
z  <- runif(n) < pi_R                                   # TRUE = recessive responder
y2 <- ifelse(z, -0.2 * rc$recessive, -0.6 * rc$additive) + rnorm(n, sd = .25)

## --- fit: corrected EM (NOTE: no centring; intercept is estimated) ----------
fit <- ser_pim_em(
  codings = list(additive = rc$additive,
                 dominant  = rc$dominant,
                 recessive = rc$recessive),
  y = y2, residual_variance = .25^2,
  n_start = 10 )


cat("\n--- estimated mixing weights pi (started at 1/3 each) ---\n")
print(round(fit$pi, 3))
fit$mu

cat(sprintf("true recessive-responder fraction = %.3f\n", mean(z)))
cat(sprintf("intercept alpha = %.3f   ELBO = %.2f\n", fit$intercept, fit$elbo))
cat("\n--- per-coding effect posterior mean mu (truth: add=0.6, rec=0.8) ---\n")
print(round(fit$mu, 3))

## responsibilities: P(individual i belongs to each coding)
cat("\n--- mean responsibility by genotype class ---\n")
for (gg in 0:2)
  cat(sprintf("  g=%d (n=%3d):  %s\n", gg, sum(g == gg),
              paste(sprintf("%s=%.2f", fit$coding_names,
                            colMeans(fit$r[g == gg, , drop = FALSE])),
                    collapse = "  ")))

## classification quality on the IDENTIFIABLE contrast (g>=1): recessive vs
## additive responder.  (At g=0 every coding is uninformative by construction.)
m   <- g >= 1
sc  <- fit$r[m, "recessive"]
lab <- as.integer(z[m])
auc <- {
  o <- order(sc); rk <- numeric(length(sc)); rk[o] <- seq_along(sc)
  n1 <- sum(lab == 1); n0 <- sum(lab == 0)
  (sum(rk[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cat(sprintf("\nAUC(recessive responsibility | g>=1) = %.3f\n", auc))

## ELBO must be monotone non-decreasing
cat(sprintf("ELBO monotone non-decreasing: %s (min step %.2e)\n",
            all(diff(fit$elbo_hist) >= -1e-7), min(diff(fit$elbo_hist))))

## --- optional: visualise responsibilities -----------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df <- data.frame(g = g, y = y2, prob_rec = fit$r[, "recessive"],
                   prob_add = fit$r[, "additive"])
  print(
    ggplot(df, aes(g, y, color = prob_rec)) +
      geom_jitter(width = .08, height = 0) +
      scale_colour_gradientn(
        colours = c("blue", "dodgerblue", "white", "orange", "red"),
        rescaler = ~ scales::rescale_mid(.x, mid = 0.5)) +
      labs(title = "P(recessive responder) per individual",
           subtitle = sprintf("pi_rec=%.2f (truth %.2f), AUC=%.2f",
                              fit$pi["recessive"], mean(z), auc)))
}

print(round(fit$pi, 3))
