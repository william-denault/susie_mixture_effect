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
# 4. SIMULATE PHENOTYPE WITH KNOWN ARCHITECTURE
# ─────────────────────────────────────────────────────────────────────────────
set.seed(1)

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
z_mix <- 0*sample(1:3, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
y <- y +
  (z_mix == 1) * 1.0 * AD[,  causal_snps[4]] +
  (z_mix == 2) * 3.0 * DOM[, causal_snps[4]] +
  (z_mix == 3) * 2.0 * REC[, causal_snps[4]]
y <- as.numeric(scale(y))

# ─────────────────────────────────────────────────────────────────────────────
#  FINE-MAP SNP that stems from different models using additive effect
# ─────────────────────────────────────────────────────────────────────────────
X_stacked <- cbind(AD, DOM, REC)
colnames(X_stacked) <- c(paste0("AD_",  1:p),
                         paste0("DOM_", 1:p),
                         paste0("REC_", 1:p))

res_susie <- susie(AD, y, L = 10, verbose = FALSE) #susie(X_stacked, y, L = 10, verbose = FALSE)
cat("\nCredible sets found:\n")
print(res_susie$sets)
causal_snps
# ─────────────────────────────────────────────────────────────────────────────
#  FINE-MAP SNP that stems from different models using Stacked matrix
# ─────────────────────────────────────────────────────────────────────────────
X_stacked <- cbind(AD, DOM, REC)
colnames(X_stacked) <- c(paste0("AD_",  1:p),
                         paste0("DOM_", 1:p),
                         paste0("REC_", 1:p))

res_susie <- susie(AD, y, L = 10, verbose = FALSE) #susie(X_stacked, y, L = 10, verbose = FALSE)
cat("\nCredible sets found:\n")
print(res_susie$sets)
causal_snps

# ─────────────────────────────────────────────────────────────────────────────
#  FINE-MAP SNP same simulation under fully additive (sanity check)
# ─────────────────────────────────────────────────────────────────────────────
set.seed(1)

causal_snps <- suitable_snps[c(5, 15, 25, 35)]
true_models <- c("additive", "dominant", "recessive", "mixture")

cat("Causal SNPs:", causal_snps, "\n")
cat("Genotype counts for each causal SNP:\n")
for (k in seq_along(causal_snps))
  cat(sprintf("  SNP %d (%s): ", causal_snps[k], true_models[k]),
      table(AD[, causal_snps[k]]), "\n")

y <- rnorm(n, sd = 0.5)
y <- y + 2.0 * scale(AD[,  causal_snps[1]])   # purely additive
y <- y + 2.0 * scale(AD[, causal_snps[2]])   # purely dominant
y <- y + 2.0 * scale(AD[, causal_snps[3]])   # purely recessive




res_susie <- susie(AD, y, L = 10, verbose = FALSE) #susie(X_stacked, y, L = 10, verbose = FALSE)
cat("\nCredible sets found:\n")
print(res_susie$sets)
causal_snps
