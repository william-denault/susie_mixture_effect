# Identifiability of the per-individual mixture-effect model

**Scope.** This note makes rigorous the claim that the *per-individual* mixture
model (Eq. 1 of CHIMERA) is **near non-identifiable at a single locus**,
separates what *is* identifiable (and powers SuSiE-Mix) from what is *not*,
derives a regularization that falls back to additive, and gives a calibrated
Bayes-factor procedure to test whether a credible set is genuinely a
per-individual mixture. All scaling laws below are confirmed numerically
(`/identifiability` checks; slopes reproduce to 2 d.p.).

Throughout, work in residual-standard-deviation units (\(\sigma=1\)) and write
the **standardized effect** \(\theta=\beta/\sigma\). For eQTL/GWAS, \(\theta\)
is small: a variant explaining \(h^2\) of variance at MAF \(f\) has
\(\theta \approx \sqrt{h^2 / (2f(1-f))}\), typically \(0.05\)–\(0.3\).

---

## 1. The model and its reduction to three genotype strata

One causal SNP, genotype \(g_i\in\{0,1,2\}\), codings
\[
x^A(g)=g,\qquad x^D(g)=\mathbf 1\{g\ge1\},\qquad x^R(g)=\mathbf 1\{g=2\}.
\]
Latent responder type \(z_i\in\{A,D,R\}\) with \(P(z_i=m)=\pi_m\),
**independent of genotype**, and
\[
y_i \mid z_i=m,\,g_i \;\sim\; N\!\big(x^m(g_i)\,\beta_m,\;\sigma^2\big),
\qquad
p(y_i\mid g_i)=\sum_{m}\pi_m\,N\!\big(x^m(g_i)\beta_m,\sigma^2\big). \tag{1}
\]

Because \(g\) has only three values, **the data are fully summarized by the
three conditional laws** \(p(y\mid g=0),p(y\mid g=1),p(y\mid g=2)\). Evaluating
the codings:

| \(g\) | additive pred. | dominant pred. | recessive pred. | conditional law of \(y\mid g\) |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | \(N(0,\sigma^2)\) — **all components coincide; uninformative** |
| 1 | \(\beta_A\) | \(\beta_D\) | 0 | \(\pi_A N(\beta_A,\sigma^2)+\pi_D N(\beta_D,\sigma^2)+\pi_R N(0,\sigma^2)\) |
| 2 | \(2\beta_A\) | \(\beta_D\) | \(\beta_R\) | \(\pi_A N(2\beta_A,\sigma^2)+\pi_D N(\beta_D,\sigma^2)+\pi_R N(\beta_R,\sigma^2)\) |

Two structural facts already follow. (i) Homozygous-reference individuals
(\(g=0\), frequency \((1-f)^2\)) carry **no information** about the mixture — at
realistic MAF that is the majority of the sample. (ii) Within each informative
stratum the law is a **finite Gaussian location-mixture with common variance**.

---

## 2. What *is* identifiable: the mean structure (order \(\theta\))

Let \(\mu_g=\mathbb E[y\mid g]\). From the table,
\[
\mu_0=0,\quad \mu_1=\pi_A\beta_A+\pi_D\beta_D,\quad
\mu_2=2\pi_A\beta_A+\pi_D\beta_D+\pi_R\beta_R .
\]
The **dominance deviation** \(d=\mu_2-2\mu_1=\pi_R\beta_R-\pi_D\beta_D\) is the
only way the mean departs from additivity, and it is \(O(\theta)\). Hence:

> The presence of **non-additivity** — and which single coding (A/D/R) best
> fits a variant — is identified at the standard rate: information \(O(\theta^2)\)
> per sample, detectable with \(n\sim\theta^{-2}\).

**This is exactly the regime SuSiE-Mix exploits.** Selecting among codings is a
statement about the *mean* structure \((\mu_1,\mu_2)\); it is well-posed,
closed-form, and is why the coding-mixture SER (the implemented `ser_mix`) works
and why heterogeneous-regulation loci are resolvable. None of the
identifiability problems below touch the coding-selection task.

---

## 3. What is *not* identifiable: the mixture interpretation (orders \(\theta^2\)–\(\theta^6\))

The per-individual claim is stronger than "non-additive": it asserts the
population is a **latent blend of distinct responder types**. That extra content
lives entirely in the **within-genotype higher moments**, not the means. Two
parametrizations with identical \((\mu_1,\mu_2)\) — one a genuine mixture, one a
single nonlinear dose–response — differ only there.

**Excess within-genotype variance.** A two-type blend (fraction \(1-\pi\)
additive with effect \(\theta\), fraction \(\pi\) recessive) at the
heterozygous stratum has mixing distribution of means \(\{(\theta,1-\pi),(0,\pi)\}\), so
\[
\boxed{\;\operatorname{Var}(y\mid g=1)=\sigma^2+\pi(1-\pi)\,\theta^2\;}
\tag{2}
\]
(verified to the digit). The mixture-specific signal is **second-order**:
\(O(\theta^2)\) in the variance, hence \(O(\theta^4)\) in information.

**Cumulants of the mixing distribution** (which become the excess cumulants of
\(y\mid g\) over a Gaussian):
\[
\kappa_2=\pi(1-\pi)\theta^2,\;\;
\kappa_3=\pi(1-\pi)(1-2\pi)\theta^3,\;\;
\kappa_4=\pi(1-\pi)\big(1-6\pi(1-\pi)\big)\theta^4 . \tag{3}
\]
Note \(\kappa_3=0\) at \(\pi=\tfrac12\) (symmetric blend): the skew signature
**vanishes** for equal mixing, leaving only the kurtosis.

### Hierarchy of nested models and the information for "is it a mixture?"

Confront the true mixture with progressively richer null classes and compute the
minimum KL (= per-sample expected log-likelihood gap = leading-order Fisher
information for the separating direction):

| Null model class (free params) | What it absorbs | leading discrepancy | min-KL scaling | verified slope |
|---|---|---|---|---|
| **Additive** \(\mu_g=\alpha g\) | linear mean | dominance deviation \(d=O(\theta)\) | \(O(\theta^2)\) | **2.00** |
| **Genotypic, homoscedastic** \(\mu_g\) free, common \(\sigma^2\) | full mean | excess variance \(\kappa_2=O(\theta^2)\) | \(O(\theta^4)\) | **3.99** |
| **Genotypic, heteroscedastic** \(\mu_g,\sigma_g^2\) free | mean + variance | skew \(\kappa_3=O(\theta^3)\) (or \(\kappa_4=O(\theta^4)\) if \(\pi=\tfrac12\)) | \(O(\theta^6)\) (resp. \(O(\theta^8)\)) | **5.98** (resp. **7.96**) |

The point of the table: **the per-individual mixture is identified *as a
mixture* only relative to a model that already explains the genotype means and
their heteroscedasticity.** Against that model the separating information is
\(O(\theta^6)\). Translated to required sample size \(n\sim c/\mathrm{KL}\):

| \(\theta\) | vs additive | vs genotypic (homosc.) | vs genotypic (heterosc.) |
|---|---|---|---|
| 0.1 | \(\sim 3\times10^4\) | \(\sim 2.6\times10^7\) | \(\sim 6\times10^{12}\) |
| 0.2 | \(\sim 7\times10^3\) | \(\sim 1.6\times10^6\) | \(\sim 2.6\times10^{10}\) |
| 0.3 | \(\sim 3\times10^3\) | \(\sim 3\times10^5\) | \(\sim 1\times10^9\) |

**Consequence for GTEx whole-blood (\(n\approx670\)):** distinguishing the
per-individual mixture from an ordinary heteroscedastic genotypic effect is off
by **3–9 orders of magnitude** in required \(n\). At realistic effect sizes it
is, for practical purposes, non-identifiable at single loci. (Simulated power of
an LRT for mixture-vs-genotypic: **0% at \(n=670\)**, ~25% only by \(n=5\times10^4\),
and a non-mixture heteroscedastic locus produces the *same* statistic — a false
"mixture" call.)

### Why this is non-regular, not just low-power

Testing \(H_0:\pi_D=\pi_R=0\) (single responder type) vs \(H_1\) (mixture) is a
**Hartigan/Davies boundary problem**: under \(H_0\) the effect sizes
\(\beta_D,\beta_R\) are *unidentified*, and the mixing weights sit on the
simplex boundary. The LRT therefore does **not** have a \(\chi^2\) limit; it
behaves like the supremum of a Gaussian process and can drift like
\(\log\log n\) (Hartigan 1985; Chen & Chen 2001; Liu & Shao 2003). Two practical
implications: (a) naive \(\chi^2\) or BIC calibration of a mixture test is
invalid; (b) any Bayes factor is **strongly prior-sensitive** in the
\((\beta_D,\beta_R)\) directions — improper priors give arbitrary-scale,
meaningless BFs. The prior is doing real inferential work, which motivates §4.

---

## 4. Regularization that falls back to additive

Because the likelihood is essentially flat in the mixture directions, a sensible
prior *should* dominate there — and the honest default is to fall back to
additive unless the data overwhelmingly object.

**Prior on the simplex.** Place \(\pi\sim\mathrm{Dirichlet}(a_A,a_D,a_R)\) with
\(a_A\gg a_D,a_R\) (additive-leaning), or a **spike-and-slab**
\[
\pi=(1,0,0)\ \text{w.p. } 1-\rho \quad(\text{pure additive}),\qquad
\pi\sim\mathrm{Dir}(1,1,1)\ \text{w.p. }\rho,
\]
which yields a clean posterior probability of "mixture vs additive" and an
automatic Occam penalty on the extra, poorly-identified parameters.

**Penalized-likelihood view.** Maximize
\(\ell(\pi,\beta,\sigma^2)-\lambda\,\mathrm{pen}(\pi)\) with \(\mathrm{pen}\)
pulling \(\pi\) toward the additive vertex (e.g. \(-\log\mathrm{Dir}\), or a
ridge on \((\pi_D,\pi_R)\)). Since the likelihood curvature in those directions
is \(O(n\theta^4)\), the penalty wins unless
\[
\boxed{\,n\,\theta^4 \;\gtrsim\; \lambda\,}\tag{4}
\]
i.e. the data override the additive default only once \(n\) exceeds
\(\lambda/\theta^4\) — exactly the scale from §3. This makes "favor additive"
quantitative rather than ad hoc.

**Structural priors that help.** Tying effects \(\beta_A=\beta_D=\beta_R\) (the
symmetric case in Fig. 1) or imposing sign/order constraints reduces the
unidentified directions — but pushes the model *toward* a genotypic mean model,
underscoring that the identifiable content is the mean, not the blend.
**Across loci**, partial pooling (a genome-wide hierarchical prior on the
mixing weights) is the only route to estimating a population-level mixture
fraction at GTEx \(n\); per-locus estimates are uninformative.

---

## 5. A calibrated Bayes-factor test for "is this CS a per-individual mixture?"

Operate on the **partial residual** \(r\) isolating one SuSiE single effect and
the genotype \(g\) of the credible set's lead variant (or marginalize over the
CS). Define three nested marginal likelihoods with **proper, scale-matched**
priors on effects (\(\beta\sim N(0,\sigma_0^2)\)) and a Dirichlet prior on
\(\pi\):

- \(M_{\mathrm{add}}\): additive, \(\mu_g=\alpha g\).
- \(M_{\mathrm{geno}}\): saturated 2-df genotypic mean (optionally heteroscedastic).
- \(M_{\mathrm{mix}}\): per-individual mixture (1), integrating over \(\pi,\beta,\sigma^2\).

Form two Bayes factors:
\[
\mathrm{BF}_{\mathrm{geno:add}}=\frac{p(r\mid M_{\mathrm{geno}})}{p(r\mid M_{\mathrm{add}})},
\qquad
\mathrm{BF}_{\mathrm{mix:geno}}=\frac{p(r\mid M_{\mathrm{mix}})}{p(r\mid M_{\mathrm{geno}})}.
\]

**Decision logic (the key separation):**

| \(\mathrm{BF}_{\mathrm{geno:add}}\) | \(\mathrm{BF}_{\mathrm{mix:geno}}\) | conclusion |
|---|---|---|
| \(\approx1\) | \(\approx1\) | additive — no action |
| large | \(\approx1\) | **non-additive but NOT a mixture**: a single nonlinear dose–response (dominance). *Do not label as mixture.* This is the case the proposal must avoid mis-calling. |
| large | large | genuine within-genotype heterogeneity — *necessary* evidence for a per-individual mixture (still confirm bimodality vs. plain heteroscedasticity). |

The discriminating quantity is \(\mathrm{BF}_{\mathrm{mix:geno}}\), **not**
\(\mathrm{BF}_{\mathrm{mix:add}}\): the latter is inflated by mean dominance and
will flag ordinary genotypic effects.

**Locally most powerful alternative (recommended at scale).** Near the additive
boundary the mixture's leading departure from \(M_{\mathrm{geno}}\) is the
excess within-genotype variance (2). The **score/dispersion test** for
within-genotype overdispersion is the locally most powerful test and avoids
fitting the unidentified mixture entirely:
\[
T=\sum_{g}\sum_{i:g_i=g}\Big[(y_i-\hat\mu_g)^2-\hat\sigma^2\Big]\ \big/\ \widehat{\mathrm{sd}},
\]
calibrated against the **genotypic null by simulation** (per O2), because the
boundary non-regularity (§3) invalidates analytic \(\chi^2\) calibration.

**Mandatory caveats baked into the procedure.** (a) Use informative,
scale-matched effect priors — improper priors make \(\mathrm{BF}_{\mathrm{mix:geno}}\)
meaningless. (b) Calibrate the threshold by null simulation under
\(M_{\mathrm{geno}}\) (including a *heteroscedastic* genotypic null), since a
plain heteroscedastic locus mimics the \(O(\theta^4)\) signal. (c) Report the
implied power: at GTEx \(n\), \(\mathrm{BF}_{\mathrm{mix:geno}}\) has essentially
no power for \(\theta\lesssim0.3\), so a per-locus null result is *expected*, not
evidence against the hypothesis — the test is only meaningful in large biobank
traits or via cross-locus pooling.

---

## 6. Bottom line

- **Coding selection** (additive vs dominant vs recessive; heterogeneous
  regulation) is a *mean-structure* problem, identified at \(O(\theta^2)\). It is
  the safe, implemented core (`MixSER::ser_mix` / `susie_mix`).
- **The per-individual mixture interpretation** is a *higher-moment* problem,
  identified only at \(O(\theta^4)\)–\(O(\theta^6)\) relative to a genotypic
  model, hence near non-identifiable at single loci for GTEx-scale data and
  confounded with ordinary heteroscedasticity.
- **Regularize toward additive** with a Dirichlet/spike-and-slab prior; the data
  override it only when \(n\theta^4\gtrsim\lambda\).
- **Test with \(\mathrm{BF}_{\mathrm{mix:geno}}\)** (not vs additive),
  simulation-calibrated against a heteroscedastic genotypic null, with proper
  effect priors. Lead empirical claims with heterogeneous regulation; treat the
  per-individual mixture as an extension contingent on this analysis.
