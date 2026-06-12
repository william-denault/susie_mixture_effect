# An implicit-weight mixture-effect SER

*A variational/empirical-Bayes reformulation that removes the Dirichlet prior on
the model weights and learns them the way `ash` learns mixture proportions —
by maximum marginal likelihood, a convex problem solved by `mixsqp`.*

This note responds to two requests:

1. *"Get rid of the Dirichlet part; let's do it implicitly and add it as default
   parameters."* — The Dirichlet prior on the mixing weights in `ser_pim_vb` is
   replaced by an **empirical-Bayes estimate of the weights**, obtained by
   maximising the marginal likelihood. This is exactly the device behind `ashr`
   (Stephens 2017) and is solved by the convex optimiser of `mixsqp`
   (Kim, Carbonetto, Stephens & Anitescu 2020).
2. *"For each SNP run additive / dominant / recessive separately, optimise their
   hyperparameters as in SuSiE, then learn the prior probability of each model
   from their Bayes factors; finally compute the posterior of each model and the
   individual probability of belonging to each."* — This is precisely the
   architecture below. The per-coding hyperparameters are the SuSiE
   empirical-Bayes prior variances; the "prior probability of each model" is the
   `ash`/`mixsqp` weight vector; the posteriors fall straight out by Bayes' rule.

Throughout, work in residual-standard-deviation units (`σ = 1`) and write the
standardised effect `θ = β/σ`.

---

## 0. The one idea

`ashr` never *tests* "signal vs. null." It fits **one** over-arching mixture
prior `g(β) = Σ_k w_k f_k(β)` over a fixed family of components `{f_k}`, and
estimates the weights `w` by maximum marginal likelihood:

> `ŵ = argmax_{w ∈ Δ} Σ_t log( Σ_k w_k L_{tk} )`,  `L_{tk} = p(data_t | f_k)`.   (★)

Because the marginal is **linear in `w`**, the log-likelihood is **concave**, so
(★) is a convex program over the simplex `Δ`. `mixsqp` solves it by sequential
quadratic programming with an active set — orders of magnitude faster than EM and
exact (Kim et al. 2020). The posterior weight on component `k` for observation
`t` is then just `w_k L_{tk} / Σ_{k'} w_{k'} L_{tk'}`.

Everything below is (★) with the components reinterpreted. **No Dirichlet, no
Occam penalty by hand** — the Occam factor is supplied by sharing one `w` across
many observations (`mixsqp` drives a useless component's weight to `0`), and by
the empirical-Bayes shrinkage of the per-component effect variances.

---

## 1. Two mixture levels — keep them apart

A genotype `g ∈ {0,1,2}` has three classical codings
`xᴬ(g)=g`, `xᴰ(g)=1{g≥1}`, `xᴿ(g)=1{g=2}`. There are two genuinely different
"mixtures of effects," and they have very different identifiability (see
`identifiability_per_individual_mixture.md`):

- **Level 1 — coding-model mixture (mean structure, identifiable at `O(θ²)`).**
  Each SNP acts through *one* coding `m ∈ {A,D,R}`, but we do not know which, and
  different SNPs in a locus may differ. This is a statement about the conditional
  *means* `E[y|g]`; it is well posed, closed-form, and is the safe core of the
  method. **This is what should drive every model decision.**

- **Level 2 — per-individual mixture (higher moment, weakly identified).**
  Within *one* SNP, individual `i` responds through a latent coding `z_i`, with
  population fractions `π`. The "individual probability of belonging to each
  model" lives here as the responsibilities `r_{im} = q(z_i = m)`. The signal is
  `O(θ²)` in a within-genotype *variance*, hence `O(θ⁴)` in information, and is
  confounded with ordinary heteroscedasticity. It is **only** estimable by
  pooling across loci.

The same convex device (★) handles both. The crucial discipline — which also
resolves the `VEB-pi` pathology recorded in memory — is to let **Level 1** make
the additive/dominant/recessive/mixture call, and to treat **Level 2** as a
*descriptive* posterior with `π` **pooled across loci**, never as a per-locus
model-selection statistic.

---

## 2. Level 1 — the coding-model mixture SER (the recommended core)

### 2.1 Per-coding hyperparameters, exactly as in SuSiE

Inside one IBSS single effect we hold a centred partial residual `y` (length `n`)
and `p` candidate variants. For variant `j` and coding `m`, fit the one-predictor
Bayesian regression

```
y = x_j^m β_{jm} + ε,   β_{jm} ~ N(0, σ²₀ₘ),   ε ~ N(0, σ² I).
```

This is the SuSiE single-effect with the **closed-form** log Bayes factor and
posterior moments of Wang et al. (2020) — already implemented in
`.ser_one_coding`. Each coding gets **its own** prior variance `σ²₀ₘ`, estimated
by empirical Bayes exactly as SuSiE estimates its prior variance: maximise the
single-effect log marginal likelihood for that coding,

```
σ̂²₀ₘ = argmax_V  log Σ_j γ_j · BF_{jm}(V),        (2.1)
```

with `γ_j` the prior over variants (uniform `1/p` by default). One scalar 1-D
optimisation per coding. (Generalisation: replace the single `σ²₀ₘ` by an `ash`
grid of scales `{σ²₀ₘₖ}` with their own weights — the math is unchanged, the
components just multiply.)

### 2.2 The model prior, learned implicitly (this replaces the Dirichlet)

Add an explicit **null** component `m = 0` with `BF_{j0} ≡ 1`. Put a weight
vector `w = (w₀, wᴬ, wᴰ, wᴿ)` on the simplex — *this* is "the prior probability
of each model." Instead of a Dirichlet prior we estimate `w` by maximum marginal
likelihood. Treating each of the `p` candidate variants as an observation, the
SER objective is exactly (★):

```
ŵ = argmax_{w ∈ Δ}  Σ_{j=1}^p  log( Σ_{m∈{0,A,D,R}} w_m · BF_{jm} ).   (2.2)
```

`BF_{jm}` is the relative likelihood of the four hypotheses for variant `j`;
multiplying any row by a constant leaves (2.2) unchanged, so using Bayes factors
(rather than raw marginal likelihoods) is exact. (2.2) is **convex** and is the
canonical `mixsqp` problem. It is self-regularising: a null variant has
`BF_{jm} ≈ 1` for every `m`, contributes `log(Σ_m w_m) = 0`, and so carries no
information about `w` — only the handful of signal variants move it.

Two scopes (the implementation exposes both):

- **Per-SER** (default for a stand-alone SER): solve (2.2) over the `p` variants
  in the current single effect. Weights are re-learned each IBSS pass.
- **Genome-wide pooled** (recommended at scale): stack the `BF_{tm}` of many
  tests/SNPs into one tall matrix and solve (★) once; freeze `ŵ` inside SuSiE.
  This is the across-loci pooling the identifiability note argues for.

### 2.3 Posteriors — by Bayes' rule, no extra machinery

With `ŵ` and the per-coding `BF`, the joint posterior over `(variant, model)` is

```
α_{jm} = γ_j w_m BF_{jm} / Σ_{j'm'} γ_{j'} w_{m'} BF_{j'm'}.    (2.3)
```

- **Variant PIP** (the SER output IBSS consumes): `α_j = Σ_m α_{jm}`.
- **Posterior over the model of variant `j`**: `P(m | j) = α_{jm} / α_j`
  — the user's "posterior of each of these models," with `m=0` being "this
  variant is null."
- **SER log Bayes factor**: `log Σ_{jm} γ_j w_m BF_{jm}` (excluding `m=0`
  appropriately), the quantity SuSiE needs for residual-variance and
  convergence bookkeeping.
- **Posterior mean fitted value** for the IBSS residual update:
  `Σ_j Σ_m α_{jm} · x_j^m · E[β_{jm} | y]`.

When only the additive coding is supplied (`M=1`, `w₀` free), (2.2)–(2.3) reduce
**exactly** to the standard additive SuSiE SER. This is the drop-in property: the
existing `ser_mix` is the special case of fixed uniform `w`; we have only made
`w` empirical-Bayes.

---

## 3. Level 2 — per-individual responsibilities without a Dirichlet

This is the rewrite of `ser_pim_vb`. For a single chosen variant with centred
codings `x^m` (length `n`):

```
z_i ~ Categorical(π),   y_i | z_i = m  ~ N(x_i^m β_m, σ²),   β_m ~ N(0, σ²₀ₘ).
```

Mean-field posterior `q({z_i}) q({β_m})`. Treat `π` and `σ²₀ₘ` as
**hyperparameters estimated by empirical Bayes** — there is no `q(π)` and hence
no Dirichlet KL term. The variational free energy is

```
F = Σ_i Σ_m r_{im} [ log π_m + E_q log N(y_i; x_i^m β_m, σ²) − log r_{im} ]
      − Σ_m KL( q(β_m) ‖ N(0, σ²₀ₘ) ).                                   (3.1)
```

### 3.1 Coordinate ascent (CAVI), Dirichlet-free

With responsibilities `r_{im} = q(z_i = m)` and `S_m = Σ_i r_{im} (x_i^m)²`:

```
q(β_m):   s²_m = ( σ²₀ₘ⁻¹ + σ⁻² S_m )⁻¹,   μ_m = s²_m σ⁻² Σ_i r_{im} x_i^m y_i.   (3.2)

q(z_i):   log r_{im} ∝ log π_m
                       − (1/2σ²)( −2 x_i^m y_i μ_m + (x_i^m)² (μ_m² + s²_m) ).    (3.3)

EB σ²₀ₘ:  σ̂²₀ₘ = μ_m² + s²_m.                                                  (3.4)

EB π   :  the stationary point of (3.1) on the simplex is the mixsqp solution    (3.5)
          π̂ = argmax_{π ∈ Δ} Σ_i log( Σ_m π_m L_{im} ),  L_{im} = exp(E_q log N(y_i; x_i^m β_m, σ²)).
```

Equation (3.5) is the key substitution. The Lagrange stationarity condition of
(3.1) in `π` is the familiar EM update `π_m = (Σ_i r_{im})/n`; **`mixsqp` solves
the same convex profile problem globally and fast**, and is the natural home for
the *pooled* estimate (stack `L_{im}` over all individuals at all loci). So the
update is literally: *drop the Dirichlet block, and where the old code computed
`Elogpi` from `digamma`, instead set `log π_m` from the `mixsqp` weights.*

`F` increases monotonically and is the model's evidence lower bound; it is what we
report for this variant's mixture fit.

### 3.2 Why pooling — and why this is the honest fix to the `VEB-pi` problem

Memory records: *point-estimating `π` always calls "mixture"; integrate `π`
(Dirichlet, `method="vb"`) for the decision.* The reason is real: within **one**
locus the mixture model nests the pure ones (`π` at a vertex), so
`max_π`(marginal) `≥` any pure marginal — a per-locus point estimate of `π`
**cannot** adjudicate pure-vs-mixture; the larger model always wins. Replacing
the Dirichlet by a *per-locus* `mixsqp` estimate of `π` would inherit exactly this
pathology.

The resolution is structural, and it is the same one that makes `ashr` valid:

1. **Decisions live at Level 1**, where the question is about identifiable mean
   structure (`O(θ²)`), is convex, and has a proper marginal likelihood. The
   additive/dominant/recessive/**mixture** posterior is the Level-1 quantity
   (2.3), in which the mixture is one component scored by *its own marginal
   likelihood with `π` integrated/EB-pooled* — never by a per-locus
   `max_π`.

2. **`π` is pooled across loci.** Estimate **one** `π` genome-wide by `mixsqp`
   (3.5 with `L_{im}` stacked over all loci). A single pure-additive locus then
   contributes its likelihood under the *shared* `π` and cannot manufacture a
   mixture on its own. This is the only regime in which the per-individual
   fractions are identified (`identifiability_per_individual_mixture.md`, §4–6),
   and it is precisely the cross-locus partial pooling that note prescribes.

3. **Responsibilities are reported, not used to select.** Given the pooled `π̂`
   and the variant's coding fits, the per-individual posterior

   ```
   r_{im} = π̂_m N(y_i; x_i^m μ_m, σ²) / Σ_{m'} π̂_{m'} N(y_i; x_i^{m'} μ_{m'}, σ²)   (3.6)
   ```

   is the "individual probability of belonging to each model." It is a clean,
   interpretable output, with the honest caveat that at GTEx-scale `n` it is
   informative only after pooling and for `θ` not too small.

In short: the Dirichlet was doing two jobs — supplying an Occam factor and
regularising `π`. The first job moves to **Level 1's pooled marginal likelihood**;
the second to **cross-locus EB pooling of `π`**. Neither needs a Dirichlet, and
both have sensible defaults (next section).

---

## 4. Defaults (the "implicit, with default parameters" part)

| quantity | default | rationale |
|---|---|---|
| coding family | `{additive, dominant, recessive}` (+ explicit null) | the three classical codings; null makes the SER honest |
| `σ²₀ₘ` (per-coding prior var.) | EB via (2.1), start `0.2`, search `[1e-6, 1e2]` | SuSiE's own prior-variance estimation, per coding |
| model weights `w` (Level 1) | **estimated** by `mixsqp` (2.2); fallback start `(w₀,wᴬ,wᴰ,wᴿ)=(0.9,0.07,0.02,0.01)` | additive-leaning so that, absent evidence, the SER behaves like additive SuSiE |
| weight scope | `per_ser` for a stand-alone SER; `pooled` genome-wide when available | pooling is the identified regime |
| `π` (Level 2) | **estimated, pooled**, by `mixsqp` (3.5); start at additive vertex | no Dirichlet; collapses to additive unless data object across loci |
| solver | `mixsqp::mixsqp`; EM fallback if the package is absent | convex either way |

The additive-leaning fallback start is the *implicit* successor to the old
`dirichlet_prior = c(1,1,1)` / `a0 < 1` knobs: it sets where the model sits when
the data are silent, but is overridden by the data through (2.2)/(3.5) rather than
fixed.

---

## 5. What ships

- `ser_mix_eb()` — Level-1 drop-in SER with per-coding EB hyperparameters and
  `mixsqp`-estimated model weights (`per_ser` / `pooled` / `fixed`). Returns
  `alpha`, `coding_post`, `pip`, `lbf_model`, `fitted`, `weights`, `V` per coding;
  identical interface to `ser_mix`, so `susie_mix` / the `susie2` adapter consume
  it unchanged. Reduces to additive SuSiE when given only the additive coding.
- `ser_pim_vb_eb()` — Level-2 per-individual VB, Dirichlet removed; `π` by
  `mixsqp` with an optional **pooled** `π` passed in (the genome-wide estimate),
  returning the responsibilities `r_{im}` of (3.6) and the free energy `F`.
- `estimate_mixture_weights()` — the shared convex solver (`mixsqp`, EM fallback)
  for (★), used by both levels and for the genome-wide pooled estimate.

See `MixSER/R/ser_mix_eb.R`, `MixSER/R/ser_pim_vb_eb.R`, and the verification
script `MixSER/tests/verify_implicit_weights.R`.

---

### References

- M. Stephens (2017). *False discovery rates: a new deal.* Biostatistics 18(2),
  275–294. doi:10.1093/biostatistics/kxw041. — the `ash` unimodal mixture and the
  convex maximum-marginal-likelihood estimation of mixture weights.
- Y. Kim, P. Carbonetto, M. Stephens & M. Anitescu (2020). *A fast algorithm for
  maximum likelihood estimation of mixture proportions using sequential quadratic
  programming.* J. Comput. Graph. Stat. 29(2), 261–273.
  doi:10.1080/10618600.2019.1689985. — `mixsqp`, the convex solver for (★).
- G. Wang, A. Sarkar, P. Carbonetto & M. Stephens (2020). *A simple new approach
  to variable selection in regression, with application to genetic fine mapping.*
  JRSS-B 82(5), 1273–1300. — the SuSiE single-effect regression and IBSS.
