# SER-Mix modelling notes

## Two distinct models — keep them separate

The CHIMERA proposal mixes two ideas. They need different machinery:

1. **Coding-model mixture (implemented here, `MixSER::ser_mix`).**
   A single effect ranges over `(variant j, coding m ∈ {A,D,R})`. Each causal
   variant has *one* coding; different variants at a locus may differ. This is
   the **heterogeneous-regulation** model (H2). It is **closed-form** — every
   `(j,m)` is scored with the exact SuSiE single-effect Bayes factor, and the
   posterior `alpha` normalises over all `p × M` hypotheses. No EM, no latent
   per-individual assignment, drops straight into IBSS. Reduces to standard
   SuSiE when only the additive coding is supplied.

2. **Per-individual mixture (Eq. 1 of the proposal, H1).**
   Within a *single* SNP, individual `i` follows additive *or* dominant *or*
   recessive. This is near non-identifiable at one locus: three genotype classes
   give only three distinguishable conditional means, but `(π_A,π_D,π_R,β_A,β_D,β_R)`
   has more free parameters. It is also observationally close to a standard 2-df
   genotypic (dominance-deviation) model. **Recommend demoting H1 to an
   extension** contingent on a formal identifiability result, and leading the
   method paper with model 1.

The old exploratory scripts (`claude_model.R`, `attempt_mixture_of_effect.R`)
implemented an EM/BIC version of model 2 plus matrix stacking. They are
superseded by the closed-form `MixSER` implementation for the locus-level
question.

## What the simulation shows (`simulation_sermix.R`)

- **S1 coding recovery** — SER-Mix labels a single causal variant's mode of
  action correctly at high rate.
- **S2 heterogeneous regulation** — two causal variants, one additive + one
  recessive. Compares: standard additive SuSiE, naive *stacked* `[A|D|R]` SuSiE
  (the tempting but wrong fix — inflates / mis-resolves CS), and SuSiE-Mix,
  which recovers both signals with the correct coding.
- **S3 null calibration** — under pure additivity SER-Mix does not spuriously
  prefer a non-additive coding.

## Running

```r
# from this folder
source("simulation_sermix.R")
```

Needs `susieR` and the `MixSER` package (auto-loaded from
`C:/Document/Serieux/Travail/Package/MixSER` via `devtools::load_all`, or
`source()` of its `R/` files as a fallback).

## Numerical validation already done

The SER-Mix log Bayes factor and posterior moments were checked against an exact
multivariate-normal marginal likelihood (agreement to ~1e-13), and the three
test scenarios above were reproduced in an independent implementation
(coding recovery and 2-CS heterogeneous regulation both correct). See
`MixSER/tests/testthat/test-ser_mix.R`.
```
```
