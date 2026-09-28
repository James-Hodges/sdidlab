# sdidlab

Synthetic difference-in-differences for a single treated unit, with the
inference and diagnostics that design needs.

`sdidlab` is built around the [synthdid](https://github.com/synth-inference/synthdid)
estimator of Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021). Where it
reproduces synthdid, it does so exactly: point estimates, unit and time
weights, and placebo standard errors match synthdid to machine precision,
and the test suite enforces it. On top of that it adds:

- **Trimmed placebo inference.** In a large, heterogeneous donor pool, many
  placebo draws are fitted badly and inflate the placebo standard error.
  `placebo_se_trimmed()` drops placebos whose pre-period fit is more than
  `trim_k` times worse than the treated unit's. It trims the *inference*
  distribution, never the donor pool.
- **Fit diagnostics that respect SDID.** SDID absorbs a constant level gap
  between the treated unit and its synthetic control, so `trend_rmse`
  measures the pre-period gap *after* removing its mean, on the
  covariate-adjusted outcome when covariates are used.
- **Two covariate routes.** Hand covariates to synthdid, or residualise the
  outcome on untreated cells first (Kranz 2022) with
  `project_out_covariates()`.
- **A bias-corrected synthetic control.** `fit_scm()` implements the
  Abadie–Diamond–Hainmueller estimator with predictor matching, the
  Abadie–L'Hour bias correction, an optional Abadie–L'Hour penalty that
  makes the weights unique, and ranked-RMSPE permutation inference through
  `scm_placebo_permutation()` and `rmspe_pvalue()`. The weights are solved
  exactly by an active-set method, and an exact pre-period fit -- where
  the weights, and therefore the estimate, are not unique -- is flagged
  rather than silently resolved.

## Installation

synthdid is not on CRAN, so install both from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("James-Hodges/sdidlab")
```

The `Remotes` field pins synthdid to the commit the parity tests were run
against.

## Example: California's Proposition 99

``` r
library(sdidlab)

data("california_prop99", package = "synthdid")
setup <- synthdid::panel.matrices(california_prop99)

set.seed(1)
fit <- fit_sdid(setup$Y, setup$N0, setup$T0, replications = 200, trim_k = 5)

fit$att                 # -15.6 packs per capita, as in Arkhangelsky et al.
c(fit$ci_lo, fit$ci_hi) # placebo-based 95% interval
fit$trend_rmse          # pre-period fit, level gap removed
fit$placebos_kept       # how many of the 200 placebos survived trimming
```

With `trim_k = Inf` the standard error is synthdid's own placebo standard
error, draw for draw under the same seed.

The synthetic-control alternative:

``` r
scm <- fit_scm(setup$Y, setup$N0, setup$T0, v_method = "mspe")
scm$att                 # bias-corrected ATT

perm <- scm_placebo_permutation(setup$Y, setup$N0, setup$T0, v_method = "mspe")
r <- perm$results
rmspe_pvalue(setNames(r$rmspe_ratio_raw, r$unit), r$att, r$unit[r$treated])
```

## Validation

The test suite checks sdidlab against independent implementations, not just
against itself:

| What | Checked against |
|---|---|
| SDID, SC and DiD estimates and weights | `synthdid`, to 1e-12 |
| Placebo standard errors (untrimmed), with and without covariates | `synthdid::vcov(method = "placebo")`, to 1e-10 under the same seed |
| Covariate residualisation | Two-way fixed-effects OLS on untreated cells, and `fixest` |
| Synthetic-control weights, plain and penalised | The KKT optimality conditions; an independent long-run solve; `Synth` with the same `V` |
| Bias correction | An independent per-period OLS computation |

## References

Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic Control Methods
for Comparative Case Studies. *Journal of the American Statistical
Association* 105(490), 493–505.

Abadie, A. and L'Hour, J. (2021). A Penalized Synthetic Control Estimator for
Disaggregated Data. *Journal of the American Statistical Association*
116(536), 1817–1834.

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W. and Wager, S.
(2021). Synthetic Difference-in-Differences. *American Economic Review*
111(12), 4088–4118.

Kranz, S. (2022). Synthetic Difference-in-Differences with Time-Varying
Covariates. Working paper.

## License

MIT © James Hodges
