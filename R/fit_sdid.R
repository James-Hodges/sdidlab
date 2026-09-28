#' Fit one synthetic-DiD specification and return its full result set
#'
#' Wraps `synthdid::synthdid_estimate()` (or `sc_estimate()` /
#' `did_estimate()`) with the inference and diagnostics a single-treated-unit
#' design needs on every specification.
#'
#' Three things it standardises:
#'
#' * **Placebo inference.** With a single treated unit the jackknife variance
#'   is undefined -- it would require dropping the treated unit -- so the
#'   placebo variance is the default. With `trim_k = Inf` it is exactly
#'   `synthdid`'s own `vcov(est, method = "placebo")`; see
#'   [placebo_se_trimmed()] for what trimming changes.
#' * **The right pre-period fit statistic.** SDID carries unit fixed effects
#'   and therefore tolerates a *constant* level gap between treated and
#'   synthetic control; absorbing it is what the DiD half does. Raw RMSE
#'   against the synthetic path penalises exactly what the estimator is
#'   designed to ignore. `trend_rmse` is the RMSE of the treated-minus-
#'   synthetic gap **after removing its pre-period mean**, which is the
#'   parallel-trends assumption; `level_gap` reports the absorbed offset
#'   separately. Both are computed on the outcome the estimator actually
#'   fits: when covariates are passed as `X`, that is `Y - X beta`, not `Y`.
#' * **A confidence interval.** When every estimate is null the bound is the
#'   result, and a bound cannot be read off a p-value.
#'
#' Placebo standard errors are simulated: the same specification moves by a
#' few percent between runs depending on the RNG state. Set a seed before
#' calling if results must be reproducible.
#'
#' @param Y Outcome matrix, units x periods, controls first and the treated
#'   unit last -- exactly what `synthdid::panel.matrices()` returns.
#' @param N0 Number of control units (the first `N0` rows of `Y`).
#' @param T0 Number of pre-treatment periods (the first `T0` columns).
#' @param X Optional unit x period x covariate array, handed to `synthdid`,
#'   which estimates the covariate coefficients jointly with the weights. To
#'   residualise on untreated cells first instead, use
#'   [project_out_covariates()] and pass the result as `Y`.
#' @param replications Placebo replications. 200 is `synthdid`'s default.
#' @param conf_level Confidence level for the interval.
#' @param trim_k Passed to [placebo_se_trimmed()]. Placebos fitting worse
#'   than `trim_k` times the treated unit's pre-period trend RMSE are dropped
#'   from the reference distribution. `Inf` reproduces `synthdid`'s untrimmed
#'   placebo SE.
#' @param se_method `"placebo"` (default) or `"none"`. Placebo inference
#'   refits the model once per replication and dominates the runtime, so
#'   sweeps that only need point estimates should pass `"none"`, which
#'   returns `NA` for `se`, `t`, `ci_lo` and `ci_hi`.
#' @param estimator `"sdid"` (default), `"sc"` for synthetic control, or
#'   `"did"` for plain difference-in-differences, all as `synthdid` defines
#'   them. The placebo distribution is refitted with the same estimator.
#' @param keep How much to return beyond the summary. `"none"` (default)
#'   returns the summary fields only. `"est"` adds the `synthdid_estimate`
#'   object, which every plot needs. `"all"` adds the individual placebo
#'   draws. Neither consumes extra randomness, so the estimate and standard
#'   error are unchanged.
#'
#' @return A list with `att`, `se`, `t`, `ci_lo`, `ci_hi`, `trend_rmse`,
#'   `level_gap`, `N0`, `T0`, `n_periods`, `omega` (unit weights, named by
#'   the control rownames of `Y`), `lambda` (time weights), `estimator`, and
#'   the placebo trimming record `trim_k`, `placebos_kept`, `placebos_total`.
#'   With `keep = "est"`, also `est`; with `keep = "all"`, also
#'   `placebo_draws`, `treated_rmse`, `placebo_rmse_median`, `unit_ids` and
#'   `periods`.
#'
#' @references Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W.
#'   and Wager, S. (2021). Synthetic Difference-in-Differences. *American
#'   Economic Review* 111(12), 4088-4118.
#'
#' @examples
#' \donttest{
#' data("california_prop99", package = "synthdid")
#' setup <- synthdid::panel.matrices(california_prop99)
#' set.seed(1)
#' r <- fit_sdid(setup$Y, setup$N0, setup$T0, replications = 50)
#' c(att = r$att, lo = r$ci_lo, hi = r$ci_hi)
#' }
#'
#' @export
fit_sdid <- function(Y, N0, T0, X = NULL, replications = 200L,
                     conf_level = 0.95, trim_k = 5,
                     se_method = c("placebo", "none"),
                     estimator = c("sdid", "sc", "did"),
                     keep = c("none", "est", "all")) {
  se_method <- match.arg(se_method)
  estimator <- match.arg(estimator)
  keep      <- match.arg(keep)
  if (N0 < 2) cli::cli_abort("Need at least 2 control units; got {N0}.")
  if (T0 < 2) cli::cli_abort("Need at least 2 pre-periods; got {T0}.")

  fit_fun <- switch(estimator,
    sdid = synthdid::synthdid_estimate,
    sc   = synthdid::sc_estimate,
    did  = synthdid::did_estimate
  )
  est <- if (is.null(X)) fit_fun(Y, N0, T0) else fit_fun(Y, N0, T0, X = X)

  trim <- if (se_method == "none") {
    NULL
  } else {
    placebo_se_trimmed(est, replications = replications, trim_k = trim_k,
                       return_draws = identical(keep, "all"))
  }
  se <- if (is.null(trim)) NA_real_ else trim$se

  w   <- attr(est, "weights")$omega
  gap <- pre_gap(adjusted_outcome(est), N0, T0, w)
  level_gap <- mean(gap)

  att <- as.numeric(est)
  z   <- stats::qnorm(1 - (1 - conf_level) / 2)

  out <- list(
    att       = att,
    se        = se,
    t         = att / se,
    ci_lo     = att - z * se,
    ci_hi     = att + z * se,
    # Parallel-trends fit, with the constant offset SDID absorbs removed.
    trend_rmse = sqrt(mean((gap - level_gap)^2)),
    level_gap = level_gap,
    N0        = N0,
    T0        = T0,
    n_periods = ncol(Y),
    omega     = stats::setNames(as.numeric(w), rownames(Y)[seq_len(N0)]),
    lambda    = as.numeric(attr(est, "weights")$lambda),
    estimator = estimator,
    trim_k    = trim_k,
    placebos_kept  = if (is.null(trim)) NA_integer_ else trim$n_kept,
    placebos_total = if (is.null(trim)) NA_integer_ else trim$n_total
  )

  # `est` and `placebo_draws` are always present as names, holding NULL when
  # not kept: `x$est` on a list without an exact `est` element would
  # partially match `estimator` and silently return a string.
  out["est"] <- list(if (keep %in% c("est", "all")) est else NULL)
  out["placebo_draws"] <- list(if (identical(keep, "all")) trim$draws else NULL)
  if (identical(keep, "all")) {
    out$treated_rmse        <- if (is.null(trim)) NA_real_ else trim$treated_rmse
    out$placebo_rmse_median <- if (is.null(trim)) NA_real_ else trim$placebo_rmse_median
    out$unit_ids            <- rownames(Y)
    out$periods             <- colnames(Y)
  }
  out
}

#' The outcome a synthdid estimate was actually fitted on
#'
#' `synthdid` estimates covariate coefficients `beta` jointly with the
#' weights and computes the effect on `Y - X beta`. Any diagnostic of the fit
#' -- the pre-period gap, the trend RMSE, a gap-path plot -- has to be
#' computed on that adjusted outcome, or it describes a fit the estimator
#' never made. Without covariates this is `Y` itself.
#'
#' @param est A `synthdid_estimate` object.
#'
#' @return A units x periods matrix.
#'
#' @examples
#' \donttest{
#' data("california_prop99", package = "synthdid")
#' setup <- synthdid::panel.matrices(california_prop99)
#' est <- synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' identical(adjusted_outcome(est), setup$Y)
#' }
#'
#' @export
adjusted_outcome <- function(est) {
  setup <- attr(est, "setup")
  beta  <- attr(est, "weights")$beta
  X     <- setup$X
  if (is.null(X) || length(dim(X)) < 3 || dim(X)[3] == 0 || !length(beta)) {
    return(setup$Y)
  }
  Xb <- matrix(0, nrow(setup$Y), ncol(setup$Y))
  for (k in seq_len(dim(X)[3])) Xb <- Xb + beta[[k]] * X[, , k]
  setup$Y - Xb
}

# Treated-minus-synthetic gap over the pre-period.
pre_gap <- function(Y, N0, T0, omega) {
  as.numeric(Y[N0 + 1L, seq_len(T0)]) -
    as.numeric(t(Y[seq_len(N0), seq_len(T0), drop = FALSE]) %*% omega)
}
