#' Placebo standard error with poorly-fitting placebos trimmed
#'
#' `synthdid`'s placebo variance permutes the control units, treats the last
#' one as if it were treated, re-estimates, and takes the spread of those
#' estimates. In a large heterogeneous donor pool most controls look nothing
#' like the treated unit, so many placebo draws are fitted badly and produce
#' large spurious "effects". Those draws fatten the tails of the reference
#' distribution and inflate the standard error.
#'
#' The symptom is measurable with in-time placebos, where the true effect is
#' zero by construction: in large heterogeneous pools `|t|` comes out far
#' more concentrated near zero than the `|N(0, 1)|` a calibrated interval
#' implies, while small homogeneous pools show no such compression.
#'
#' The fix follows the spirit of Abadie, Diamond and Hainmueller (2010):
#' restrict the placebo reference distribution to units the estimator can
#' actually fit, by discarding placebos whose pre-treatment fit is worse than
#' `trim_k` times the treated unit's. A placebo that cannot be fitted says
#' nothing about how variable a well-fitted unit is -- it is evidence of a
#' bad match, not of uncertainty.
#'
#' **This trims the inference distribution, not the donor pool.** The
#' counterfactual is still built from every donor, and the criterion is
#' pre-treatment fit, evaluated without reference to any unit's estimated
#' effect, so it cannot be used to select donors that produce a preferred
#' estimate.
#'
#' Fit is measured the way [fit_sdid()] measures it: RMSE of the treated-
#' minus-synthetic gap after removing its pre-period mean, on the
#' covariate-adjusted outcome (see [adjusted_outcome()]).
#'
#' # Relationship to synthdid
#'
#' Each placebo is refitted exactly as `synthdid:::placebo_se()` refits it:
#' `synthdid_estimate()` on the permuted panel, with the original fit's
#' options and its unit weights renormalised over the placebo controls as
#' the starting point. Because an `sc_estimate()` or `did_estimate()` stores
#' its fixed weights and options on the estimate, this reproduces those
#' estimators too. The draws consume the RNG exactly as `synthdid` does, so
#' under the same seed `trim_k = Inf` returns `synthdid`'s
#' `sqrt(vcov(est, method = "placebo"))` to machine precision.
#'
#' @param est A `synthdid_estimate` object.
#' @param replications Number of placebo permutations.
#' @param trim_k Keep placebos whose pre-period trend RMSE is at most
#'   `trim_k` times the treated unit's. `Inf` disables trimming.
#' @param return_draws Also return the per-draw estimates. Consumes no extra
#'   randomness, so a run with and without it gives the same `se`.
#'
#' @return A list: `se`, `n_kept`, `n_total`, `treated_rmse`,
#'   `placebo_rmse_median`, and -- when `return_draws` is `TRUE` -- `draws`,
#'   a data frame of `draw`, `att`, `rmse`, `keep` and `placebo_unit` (the
#'   control that stood in for the treated unit).
#'
#' @references Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic
#'   Control Methods for Comparative Case Studies. *Journal of the American
#'   Statistical Association* 105(490), 493-505.
#'
#' @examples
#' \donttest{
#' data("california_prop99", package = "synthdid")
#' setup <- synthdid::panel.matrices(california_prop99)
#' est <- synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' set.seed(1)
#' placebo_se_trimmed(est, replications = 50, trim_k = 5)$se
#' }
#'
#' @export
placebo_se_trimmed <- function(est, replications = 200L, trim_k = 5,
                               return_draws = FALSE) {
  setup   <- attr(est, "setup")
  opts    <- attr(est, "opts")
  weights <- attr(est, "weights")
  N0 <- setup$N0
  N1 <- nrow(setup$Y) - N0
  if (N0 <= N1) {
    cli::cli_abort("Need more controls than treated units for a placebo SE.")
  }
  if (N1 != 1L) {
    cli::cli_abort(c(
      "Trimmed placebo inference is defined for a single treated unit; got {N1}.",
      i = "Use {.code vcov(est, method = \"placebo\")} from synthdid instead."
    ))
  }

  # Mirrors synthdid:::sum_normalize, reimplemented to avoid an unexported
  # dependency.
  norm1 <- function(x) if (sum(x) != 0) x / sum(x) else rep(1 / length(x), length(x))

  trend_rmse <- function(e) {
    s   <- attr(e, "setup")
    gap <- pre_gap(adjusted_outcome(e), s$N0, s$T0, attr(e, "weights")$omega)
    sqrt(mean((gap - mean(gap))^2))
  }
  treated_rmse <- trend_rmse(est)

  draws <- vapply(seq_len(replications), function(i) {
    ind <- sample(seq_len(N0))
    n0p <- length(ind) - N1
    wb  <- weights
    wb$omega <- norm1(weights$omega[ind[seq_len(n0p)]])
    e <- do.call(synthdid::synthdid_estimate, c(
      list(Y = setup$Y[ind, , drop = FALSE], N0 = n0p, T0 = setup$T0,
           X = setup$X[ind, , , drop = FALSE], weights = wb), opts
    ))
    # The permutation puts the stand-in treated unit last; recording it
    # costs no randomness and lets the variance be attributed to donors.
    c(att = as.numeric(e), rmse = trend_rmse(e), unit = ind[n0p + 1L])
  }, numeric(3))

  att  <- draws["att", ]
  rmse <- draws["rmse", ]
  unit <- draws["unit", ]
  keep <- is.finite(att) & is.finite(rmse) & rmse <= trim_k * treated_rmse
  if (sum(keep) < 5) {
    cli::cli_warn(
      "Only {sum(keep)} of {replications} placebos survive trim_k = {trim_k}; falling back to untrimmed."
    )
    keep <- is.finite(att)
  }

  n <- sum(keep)
  out <- list(
    # synthdid's finite-sample scaling, sqrt((n - 1) / n) * sd.
    se = sqrt((n - 1) / n) * stats::sd(att[keep]),
    n_kept = n,
    n_total = replications,
    treated_rmse = treated_rmse,
    placebo_rmse_median = stats::median(rmse, na.rm = TRUE)
  )
  if (return_draws) {
    ids <- rownames(setup$Y)
    out$draws <- data.frame(
      draw = seq_len(replications), att = att, rmse = rmse, keep = keep,
      placebo_unit = if (is.null(ids)) unit else ids[unit]
    )
  }
  out
}
