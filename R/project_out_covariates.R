#' Residualise the outcome on covariates using untreated cells only
#'
#' The covariate adjustment of Kranz (2022): regress the outcome on the
#' covariates with unit and time fixed effects, using only the **untreated**
#' unit-periods -- every control-unit cell plus the treated unit's
#' pre-period -- then subtract the covariates' fitted contribution from the
#' whole outcome matrix. Fitting on untreated cells keeps the treatment
#' effect out of the covariate coefficients. The residualised matrix is then
#' passed to [fit_sdid()] as `Y` with no `X`.
#'
#' This is the alternative to handing `X` to `synthdid`, which estimates the
#' coefficients jointly with the synthetic-control weights. The two can give
#' materially different answers when the covariates are strongly related to
#' the outcome; comparing them is a useful robustness check.
#'
#' The coefficients are estimated once, on the real treatment assignment.
#' Placebo inference on the residualised matrix therefore holds them fixed
#' rather than re-estimating them for each placebo.
#'
#' @param Y Outcome matrix, units x periods, controls first and the treated
#'   unit last.
#' @param X Units x periods x covariates array aligned with `Y` (a units x
#'   periods matrix is accepted for a single covariate).
#' @param N0 Number of control units.
#' @param T0 Number of pre-treatment periods.
#'
#' @return `Y` with the covariates' estimated contribution removed, with the
#'   coefficients attached as the `"beta"` attribute. A covariate that is
#'   collinear with the fixed effects gets an `NA` coefficient and is left
#'   out of the adjustment.
#'
#' @references Kranz, S. (2022). Synthetic Difference-in-Differences with
#'   Time-Varying Covariates. Working paper.
#'
#' @examples
#' \donttest{
#' data("california_prop99", package = "synthdid")
#' setup <- synthdid::panel.matrices(california_prop99)
#' set.seed(1)
#' X <- setup$Y * 0.01 + rnorm(length(setup$Y))
#' Yres <- project_out_covariates(setup$Y, X, setup$N0, setup$T0)
#' attr(Yres, "beta")
#' }
#'
#' @export
project_out_covariates <- function(Y, X, N0, T0) {
  if (length(dim(X)) == 2L) dim(X) <- c(dim(X), 1L)
  if (!identical(dim(X)[1:2], dim(Y))) {
    cli::cli_abort("{.arg X} is {dim(X)[1]} x {dim(X)[2]} per covariate but {.arg Y} is {nrow(Y)} x {ncol(Y)}.")
  }
  K <- dim(X)[3]
  grid <- expand.grid(i = seq_len(nrow(Y)), t = seq_len(ncol(Y)))
  cells <- cbind(grid$i, grid$t)
  df <- data.frame(
    y = Y[cells],
    unit = factor(grid$i), time = factor(grid$t)
  )
  xs <- paste0("x", seq_len(K))
  for (k in seq_len(K)) df[[xs[k]]] <- X[, , k][cells]

  untreated <- grid$i <= N0 | grid$t <= T0
  form <- stats::as.formula(paste("y ~", paste(xs, collapse = " + "), "+ unit + time"))
  coefs <- stats::coef(stats::lm(form, data = df[untreated, , drop = FALSE]))
  beta <- stats::setNames(unname(coefs[xs]), dimnames(X)[[3]] %||% xs)

  out <- Y
  for (k in seq_len(K)) if (!is.na(beta[[k]])) out <- out - beta[[k]] * X[, , k]
  attr(out, "beta") <- beta
  out
}
