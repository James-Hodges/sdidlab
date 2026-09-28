#' Permute treatment in space across every donor
#'
#' Refits [fit_scm()] once per unit, each time pretending that unit is the
#' treated one and that the rest of the pool are its donors. The result is
#' the reference distribution [rmspe_pvalue()] ranks against, and the set of
#' grey placebo paths behind an Abadie-style gap figure.
#'
#' The treated unit is always excluded from the donor pool of a placebo fit.
#' If it were left in, a placebo could reproduce itself out of the genuinely
#' treated unit's post-treatment trajectory, which would import the effect
#' being tested into the null distribution.
#'
#' Everything is deterministic: the simplex solve is Frank-Wolfe from a
#' fixed start and the bias correction is a regression, so no seed is
#' involved and repeated calls agree exactly. That is the difference between
#' this and the placebo variance in [fit_sdid()], which is simulated.
#'
#' @param Y Outcome matrix, units x periods, controls first and the treated
#'   unit last.
#' @param N0 Number of control units.
#' @param T0 Number of pre-treatment periods.
#' @param X Optional units x predictors matrix, rows aligned with `Y`.
#' @param ... Passed to [fit_scm()] -- `outcome_periods`, `bias_correct`,
#'   `v_method`, `bc_donors` and friends. The same arguments are used for the treated fit
#'   and every placebo, which is the point.
#'
#' @return A list with `results`, a tibble of one row per unit (`unit`,
#'   `treated`, `att`, `pre_rmspe`, `post_rmspe`, `rmspe_ratio` and their
#'   `_raw` counterparts computed before bias correction), `gaps`, a
#'   units x periods matrix of the corrected gap paths, and `gaps_raw`, the
#'   same before correction.
#'
#' @seealso [rmspe_pvalue()] for the ranking that turns this into a
#'   p-value, and [fit_scm()] for the fit it repeats.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(6 * 8), nrow = 6,
#'             dimnames = list(c(paste0("d", 1:5), "treated"), 2010:2017))
#' Y[6, 7:8] <- Y[6, 7:8] + 3
#' scm_placebo_permutation(Y, N0 = 5, T0 = 6)$results
#'
#' @export
scm_placebo_permutation <- function(Y, N0, T0, X = NULL, ...) {
  units <- rownames(Y)
  if (is.null(units)) {
    units <- paste0("unit", seq_len(nrow(Y)))
    rownames(Y) <- units
  }
  n <- nrow(Y)

  fits <- purrr::map(seq_len(n), function(i) {
    # Reorder so unit i is last (the treated slot) and everyone else is a
    # donor, preserving the row alignment of X. For a placebo the genuinely
    # treated unit (row n) is dropped from the donors -- see above.
    donors <- setdiff(seq_len(n), if (i == n) n else c(i, n))
    ord <- c(donors, i)
    # An exact placebo fit is reported through its NA ratio, which
    # rmspe_pvalue() then refuses to rank; warning once per placebo would
    # only bury that.
    fit <- withCallingHandlers(
      fit_scm(Y[ord, , drop = FALSE], N0 = length(donors), T0 = T0,
              X = if (is.null(X)) NULL else X[ord, , drop = FALSE], ...),
      sdidlab_exact_fit = function(w) invokeRestart("muffleWarning")
    )
    list(
      row = tibble::tibble(
        unit        = units[i],
        treated     = i == n,
        att         = fit$att,
        pre_rmspe   = fit$pre_rmspe,
        post_rmspe  = fit$post_rmspe,
        rmspe_ratio = fit$rmspe_ratio,
        pre_rmspe_raw   = fit$pre_rmspe_raw,
        post_rmspe_raw  = fit$post_rmspe_raw,
        rmspe_ratio_raw = fit$rmspe_ratio_raw,
        rmspe_ratio_hybrid = fit$rmspe_ratio_hybrid,
        exact_fit   = fit$exact_fit
      ),
      gap = fit$gap,
      gap_raw = fit$gap_raw
    )
  })

  stack <- function(part) {
    m <- do.call(rbind, purrr::map(fits, part))
    dimnames(m) <- list(units, colnames(Y))
    m
  }
  gaps <- stack("gap")
  gaps_raw <- stack("gap_raw")

  list(
    results = purrr::list_rbind(purrr::map(fits, "row")),
    gaps    = gaps,
    gaps_raw = gaps_raw
  )
}
