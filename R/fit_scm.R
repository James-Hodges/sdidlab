#' Fit a bias-corrected synthetic control
#'
#' The Abadie, Diamond and Hainmueller synthetic control, matched on an
#' arbitrary predictor matrix rather than on pre-period outcomes alone, with
#' the Abadie-L'Hour / Ben-Michael-Feller-Rothstein bias correction applied
#' on top. This is the estimator Orlova, Rury and Wiltshire (2023) use as
#' their preferred specification, via the Stata package `allsynth`
#' (Wiltshire 2022); [fit_sdid()] covers the synthetic-DiD and plain-DiD
#' comparisons they report alongside it.
#'
#' The weight solve is checked against an exact quadratic-programming
#' solution in the test suite: at convergence it reaches the same objective,
#' which is at least as low as `Synth`'s interior-point solve.
#'
#' It is deliberately separate from [fit_sdid()] rather than another
#' `estimator =` branch there. `fit_sdid()` delegates to `synthdid`, whose
#' `sc_estimate()` matches on pre-period outcomes only and has no bias
#' correction; bolting predictor matching onto it would change the meaning
#' of the existing `"sc"` option and put every stored run's numbers in play.
#'
#' # What the two steps do
#'
#' **Matching.** Donor weights `w` sit on the simplex and minimise
#' \eqn{(X_1 - X_0 w)' V (X_1 - X_0 w)}, where the columns of `X` are the
#' predictors and `V` is diagonal. Predictors are standardised by their
#' cross-unit standard deviation first, so `V = I` means "weight every
#' predictor equally after standardisation" -- the special case of ADH that
#' does not require the nested optimisation. `v_method = "mspe"` instead
#' searches over diagonal `V` for the one minimising pre-period MSPE of the
#' outcome, which is what `Synth` does by default.
#'
#' **Neither search dominates, which is why `"synth"` runs both.** `Synth`
#' optimises `V` with `optimx` from a single equal-weight start and solves
#' for `w` with an interior-point QP; this package uses several
#' deterministic starts and an exact active-set solve (see
#' [sc_weights_simplex()]). The `V` surface
#' is rough, so the two searches can land in *different* local minima, and
#' either can find the better one. The objective is the same and is directly
#' comparable, so taking the better of the two is free accuracy and the
#' only defensible default when the surface is this rough.
#'
#' When several convex combinations of donors fit equally well the
#' minimiser is not unique, and different solvers can return very different
#' comparison groups with the same loss. `penalty` resolves such ties;
#' without it, report `omega` rather than assuming it.
#'
#' **Bias correction.** Even a good synthetic control leaves a residual
#' discrepancy in the predictors, and that discrepancy biases the estimate.
#' The correction regresses the outcome on the predictors *among the donor
#' units*, one regression per period, and subtracts the model's prediction
#' of the discrepancy. A ridge penalty keeps it defined when the positive-
#' weight donors are fewer than the predictors.
#'
#' Note the diagnostic consequence, which is a feature and not a bug: **in
#' any period whose outcome is itself a predictor, the corrected gap is
#' mechanically zero**, because the period's regression fits that column
#' exactly. Seeing a row of exact zeros in the gap path is confirmation that
#' the correction ran, not evidence of a perfect match. It is what produces
#' the flat segments at matched years in bias-corrected SCM gap figures.
#'
#' @param Y Outcome matrix, units x periods, controls first and the treated
#'   unit last -- the layout `synthdid::panel.matrices()` returns.
#' @param N0 Number of control units (the first `N0` rows of `Y`).
#' @param T0 Number of pre-treatment periods (the first `T0` columns).
#' @param X Optional units x predictors matrix of time-invariant predictors,
#'   rows in the same order as `Y`. Columns that are constant across units
#'   are dropped with a warning: they cannot discriminate between donors.
#' @param outcome_periods Integer indices of pre-treatment columns of `Y` to
#'   include as predictors, on top of `X`. `NULL` (the default) uses every
#'   pre-treatment period. Pass a subset to reproduce a specification that
#'   matches on selected pre-treatment years.
#' @param bias_correct Apply the bias correction. `TRUE` by default.
#' @param bc_ridge Ridge penalty for the bias-correction regressions, on
#'   standardised predictors. Small and positive; only binds when the
#'   positive-weight donors are few.
#' @param bc_donors Which donors the bias-correction regressions are fitted
#'   on. `"positive"` uses only donors carrying positive weight, as
#'   `allsynth` does; `"all"` uses the whole pool, as augmented SCM does;
#'   `"auto"` (the default) uses the positive-weight donors when there are
#'   comfortably more of them than predictors and falls back to the whole
#'   pool otherwise, reporting when it does. The fallback matters: a sparse
#'   synthetic control can easily name fewer donors than there are
#'   predictors, and the regression is then rank-deficient exactly when the
#'   correction is being asked to do the most work.
#' @param v_method How the diagonal of `V` is chosen. `"equal"` (default)
#'   weights standardised predictors equally. `"mspe"` searches for the `V`
#'   minimising pre-period MSPE, from `v_starts` deterministic starting
#'   points. `"synth"` does the `"mspe"` search AND asks the `Synth` package
#'   for its own solution, then keeps whichever fits the pre-period better --
#'   see the note below on why that is not redundant. `"synth"` needs the
#'   suggested `Synth` package and falls back to `"mspe"` with a warning if
#'   it is missing.
#' @param v_starts Restarts for the `"mspe"` search, the first from equal
#'   weights and the rest scattered. One start is not enough -- see the
#'   note in the code. Ignored when `v_method = "equal"`.
#' @param v_maxit Nelder-Mead iterations per restart.
#' @param penalty Abadie-L'Hour (2021) penalty `lambda`, non-negative.
#'   Zero (the default) is the classic synthetic control. A positive value
#'   adds a charge for each donor's own distance from the treated unit,
#'   which makes the weights **unique** and favours donors that resemble the
#'   treated unit individually rather than only on average. Its main use
#'   here is breaking ties when the pre-period is fitted exactly: any small
#'   positive value then picks, among the exact fits, the one built from the
#'   nearest donors. Larger values trade pre-period fit for locality, and as
#'   `penalty` grows the synthetic control becomes the single nearest
#'   neighbour. The scale is that of the standardised predictors.
#'
#' @return A list with `att` (mean corrected gap over the post-treatment
#'   periods), `gap` (the corrected gap in every period, named by the
#'   columns of `Y`), `gap_raw` (before correction), `pre_rmspe`,
#'   `post_rmspe`, `rmspe_ratio` (`post_rmspe / pre_rmspe`), the same three
#'   computed on the uncorrected gaps with a `_raw` suffix --
#'   plus `rmspe_ratio_hybrid` (corrected post over uncorrected pre) -- see
#'   the code for which to hand [rmspe_pvalue()] -- `omega` (donor weights, named by the control
#'   rownames), `v` (the predictor weights used), `n_predictors`,
#'   `bias_correct`, `N0` and `T0`, plus `exact_fit`: `TRUE` when the
#'   uncorrected synthetic control reproduces the pre-period exactly. An
#'   exact fit warns (class `sdidlab_exact_fit`), because the weights are
#'   then generally not unique and the estimate depends on which optimal
#'   weights the solver returns; the RMSPE ratios are `NA`.
#'
#' @references
#' Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic Control
#' Methods for Comparative Case Studies. *Journal of the American Statistical
#' Association* 105(490), 493-505.
#'
#' Abadie, A. and L'Hour, J. (2021). A Penalized Synthetic Control Estimator
#' for Disaggregated Data. *Journal of the American Statistical Association*
#' 116(536), 1817-1834.
#'
#' Ben-Michael, E., Feller, A. and Rothstein, J. (2021). The Augmented
#' Synthetic Control Method. *Journal of the American Statistical
#' Association* 116(536), 1789-1803.
#'
#' @seealso [rmspe_pvalue()] for the permutation inference this feeds,
#'   [fit_sdid()] for the synthetic-DiD and DiD comparisons.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(11 * 8), nrow = 11,
#'             dimnames = list(c(paste0("d", 1:10), "treated"), 2010:2017))
#' Y[11, 7:8] <- Y[11, 7:8] + 3
#' fit_scm(Y, N0 = 10, T0 = 6)$att
#'
#' @export
fit_scm <- function(Y, N0, T0, X = NULL, outcome_periods = NULL,
                    bias_correct = TRUE, bc_ridge = 1e-4,
                    bc_donors = c("auto", "positive", "all"),
                    v_method = c("equal", "mspe", "synth"),
                    v_starts = 8L, v_maxit = 1000L, penalty = 0) {
  v_method  <- match.arg(v_method)
  bc_donors <- match.arg(bc_donors)
  if (N0 < 2) cli::cli_abort("Need at least 2 control units; got {N0}.")
  if (T0 < 2) cli::cli_abort("Need at least 2 pre-periods; got {T0}.")
  if (ncol(Y) <= T0) cli::cli_abort("No post-treatment periods: T0 = {T0}, ncol(Y) = {ncol(Y)}.")

  pre  <- seq_len(T0)
  post <- seq.int(T0 + 1L, ncol(Y))
  keep <- if (is.null(outcome_periods)) pre else intersect(outcome_periods, pre)

  # Predictors: the requested pre-period outcomes, then anything in X.
  Z <- Y[, keep, drop = FALSE]
  colnames(Z) <- paste0("y_", colnames(Y)[keep])
  if (!is.null(X)) {
    if (nrow(X) != nrow(Y)) {
      cli::cli_abort("X has {nrow(X)} rows but Y has {nrow(Y)}.")
    }
    Z <- cbind(Z, as.matrix(X))
  }

  # Standardise by cross-unit sd. A predictor with no variation carries no
  # information about which donors to pick and would divide by zero.
  sds <- apply(Z, 2, stats::sd)
  flat <- !is.finite(sds) | sds < .Machine$double.eps^0.5
  if (any(flat)) {
    cli::cli_warn(
      "Dropping {sum(flat)} constant predictor{?s}: {.field {colnames(Z)[flat]}}."
    )
    Z <- Z[, !flat, drop = FALSE]
    sds <- sds[!flat]
  }
  if (!ncol(Z)) cli::cli_abort("No usable predictors left after dropping constants.")
  Z <- sweep(Z, 2, sds, "/")

  Z0 <- Z[seq_len(N0), , drop = FALSE]
  Z1 <- Z[N0 + 1L, ]
  Y0 <- Y[seq_len(N0), , drop = FALSE]
  Y1 <- Y[N0 + 1L, ]

  mspe_for <- function(v) {
    w <- sc_weights_simplex(Z0, Z1, v, penalty = penalty)
    mean((Y1[pre] - as.numeric(crossprod(Y0[, pre, drop = FALSE], w)))^2)
  }
  if (identical(v_method, "synth") && !requireNamespace("Synth", quietly = TRUE)) {
    cli::cli_warn(c(
      'Package {.pkg Synth} is not installed; falling back to {.val mspe}.',
      i = '{.code install.packages("Synth")}'
    ))
    v_method <- "mspe"
  }

  v <- rep(1, ncol(Z))
  if (v_method %in% c("mspe", "synth") && ncol(Z) > 1L) {
    # The pre-period MSPE is a rough function of V with plenty of local
    # minima, and a single Nelder-Mead run from equal weights lands in one --
    # it can put nearly all the weight on a single covariate and essentially
    # none on the pre-period outcomes, which is the V-degeneracy the
    # synthetic control literature warns about. Restarting from scattered
    # points and keeping the best guards against it. Optimise on the log scale so
    # weights stay positive; only relative sizes matter to the solve.
    # Starts are DETERMINISTIC -- a van der Corput low-discrepancy sequence
    # stretched over [-3, 3] in log space, not draws. Nothing here may touch
    # the RNG: fit_scm() is called inside placebo loops that also drive the
    # simulated standard errors elsewhere in this package, and consuming
    # randomness would silently shift their stream.
    starts <- c(
      list(rep(0, ncol(Z))),
      lapply(seq_len(max(0L, v_starts - 1L)), function(i) {
        vapply(seq_len(ncol(Z)),
               function(k) 6 * van_der_corput(i, base = k + 1L) - 3,
               numeric(1))
      })
    )
    fits <- purrr::compact(lapply(starts, function(s) {
      tryCatch(
        stats::optim(s, function(lv) mspe_for(exp(lv)), method = "Nelder-Mead",
                     control = list(maxit = v_maxit, reltol = 1e-10)),
        error = function(e) NULL
      )
    }))
    if (length(fits)) {
      v <- exp(fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]$par)
    }

    if (identical(v_method, "synth")) {
      # Synth scales predictors by their cross-unit sd exactly as this
      # function does (up to a constant common to every predictor, which
      # cancels), so its solution.v transfers directly onto our own
      # standardised Z and the two are comparable on the same objective.
      v_syn <- tryCatch({
        fit <- Synth::synth(
          X1 = matrix(Z1, ncol = 1), X0 = t(Z0),
          Z1 = matrix(Y1[pre], ncol = 1), Z0 = t(Y0[, pre, drop = FALSE]),
          verbose = FALSE
        )
        as.numeric(fit$solution.v)
      }, error = function(e) {
        cli::cli_warn("Synth failed ({conditionMessage(e)}); keeping the mspe search.")
        NULL
      })
      if (length(v_syn) == ncol(Z) && all(is.finite(v_syn)) && sum(v_syn) > 0) {
        if (mspe_for(v_syn) < mspe_for(v)) v <- v_syn
      }
    }
  }
  v <- v / sum(v) * length(v)
  names(v) <- colnames(Z)

  w <- sc_weights_simplex(Z0, Z1, v, penalty = penalty)
  names(w) <- rownames(Y)[seq_len(N0)]

  gap_raw <- Y1 - as.numeric(crossprod(Y0, w))
  names(gap_raw) <- colnames(Y)

  gap <- gap_raw
  if (bias_correct) {
    # Regress the outcome on the predictors among donors, one period at a
    # time, and remove the model's read on the residual predictor
    # discrepancy. Which donors that regression sees is `bc_donors`.
    positive <- which(w > 1e-8)
    enough   <- length(positive) >= ncol(Z) + 2L
    use <- switch(bc_donors,
      all      = seq_len(N0),
      positive = positive,
      auto     = if (enough) positive else seq_len(N0)
    )
    if (identical(bc_donors, "auto") && !enough && length(positive) < N0) {
      cli::cli_inform(c(
        "Bias correction fitted on all {N0} donors, not the {length(positive)} with positive weight.",
        i = "{ncol(Z)} predictors need more donors than that to be identified."
      ))
    }
    if (length(use) >= 2L) {
      Zu <- cbind(1, Z0[use, , drop = FALSE])
      wu <- w[use]
      pen <- diag(c(0, rep(bc_ridge, ncol(Z))))   # never penalise the intercept
      xtx <- crossprod(Zu) + pen
      beta <- tryCatch(solve(xtx, crossprod(Zu, Y0[use, , drop = FALSE])),
                       error = function(e) NULL)
      if (is.null(beta)) {
        cli::cli_warn("Bias-correction regression is singular; reporting uncorrected gaps.")
      } else {
        # Discrepancy between the treated unit's predictors and the
        # weighted donors', mapped through the outcome model.
        disc <- c(0, Z1 - as.numeric(crossprod(Z0[use, , drop = FALSE], wu / sum(wu))))
        gap <- gap_raw - as.numeric(crossprod(beta, disc))
        names(gap) <- colnames(Y)
      }
    } else {
      cli::cli_warn(
        "Only {length(use)} donor{?s} carry positive weight; reporting uncorrected gaps."
      )
    }
  }

  pre_rmspe  <- sqrt(mean(gap[pre]^2))
  post_rmspe <- sqrt(mean(gap[post]^2))
  pre_raw    <- sqrt(mean(gap_raw[pre]^2))
  post_raw   <- sqrt(mean(gap_raw[post]^2))

  # An RMSPE this small relative to the outcome's scale is an exact fit up
  # to floating point, not a very good one. Treating 1e-9 as "positive"
  # makes the post/pre ratio astronomically large and ranks the unit first
  # in any permutation test by construction, so exact fits get an
  # undefined ratio instead -- which rmspe_pvalue() refuses to rank.
  exact_tol <- 1e-7 * (1 + sqrt(mean(Y1[pre]^2)))
  ratio <- function(post, pre) if (pre > exact_tol) post / pre else NA_real_
  exact_fit <- pre_raw <= exact_tol
  # With a positive penalty the weights are unique, so an exact fit is no
  # longer ambiguous -- only the ratio stays undefined.
  if (exact_fit && penalty == 0) {
    cli::cli_warn(c(
      "The synthetic control reproduces the pre-period exactly (RMSPE {signif(pre_raw, 2)}).",
      i = "With an exact fit the weights are generally not unique, and different optimal weights imply different effects.",
      i = "Treat the estimate as unidentified, or pass a small {.arg penalty} to select the unique penalised solution (Abadie and L'Hour 2021)."
    ), class = "sdidlab_exact_fit")
  }

  list(
    att          = mean(gap[post]),
    gap          = gap,
    gap_raw      = gap_raw,
    pre_rmspe    = pre_rmspe,
    post_rmspe   = post_rmspe,
    pre_rmspe_raw  = pre_raw,
    post_rmspe_raw = post_raw,
    # Zero pre-period RMSPE means the predictors pinned every pre-period
    # exactly; the ratio is undefined rather than infinite.
    rmspe_ratio  = ratio(post_rmspe, pre_rmspe),
    # The ratio on UNCORRECTED gaps, which is the one to rank on. The
    # correction pins the gap to zero in every period used as a predictor,
    # so a corrected pre-period RMSPE is computed over a mixture of real
    # residuals and structural zeros -- and it is a mixture whose
    # composition differs between the treated unit and a placebo only by
    # accident. It also inherits the correction's linear extrapolation,
    # which on a large predictor discrepancy can be worse than no
    # correction at all. De-biasing the reported effect and measuring
    # pre-period fit are different jobs.
    rmspe_ratio_raw = ratio(post_raw, pre_raw),
    # Corrected post-period gaps over the UNCORRECTED pre-period fit. This
    # is the ranking statistic to prefer when the predictor set pins most of
    # the pre-period: when most pre-period outcomes are predictors the
    # correction fits them exactly, the corrected pre-period RMSPE collapses
    # towards zero, and `rmspe_ratio` explodes, which ranks noise.
    # Using the uncorrected pre-fit keeps the denominator a measure of MATCH
    # QUALITY, which is what it is for in Abadie's ratio, while the
    # numerator still reflects the de-biased effect.
    rmspe_ratio_hybrid = ratio(post_rmspe, pre_raw),
    exact_fit    = exact_fit,
    omega        = w,
    v            = v,
    n_predictors = ncol(Z),
    bias_correct = bias_correct,
    N0           = N0,
    T0           = T0
  )
}

#' Simplex weights minimising a V-weighted predictor discrepancy
#'
#' Solves the synthetic-control weight problem exactly:
#' \deqn{\min_w (Z_1 - Z_0' w)' V (Z_1 - Z_0' w) + \lambda \sum_j w_j
#'   (Z_1 - Z_{0j})' V (Z_1 - Z_{0j})}
#' over the unit simplex (`w >= 0`, `sum(w) = 1`), where `V = diag(v)` and
#' `lambda` is the Abadie-L'Hour `penalty` (zero by default, the classic
#' synthetic control).
#'
#' # Method
#'
#' A primal active-set method -- Lawson and Hanson's non-negative least
#' squares algorithm, adapted to the simplex. It keeps a working set of
#' donors, solves the equality-constrained problem on that set exactly
#' (through a pseudo-inverse, so a singular set is handled), steps back to
#' feasibility when a weight would turn negative, and adds the donor whose
#' gradient most violates optimality. It stops at the exact KKT point: every
#' donor with positive weight has the same gradient and no donor has a
#' lower one. Starting from the nearest donor, it typically needs about two
#' steps per donor in the final support.
#'
#' Two alternatives were tried and rejected. A first-order (Frank-Wolfe)
#' solve needs over 12,000 iterations on the highly collinear pre-period
#' outcomes of the Proposition 99 panel, and at any practical iteration cap
#' returns weights well short of the optimum. A dual QP solver needs a ridge
#' to make the singular quadratic term invertible, which with a penalty
#' makes the problem nearly linear and the solve numerically unreliable.
#'
#' # Uniqueness
#'
#' When the treated unit lies inside the donors' convex hull, many weight
#' vectors can fit the predictors equally well, and they imply different
#' counterfactuals. Which of them is returned then depends on the order in
#' which donors enter the working set -- deterministic, but arbitrary.
#' [fit_scm()] warns when the fit is exact. A positive `penalty` makes the
#' solution unique (Abadie and L'Hour 2021).
#'
#' @param Z0 Donors x predictors matrix, already standardised.
#' @param Z1 The treated unit's predictor vector.
#' @param v Non-negative predictor weights (the diagonal of `V`).
#' @param penalty The Abadie-L'Hour penalty `lambda`, non-negative. The
#'   scale is that of the standardised predictors. See [fit_scm()] for when
#'   to use it.
#'
#' @return A numeric vector of donor weights, non-negative and summing to 1.
#'
#' @references Lawson, C. L. and Hanson, R. J. (1974). *Solving Least
#'   Squares Problems*. Prentice-Hall, chapter 23.
#'
#' @examples
#' Z0 <- matrix(c(0, 1, 2, 0, 1, 2), nrow = 3)
#' sc_weights_simplex(Z0, c(1, 1), v = c(1, 1))
#'
#' @export
sc_weights_simplex <- function(Z0, Z1, v = rep(1, ncol(Z0)), penalty = 0) {
  if (!is.numeric(penalty) || length(penalty) != 1L || !is.finite(penalty) || penalty < 0) {
    cli::cli_abort("{.arg penalty} must be a single non-negative number.")
  }
  Z0 <- as.matrix(Z0)
  n <- nrow(Z0)
  if (n == 1L) return(1)
  if (length(Z1) != ncol(Z0) || length(v) != ncol(Z0)) {
    cli::cli_abort("{.arg Z1} and {.arg v} must have one entry per predictor ({ncol(Z0)}).")
  }

  # Objective w'Gw - 2c'w + penalty * dist'w (+ constant), with
  #   G = Z0 V Z0',  c = Z0 V Z1,  dist_j = (Z1 - Z0j)' V (Z1 - Z0j).
  vZ0t <- t(Z0) * v
  G    <- Z0 %*% vZ0t
  cc   <- as.numeric(Z0 %*% (v * Z1))
  dist <- as.numeric((sweep(Z0, 2, Z1, "-")^2) %*% v)
  lin  <- 2 * cc - penalty * dist            # the gradient is 2 G w - lin
  tol  <- 1e-11 * max(1, abs(lin), max(abs(G)))

  # Null space of the working set's quadratic term within the simplex's
  # tangent space: directions p with G_SS p = 0 and sum(p) = 0. Along them
  # the objective is linear, so if `lin` has a component there the
  # restricted problem has no interior minimum.
  null_dirs <- function(S) {
    Q <- rbind(G[S, S, drop = FALSE], 1)
    sv <- svd(Q, nu = 0, nv = length(S))
    d <- c(sv$d, rep(0, length(S) - length(sv$d)))
    sv$v[, d <= max(sv$d) * 1e-10, drop = FALSE]
  }
  # Stationary point of the working-set problem: the KKT system
  #   [2 G_SS  1] [w]   [lin_S]
  #   [1'      0] [nu] = [1    ]
  # by pseudo-inverse (minimum norm when the system is singular).
  stationary <- function(S) {
    m <- length(S)
    M <- rbind(cbind(2 * G[S, S, drop = FALSE], 1), c(rep(1, m), 0))
    sv <- svd(M)
    keep <- sv$d > max(sv$d) * 1e-12
    z <- sv$v[, keep, drop = FALSE] %*% ((t(sv$u[, keep, drop = FALSE]) %*% c(lin[S], 1)) / sv$d[keep])
    z[seq_len(m)]
  }
  # Move w_S towards `target` (or along direction `p`) until the first
  # weight reaches zero, then drop the zeros from the working set.
  drop_zeros <- function(S) S[w[S] > tol]

  w <- numeric(n)
  S <- which.min(dist)                       # start at the nearest donor
  w[S] <- 1
  for (outer in seq_len(20L * n + 100L)) {
    # Solve the working-set problem to optimality, shrinking S as needed.
    for (inner in seq_len(4L * n + 10L)) {
      N <- if (length(S) > 1L) null_dirs(S) else NULL
      if (!is.null(N) && ncol(N) && sqrt(sum((crossprod(N, lin[S]))^2)) > tol) {
        # Ray step: objective decreases linearly along p = N N' lin_S.
        p <- as.numeric(N %*% crossprod(N, lin[S]))
        neg <- p < 0
        w[S] <- w[S] + min(w[S][neg] / -p[neg]) * p
        w[w < tol] <- 0
        S <- drop_zeros(S)
        next
      }
      z <- stationary(S)
      if (all(z > tol)) {
        w[] <- 0
        w[S] <- z
        break
      }
      # Step from the current feasible point towards z only as far as
      # feasibility allows, then drop whatever reached zero.
      neg <- z <= tol
      alpha <- min(w[S][neg] / (w[S][neg] - z[neg]))
      w[S] <- w[S] + alpha * (z - w[S])
      w[w < tol] <- 0
      S <- drop_zeros(S)
    }
    grad <- 2 * as.numeric(G %*% w) - lin
    mu <- mean(grad[S])                      # common gradient on the support
    out <- setdiff(seq_len(n), S)
    if (!length(out)) break
    j <- out[which.min(grad[out])]
    if (grad[j] >= mu - tol) break           # KKT: no donor outside S does better
    S <- c(S, j)
  }
  w[w < 0] <- 0
  w / sum(w)
}

#' A deterministic low-discrepancy point in [0, 1)
#'
#' The van der Corput sequence: write `i` in the given base and mirror its
#' digits about the radix point. Successive points fill the interval far
#' more evenly than independent draws, which is what makes it a good source
#' of restart positions for an optimiser -- and, unlike draws, it consumes
#' no randomness and returns the same points on every machine.
#'
#' @param i Index, a non-negative integer.
#' @param base Radix, an integer of at least 2. Using a different base per
#'   coordinate is what keeps a multi-dimensional start from lying on a
#'   diagonal.
#'
#' @return A single number in `[0, 1)`.
#'
#' @examples
#' vapply(1:4, van_der_corput, numeric(1), base = 2)
#'
#' @export
van_der_corput <- function(i, base = 2L) {
  if (base < 2L) cli::cli_abort("{.arg base} must be at least 2; got {base}.")
  out <- 0
  f <- 1 / base
  while (i > 0) {
    out <- out + f * (i %% base)
    i <- i %/% base
    f <- f / base
  }
  out
}
