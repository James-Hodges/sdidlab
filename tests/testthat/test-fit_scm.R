test_that("sc_weights_simplex returns a point on the simplex", {
  Z0 <- matrix(c(0, 1, 2, 5, 0, 1, 2, 5), nrow = 4)
  w <- sc_weights_simplex(Z0, c(1, 1), v = c(1, 1))
  expect_equal(sum(w), 1)
  expect_true(all(w >= 0))
})

# The synthetic-control objective has no unique argmin whenever several
# convex combinations hit the target exactly, which is the normal case with
# more donors than predictors. These tests therefore check the OBJECTIVE,
# never a particular weight vector.
objective <- function(Z0, Z1, w, v = rep(1, ncol(Z0))) {
  sum(v * (Z1 - as.numeric(crossprod(Z0, w)))^2)
}

test_that("sc_weights_simplex reaches an exactly attainable target", {
  Z0 <- matrix(c(-5, 0, 5, 9,
                 -5, 0, 5, 9), nrow = 4)
  w <- sc_weights_simplex(Z0, c(0, 0), v = c(1, 1))
  expect_lt(objective(Z0, c(0, 0), w), 1e-4)
})

test_that("sc_weights_simplex mixes donors when no single one will do", {
  Z0 <- matrix(c(0, 10, 100), nrow = 3)
  w <- sc_weights_simplex(Z0, 5, v = 1)
  expect_lt(objective(Z0, 5, w), 1e-4)
  expect_gt(sum(w > 1e-6), 1L)          # genuinely interior, not a corner
})

test_that("sc_weights_simplex gets as close as it can to an unattainable target", {
  # The target sits outside the donors' convex hull: the best it can do is
  # put all its weight on the nearest donor, and that solution IS unique.
  Z0 <- matrix(c(0, 10, 100), nrow = 3)
  w <- sc_weights_simplex(Z0, -50, v = 1)
  expect_equal(unname(w[1]), 1, tolerance = 1e-3)
})

test_that("fit_scm recovers a known effect on a parallel-trends panel", {
  fit <- fit_scm(toy(effect = 5), N0 = 5, T0 = 6, bias_correct = FALSE)
  expect_equal(fit$att, 5, tolerance = 0.5)
  expect_equal(sum(fit$omega), 1)
  expect_length(fit$gap, 8)
})

test_that("a period used as a predictor gets a mechanically zero gap", {
  # Columns 3 and 4 are 2012 and 2013 -- the years the replicated paper
  # matches on, and the reason its figures show flat segments there. With
  # no ridge the outcome model fits those columns exactly, so the
  # correction cancels the raw gap to machine precision. It needs enough
  # donors to be identified, which on a five-donor toy means all of them.
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                 bias_correct = TRUE, bc_ridge = 0, bc_donors = "all")
  expect_equal(unname(fit$gap[3]), 0, tolerance = 1e-8)
  expect_equal(unname(fit$gap[4]), 0, tolerance = 1e-8)
  # Uncorrected gaps in those periods are generally NOT zero.
  expect_gt(abs(fit$gap_raw[3]), 1e-3)
})

test_that("the ridge penalty leaves the zeroing near-exact, not exact", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                 bias_correct = TRUE, bc_ridge = 1e-4, bc_donors = "all")
  expect_lt(abs(fit$gap[3]), 1e-2)
  expect_gt(abs(fit$gap[3]), 0)
})

test_that("bias correction leaves the raw gaps available for comparison", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                 bc_donors = "all")
  expect_false(isTRUE(all.equal(unname(fit$gap), unname(fit$gap_raw))))
  expect_true(fit$bias_correct)
})

test_that("constant predictors are dropped with a warning, not divided by zero", {
  Y <- toy()
  X <- cbind(useful = c(1, 2, 3, 4, 5, 3), flat = rep(7, 6))
  expect_warning(fit <- fit_scm(Y, N0 = 5, T0 = 6, X = X), "constant predictor")
  expect_false("flat" %in% names(fit$v))
  expect_true("useful" %in% names(fit$v))
})

test_that("the RMSPE ratio is the post/pre ratio of the corrected gaps", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, bias_correct = FALSE)
  expect_equal(fit$rmspe_ratio, fit$post_rmspe / fit$pre_rmspe)
  expect_gt(fit$rmspe_ratio, 1)
})

test_that("fit_scm rejects panels it cannot estimate on", {
  Y <- toy()
  expect_error(fit_scm(Y, N0 = 1, T0 = 6), "at least 2 control")
  expect_error(fit_scm(Y, N0 = 5, T0 = 1), "at least 2 pre-periods")
  expect_error(fit_scm(Y, N0 = 5, T0 = 8), "No post-treatment periods")
})

test_that("the mspe V search does no worse than equal weights", {
  Y <- toy()
  X <- cbind(noise = c(9, -3, 4, 0, 7, 1), signal = c(10, 12, 14, 16, 18, 13))
  equal <- fit_scm(Y, N0 = 5, T0 = 6, X = X, v_method = "equal",
                   bias_correct = FALSE)
  mspe  <- fit_scm(Y, N0 = 5, T0 = 6, X = X, v_method = "mspe",
                   bias_correct = FALSE)
  expect_lte(mspe$pre_rmspe, equal$pre_rmspe + 1e-6)
})


test_that("auto falls back to the full pool when the support is too small", {
  # Two predictors against a sparse synthetic control: "positive" would fit
  # the regression on a handful of donors, "auto" says so and uses them all.
  expect_message(
    auto <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                    bc_donors = "auto"),
    "not the"
  )
  all_donors <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                        bc_donors = "all")
  expect_equal(auto$gap, all_donors$gap)
})

test_that("van_der_corput mirrors the digits of its index", {
  # 1, 2, 3 in base 2 are 0.1, 0.01, 0.11 reflected: 1/2, 1/4, 3/4.
  expect_equal(vapply(1:3, van_der_corput, numeric(1), base = 2),
               c(0.5, 0.25, 0.75))
  expect_equal(van_der_corput(0), 0)
  expect_equal(van_der_corput(1, base = 3), 1 / 3)
  expect_error(van_der_corput(1, base = 1), "at least 2")
})

test_that("the mspe restarts consume no randomness", {
  # fit_scm() runs inside placebo loops that share the RNG stream with the
  # simulated standard errors elsewhere in this package.
  set.seed(99)
  before <- .Random.seed
  fit_scm(toy(), N0 = 5, T0 = 6, v_method = "mspe", v_starts = 4L,
          v_maxit = 50L)
  expect_identical(.Random.seed, before)
})

test_that("the raw RMSPE ratio ignores the correction's structural zeros", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                 bc_donors = "all", bc_ridge = 0)
  # Two of six corrected pre-periods are pinned to zero by construction,
  # and the remaining four carry the correction's extrapolation -- which can
  # make them LARGER than the raw residuals, not smaller. That is precisely
  # why the corrected pre-period RMSPE is not a fit statistic.
  expect_equal(unname(fit$gap[3]), 0, tolerance = 1e-8)
  expect_equal(unname(fit$gap[4]), 0, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(fit$pre_rmspe, fit$pre_rmspe_raw)))
  # The raw ratio is computed entirely from gap_raw.
  expect_equal(fit$pre_rmspe_raw, sqrt(mean(fit$gap_raw[1:6]^2)))
  expect_equal(fit$post_rmspe_raw, sqrt(mean(fit$gap_raw[7:8]^2)))
  expect_equal(fit$rmspe_ratio_raw, fit$post_rmspe_raw / fit$pre_rmspe_raw)
})

test_that("the hybrid ratio pairs corrected post gaps with the uncorrected pre-fit", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, outcome_periods = c(3L, 4L),
                 bc_donors = "all")
  expect_equal(fit$rmspe_ratio_hybrid, fit$post_rmspe / fit$pre_rmspe_raw)
  # It is finite where the fully-corrected ratio can degenerate.
  expect_true(is.finite(fit$rmspe_ratio_hybrid))
})

test_that("v_method = 'synth' is never worse than the mspe search", {
  skip_if_not_installed("Synth")
  Y <- toy()
  X <- cbind(a = c(1, 2, 3, 4, 5, 3), b = c(9, -3, 4, 0, 7, 1))
  mspe_of <- function(f) {
    mean((Y[6, 1:6] - crossprod(Y[1:5, 1:6, drop = FALSE], f$omega))^2)
  }
  a <- fit_scm(Y, N0 = 5, T0 = 6, X = X, v_method = "mspe",
               v_starts = 4L, v_maxit = 200L, bc_donors = "all")
  b <- fit_scm(Y, N0 = 5, T0 = 6, X = X, v_method = "synth",
               v_starts = 4L, v_maxit = 200L, bc_donors = "all")
  # "synth" keeps whichever V fits better, so it can only improve.
  expect_lte(mspe_of(b), mspe_of(a) + 1e-8)
})

test_that("v_method = 'synth' degrades to 'mspe' when Synth is absent", {
  # Simulated by asking for a predictor set Synth cannot handle, which the
  # tryCatch must swallow rather than propagate.
  Y <- toy()
  expect_no_error(
    fit_scm(Y, N0 = 5, T0 = 6, v_method = "synth", v_starts = 2L,
            v_maxit = 50L, bc_donors = "all")
  )
})

test_that("bias-corrected gaps match an independent Abadie-L'Hour computation", {
  # mu_t fitted on the donors by OLS in each period; the corrected gap is the
  # raw gap minus mu_t(Z_treated) - sum_j w_j mu_t(Z_j).
  data("california_prop99", package = "synthdid", envir = environment())
  s <- synthdid::panel.matrices(california_prop99)
  fit <- fit_scm(s$Y, s$N0, s$T0, bias_correct = TRUE, bc_donors = "all",
                 v_method = "equal", bc_ridge = 0)
  Z <- s$Y[, 1:s$T0]; Z <- sweep(Z, 2, apply(Z, 2, stats::sd), "/")
  Xd <- cbind(1, Z[1:s$N0, ]); w <- fit$omega
  corr <- vapply(seq_len(ncol(s$Y)), function(t) {
    b <- qr.coef(qr(Xd), s$Y[1:s$N0, t]); b[is.na(b)] <- 0
    sum(b * c(1, Z[s$N0 + 1, ])) - sum(w * as.numeric(Xd %*% b))
  }, numeric(1))
  expect_equal(unname(fit$gap), unname(fit$gap_raw - corr), tolerance = 1e-6)
})

test_that("an exact pre-period fit warns and leaves the RMSPE ratio undefined", {
  # The treated unit is an exact mixture of two donors in every period, so
  # the pre-period fits to floating point. A 1e-9 denominator must not
  # produce an astronomically large ratio that ranks first by construction.
  Y <- toy()
  Y["treated", ] <- 0.5 * Y["d1", ] + 0.5 * Y["d3", ]
  Y["treated", 7:8] <- Y["treated", 7:8] + 5
  expect_warning(fit <- fit_scm(Y, N0 = 5, T0 = 6, bias_correct = FALSE),
                 class = "sdidlab_exact_fit")
  expect_true(fit$exact_fit)
  expect_true(is.na(fit$rmspe_ratio_raw))
  expect_true(is.na(fit$rmspe_ratio_hybrid))
  expect_equal(fit$att, 5, tolerance = 1e-6)
})

test_that("an ordinary fit is not flagged as exact", {
  fit <- fit_scm(toy(), N0 = 5, T0 = 6, bias_correct = FALSE)
  expect_false(fit$exact_fit)
  expect_true(is.finite(fit$rmspe_ratio_raw))
})

test_that("a permutation with an exact treated fit cannot be ranked", {
  Y <- toy()
  Y["treated", ] <- 0.5 * Y["d1", ] + 0.5 * Y["d3", ]
  perm <- scm_placebo_permutation(Y, N0 = 5, T0 = 6, bias_correct = FALSE)
  r <- perm$results
  expect_true(r$exact_fit[r$treated])
  expect_error(rmspe_pvalue(stats::setNames(r$rmspe_ratio_raw, r$unit), r$att,
                            r$unit[r$treated]), "non-finite")
})

test_that("an exact predictor match is flagged even when the pre-period is not reproduced", {
  data("california_prop99", package = "synthdid", envir = environment())
  s <- synthdid::panel.matrices(california_prop99)
  expect_warning(fit <- fit_scm(s$Y, s$N0, s$T0, outcome_periods = c(5, 12, 19), bias_correct = FALSE),
                 class = "sdidlab_exact_fit")
  expect_true(fit$exact_match)
  expect_false(fit$exact_fit)
  expect_true(is.finite(fit$rmspe_ratio_raw))
  # The penalty makes the weights unique and silences the warning.
  expect_no_warning(fit_scm(s$Y, s$N0, s$T0, outcome_periods = c(5, 12, 19), bias_correct = FALSE,
                            penalty = 1e-3))
})
