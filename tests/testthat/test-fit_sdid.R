skip_if_not_installed("synthdid")

# A panel with a known answer: controls on a common trend at their own
# levels, the treated unit parallel to them, then a jump in the last periods.
#
# The noise is not decoration. synthdid derives its regularisation parameter
# from the SD of control first differences, so a perfectly deterministic
# panel gives zeta = 0 and the Frank-Wolfe solver returns NaN.
toy_Y <- function(effect = 5, n_ctrl = 12, T_pre = 8, T_post = 3,
                  sigma = 0.4, seed = 42) {
  set.seed(seed)
  levels <- seq(50, by = 2, length.out = n_ctrl + 1)
  trend  <- seq_len(T_pre + T_post)
  Y <- outer(levels, trend, function(l, t) l + 0.5 * t)
  Y <- Y + matrix(stats::rnorm(length(Y), sd = sigma), nrow = nrow(Y))
  post <- (T_pre + 1):(T_pre + T_post)
  Y[n_ctrl + 1, post] <- Y[n_ctrl + 1, post] + effect
  rownames(Y) <- c(paste0("c", seq_len(n_ctrl)), "treated")
  colnames(Y) <- trend
  Y
}

test_that("a known effect is recovered", {
  Y <- toy_Y(effect = 5)
  r <- fit_sdid(Y, N0 = 12, T0 = 8, replications = 20)
  expect_equal(r$att, 5, tolerance = 0.5)
  expect_lt(r$ci_lo, 5)
  expect_gt(r$ci_hi, 5)
})

test_that("a zero effect returns approximately zero", {
  Y <- toy_Y(effect = 0)
  r <- fit_sdid(Y, N0 = 12, T0 = 8, replications = 20)
  expect_lt(abs(r$att), 0.6)
})

test_that("the confidence interval brackets the estimate and matches the SE", {
  r <- fit_sdid(toy_Y(), N0 = 12, T0 = 8, replications = 20)
  expect_lt(r$ci_lo, r$att)
  expect_gt(r$ci_hi, r$att)
  expect_equal(r$ci_hi - r$att, qnorm(0.975) * r$se, tolerance = 1e-8)
  expect_equal(r$t, r$att / r$se, tolerance = 1e-8)
})

test_that("trend_rmse ignores a constant level gap but level_gap records it", {
  # Shifting the treated unit up by a constant in EVERY period is exactly
  # what SDID's unit fixed effects absorb: the parallel-trends fit must not
  # degrade, and the effect must not change.
  Y <- toy_Y(effect = 5)
  shifted <- Y
  shifted["treated", ] <- shifted["treated", ] + 20
  base <- fit_sdid(Y, 12, 8, replications = 20)
  shft <- fit_sdid(shifted, 12, 8, replications = 20)
  expect_equal(shft$trend_rmse, base$trend_rmse, tolerance = 1e-6)
  expect_equal(shft$att, base$att, tolerance = 1e-6)
  expect_gt(shft$level_gap - base$level_gap, 19)
})

test_that("unit weights are named, non-negative and sum to one", {
  r <- fit_sdid(toy_Y(), N0 = 12, T0 = 8, replications = 20)
  expect_length(r$omega, 12)
  expect_equal(names(r$omega), paste0("c", 1:12))
  expect_true(all(r$omega >= 0))
  expect_equal(sum(r$omega), 1, tolerance = 1e-6)
})

test_that("degenerate panels abort with a specific message", {
  Y <- toy_Y()
  expect_error(fit_sdid(Y, N0 = 1, T0 = 8), "control")
  expect_error(fit_sdid(Y, N0 = 12, T0 = 1), "pre-period")
})

test_that("se_method = 'none' skips inference but keeps the point estimate", {
  Y <- toy_Y(effect = 5)
  full <- fit_sdid(Y, 12, 8, replications = 20)
  fast <- fit_sdid(Y, 12, 8, se_method = "none")
  expect_equal(fast$att, full$att, tolerance = 1e-8)
  expect_true(is.na(fast$se))
  expect_true(is.na(fast$t))
  expect_true(is.na(fast$ci_lo))
  # Diagnostics that do not depend on inference must still be present.
  expect_equal(fast$trend_rmse, full$trend_rmse, tolerance = 1e-8)
  expect_length(fast$omega, 12)
})
