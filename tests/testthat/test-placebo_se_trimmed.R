skip_if_not_installed("synthdid")

# Heterogeneous controls: two thirds track a common trend closely, one third
# is deliberately erratic. The erratic ones are what inflate an untrimmed
# placebo distribution, so trimming should bite here and not elsewhere.
mixed_Y <- function(n_good = 16, n_bad = 8, T_pre = 8, T_post = 3,
                    effect = 3, seed = 7) {
  set.seed(seed)
  n_ctrl <- n_good + n_bad
  trend  <- seq_len(T_pre + T_post)
  lev    <- seq(50, by = 1.5, length.out = n_ctrl + 1)
  Y <- outer(lev, trend, function(l, t) l + 0.5 * t)
  Y <- Y + matrix(stats::rnorm(length(Y), sd = 0.3), nrow = nrow(Y))
  bad <- (n_good + 1):n_ctrl
  Y[bad, ] <- Y[bad, ] + matrix(stats::rnorm(length(bad) * ncol(Y), sd = 6),
                                nrow = length(bad))
  post <- (T_pre + 1):(T_pre + T_post)
  Y[n_ctrl + 1, post] <- Y[n_ctrl + 1, post] + effect
  rownames(Y) <- c(paste0("c", seq_len(n_ctrl)), "treated")
  colnames(Y) <- trend
  list(Y = Y, N0 = n_ctrl, T0 = T_pre)
}

test_that("untrimmed reproduces synthdid's own placebo SE exactly", {
  # The reimplementation must match the package before trimming means
  # anything -- otherwise a difference is a bug, not a method. The draws
  # consume the RNG exactly as synthdid's do, so the match is exact.
  d <- mixed_Y()
  est <- synthdid::synthdid_estimate(d$Y, d$N0, d$T0)
  set.seed(99); mine <- placebo_se_trimmed(est, replications = 400,
                                           trim_k = Inf)$se
  set.seed(99); theirs <- sqrt(stats::vcov(est, method = "placebo",
                                           replications = 400))
  expect_equal(mine, as.numeric(theirs), tolerance = 1e-10)
})

test_that("trimming shrinks the SE when the pool is heterogeneous", {
  d <- mixed_Y()
  est <- synthdid::synthdid_estimate(d$Y, d$N0, d$T0)
  set.seed(1); loose <- placebo_se_trimmed(est, replications = 300, trim_k = Inf)
  set.seed(1); tight <- placebo_se_trimmed(est, replications = 300, trim_k = 2)
  expect_lt(tight$se, loose$se)
  expect_lt(tight$n_kept, tight$n_total)
  expect_equal(loose$n_kept, loose$n_total)
})

test_that("trimming barely matters when every donor fits well", {
  # All controls on the same trend: nothing to trim, so the SE should move
  # very little. This is the guard against trimming being a free lunch.
  set.seed(3)
  n <- 20; trend <- 1:11
  Y <- outer(seq(50, by = 1.5, length.out = n + 1), trend, \(l, t) l + 0.5 * t)
  Y <- Y + matrix(stats::rnorm(length(Y), sd = 0.3), nrow = nrow(Y))
  Y[n + 1, 9:11] <- Y[n + 1, 9:11] + 3
  rownames(Y) <- c(paste0("c", seq_len(n)), "treated")
  est <- synthdid::synthdid_estimate(Y, n, 8)
  set.seed(5); loose <- placebo_se_trimmed(est, replications = 200, trim_k = Inf)$se
  set.seed(5); tight <- placebo_se_trimmed(est, replications = 200, trim_k = 5)$se
  expect_equal(tight, loose, tolerance = 0.25)
})

test_that("an impossibly tight trim warns and falls back rather than failing", {
  d <- mixed_Y()
  est <- synthdid::synthdid_estimate(d$Y, d$N0, d$T0)
  expect_warning(r <- placebo_se_trimmed(est, replications = 30, trim_k = 1e-9),
                 "survive")
  expect_true(is.finite(r$se))
})

test_that("fit_sdid carries the trimming record through", {
  d <- mixed_Y()
  r <- fit_sdid(d$Y, d$N0, d$T0, replications = 50, trim_k = 2)
  expect_equal(r$trim_k, 2)
  expect_lte(r$placebos_kept, r$placebos_total)
  expect_equal(r$placebos_total, 50)
  expect_equal(r$t, r$att / r$se, tolerance = 1e-8)
})
