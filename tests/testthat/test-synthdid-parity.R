# sdidlab must agree with synthdid wherever it claims to reproduce it. These
# tests pin that down on the package's own reference panel, California's
# Proposition 99 (Abadie, Diamond and Hainmueller 2010), which synthdid ships.

prop99 <- function() {
  data("california_prop99", package = "synthdid", envir = environment())
  synthdid::panel.matrices(california_prop99)
}

# A time-varying covariate correlated with the outcome, fixed by seed.
prop99_X <- function(Y) {
  X <- withr::with_seed(99, matrix(stats::rnorm(length(Y)), nrow(Y), ncol(Y))) + 0.02 * Y
  array(X, dim = c(dim(Y), 1))
}

test_that("point estimates and weights match synthdid for sdid, sc and did", {
  s <- prop99()
  refs <- list(sdid = synthdid::synthdid_estimate, sc = synthdid::sc_estimate,
               did = synthdid::did_estimate)
  for (name in names(refs)) {
    ref <- refs[[name]](s$Y, s$N0, s$T0)
    got <- fit_sdid(s$Y, s$N0, s$T0, se_method = "none", estimator = name)
    w <- attr(ref, "weights")
    expect_equal(got$att, as.numeric(ref), tolerance = 1e-12, label = paste(name, "att"))
    expect_equal(unname(got$omega), as.numeric(w$omega), tolerance = 1e-12)
    expect_equal(got$lambda, as.numeric(w$lambda), tolerance = 1e-12)
  }
  # The published SDID estimate for Proposition 99.
  expect_equal(fit_sdid(s$Y, s$N0, s$T0, se_method = "none")$att, -15.6, tolerance = 0.01)
})

test_that("untrimmed placebo SE is synthdid's placebo SE, draw for draw", {
  s <- prop99()
  for (name in c("sdid", "sc", "did")) {
    ref <- switch(name, sdid = synthdid::synthdid_estimate, sc = synthdid::sc_estimate,
                  did = synthdid::did_estimate)(s$Y, s$N0, s$T0)
    set.seed(42); theirs <- sqrt(stats::vcov(ref, method = "placebo", replications = 40))
    set.seed(42); mine <- fit_sdid(s$Y, s$N0, s$T0, replications = 40, trim_k = Inf,
                                   estimator = name)
    expect_equal(mine$se, as.numeric(theirs), tolerance = 1e-10, label = paste(name, "se"))
    expect_equal(mine$placebos_kept, 40L)
  }
})

test_that("covariates passed as X match synthdid, estimate and placebo SE", {
  s <- prop99(); X <- prop99_X(s$Y)
  ref <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0, X = X)
  set.seed(3); theirs <- sqrt(stats::vcov(ref, method = "placebo", replications = 20))
  set.seed(3); mine <- fit_sdid(s$Y, s$N0, s$T0, X = X, replications = 20, trim_k = Inf)
  expect_equal(mine$att, as.numeric(ref), tolerance = 1e-12)
  expect_equal(mine$se, as.numeric(theirs), tolerance = 1e-10)
})

test_that("fit diagnostics are computed on the covariate-adjusted outcome", {
  s <- prop99(); X <- prop99_X(s$Y)
  r <- fit_sdid(s$Y, s$N0, s$T0, X = X, se_method = "none", keep = "est")
  beta <- attr(r$est, "weights")$beta
  Yadj <- s$Y - beta * X[, , 1]
  expect_equal(adjusted_outcome(r$est), Yadj)
  gap <- as.numeric(Yadj[s$N0 + 1, 1:s$T0]) -
    as.numeric(t(Yadj[1:s$N0, 1:s$T0]) %*% attr(r$est, "weights")$omega)
  expect_equal(r$trend_rmse, sqrt(mean((gap - mean(gap))^2)), tolerance = 1e-10)
  expect_equal(r$level_gap, mean(gap), tolerance = 1e-10)
})

test_that("adjusted_outcome is Y itself without covariates", {
  s <- prop99()
  est <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0)
  expect_identical(adjusted_outcome(est), s$Y)
})

test_that("SDID's per-period effect is the post gap minus the lambda-weighted pre gap", {
  # The identity the gap-path figure relies on to average to the ATT.
  s <- prop99()
  est <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0)
  w <- attr(est, "weights")
  gap <- as.numeric(s$Y[s$N0 + 1, ]) - as.numeric(t(s$Y[1:s$N0, ]) %*% w$omega)
  post <- (s$T0 + 1):ncol(s$Y)
  manual <- gap[post] - sum(w$lambda * gap[1:s$T0])
  expect_equal(as.numeric(synthdid::synthdid_effect_curve(est)), manual, tolerance = 1e-10)
  expect_equal(mean(manual), as.numeric(est), tolerance = 1e-10)
})
