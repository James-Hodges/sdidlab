# The simplex solve must reach the optimum of the quadratic program it
# claims to solve. quadprog solves the same program exactly.

objective <- function(Z0, Z1, v, w) sum(v * (Z1 - as.numeric(crossprod(Z0, w)))^2)

qp_weights <- function(Z0, Z1, v) {
  n <- nrow(Z0)
  D <- 2 * (Z0 %*% (v * t(Z0))) + diag(1e-10, n)
  d <- 2 * as.numeric(Z0 %*% (v * Z1))
  w <- quadprog::solve.QP(D, d, cbind(1, diag(n)), c(1, rep(0, n)), meq = 1)$solution
  w <- pmax(w, 0)
  w / sum(w)
}

test_that("Frank-Wolfe reaches the exact QP optimum on random problems", {
  skip_if_not_installed("quadprog")
  withr::local_seed(2024)
  excess <- vapply(1:100, function(rep) {
    n <- sample(5:60, 1); k <- sample(2:15, 1)
    Z0 <- matrix(stats::rnorm(n * k), n, k)
    Z1 <- if (rep %% 2) stats::rnorm(k) else as.numeric(crossprod(Z0, prop.table(stats::runif(n))))
    v  <- if (rep %% 3) rep(1, k) else exp(stats::runif(k, -3, 3))
    objective(Z0, Z1, v, sc_weights_simplex(Z0, Z1, v)) - objective(Z0, Z1, v, qp_weights(Z0, Z1, v))
  }, numeric(1))
  expect_lt(max(excess), 1e-6)
})

prop99_predictors <- function() {
  data("california_prop99", package = "synthdid", envir = environment())
  s <- synthdid::panel.matrices(california_prop99)
  Z <- s$Y[, 1:s$T0]
  Z <- sweep(Z, 2, apply(Z, 2, stats::sd), "/")
  list(Z0 = Z[1:s$N0, ], Z1 = Z[s$N0 + 1, ], s = s)
}

test_that("the Proposition 99 weights converge to the QP optimum", {
  # This panel needs ~12,500 iterations; the old cap of 500 stopped 2.7%
  # above the optimum and moved the synthetic-control ATT by half a pack.
  skip_if_not_installed("quadprog")
  p <- prop99_predictors(); v <- rep(1, ncol(p$Z0))
  f_fw <- objective(p$Z0, p$Z1, v, sc_weights_simplex(p$Z0, p$Z1, v))
  f_qp <- objective(p$Z0, p$Z1, v, qp_weights(p$Z0, p$Z1, v))
  expect_equal(f_fw, f_qp, tolerance = 1e-7)
})

test_that("hitting the iteration cap warns when asked to", {
  p <- prop99_predictors()
  expect_warning(sc_weights_simplex(p$Z0, p$Z1, iter = 50L, warn = TRUE), "did not converge")
  expect_no_warning(sc_weights_simplex(p$Z0, p$Z1, iter = 50L))
})

test_that("equal-V SCM fits at least as well as Synth with the same V", {
  skip_if_not_installed("Synth")
  p <- prop99_predictors(); s <- p$s; pre <- 1:s$T0
  fit <- fit_scm(s$Y, s$N0, s$T0, bias_correct = FALSE, v_method = "equal")
  syn <- Synth::synth(X1 = matrix(p$Z1, ncol = 1), X0 = t(p$Z0),
                      Z1 = matrix(s$Y[s$N0 + 1, pre], ncol = 1), Z0 = t(s$Y[1:s$N0, pre]),
                      custom.v = rep(1, s$T0), verbose = FALSE)
  v <- rep(1, s$T0)
  expect_lte(objective(p$Z0, p$Z1, v, fit$omega),
             objective(p$Z0, p$Z1, v, as.numeric(syn$solution.w)) * (1 + 1e-8))
})
