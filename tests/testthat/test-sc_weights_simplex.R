# The weight solve must be optimal for the program it claims to solve. That
# is checked here against the KKT conditions of the simplex-constrained
# problem, which do not depend on how the solution was found:
#   for every donor with positive weight the gradient equals a common value
#   mu, and no donor has a gradient below mu.

objective <- function(Z0, Z1, v, w, penalty = 0) {
  d <- as.numeric((sweep(Z0, 2, Z1, "-")^2) %*% v)
  sum(v * (Z1 - as.numeric(crossprod(Z0, w)))^2) + penalty * sum(w * d)
}

kkt_violation <- function(Z0, Z1, v, w, penalty = 0) {
  d <- as.numeric((sweep(Z0, 2, Z1, "-")^2) %*% v)
  resid <- Z1 - as.numeric(crossprod(Z0, w))
  grad <- -2 * as.numeric((t(t(Z0) * v)) %*% resid) + penalty * d
  active <- w > 1e-6
  mu <- mean(grad[active])
  scale <- max(1, abs(grad))
  max(abs(grad[active] - mu), pmax(mu - grad, 0)) / scale
}

random_problem <- function(rep) {
  n <- sample(5:80, 1); k <- sample(2:20, 1)
  Z0 <- matrix(stats::rnorm(n * k), n, k)
  Z1 <- if (rep %% 2) stats::rnorm(k) else as.numeric(crossprod(Z0, prop.table(stats::runif(n))))
  v  <- if (rep %% 3) rep(1, k) else exp(stats::runif(k, -3, 3))
  list(Z0 = Z0, Z1 = Z1, v = v)
}

test_that("weights satisfy the KKT conditions on random problems", {
  withr::local_seed(2024)
  worst <- max(vapply(1:150, function(rep) {
    p <- random_problem(rep)
    w <- sc_weights_simplex(p$Z0, p$Z1, p$v)
    expect_equal(sum(w), 1)
    expect_true(all(w >= 0))
    kkt_violation(p$Z0, p$Z1, p$v, w)
  }, numeric(1)))
  expect_lt(worst, 1e-5)
})

test_that("penalised weights satisfy the penalised KKT conditions", {
  withr::local_seed(7)
  worst <- max(vapply(1:100, function(rep) {
    p <- random_problem(rep); lambda <- 10^stats::runif(1, -4, 1)
    kkt_violation(p$Z0, p$Z1, p$v, sc_weights_simplex(p$Z0, p$Z1, p$v, penalty = lambda), lambda)
  }, numeric(1)))
  expect_lt(worst, 1e-5)
})

# An independent reference: a plain pairwise Frank-Wolfe solve of the same
# (penalised) program, run for as long as it takes. Slow, but it shares no
# code with the solver under test.
reference_fw <- function(Z0, Z1, v, penalty = 0, iter = 2e5) {
  d <- as.numeric((sweep(Z0, 2, Z1, "-")^2) %*% v)
  w <- numeric(nrow(Z0)); w[which.min(d)] <- 1
  for (t in seq_len(iter)) {
    r <- Z1 - as.numeric(crossprod(Z0, w))
    g <- -2 * as.numeric(Z0 %*% (v * r)) + penalty * d
    s <- which.min(g); act <- which(w > 0); a <- act[which.max(g[act])]
    if (s == a || g[a] - g[s] < 1e-12) break
    u <- Z0[s, ] - Z0[a, ]
    step <- min(w[a], max(0, (sum(v * u * r) - penalty * (d[s] - d[a]) / 2) / sum(v * u * u)))
    if (step <= 0) break
    w[s] <- w[s] + step; w[a] <- w[a] - step
  }
  w
}

test_that("the objective matches an independent long-run reference solve", {
  withr::local_seed(31)
  excess <- vapply(1:40, function(rep) {
    p <- random_problem(rep); lambda <- if (rep %% 2) 0 else 10^stats::runif(1, -4, 0)
    f_ref <- objective(p$Z0, p$Z1, p$v, reference_fw(p$Z0, p$Z1, p$v, lambda), lambda)
    f_got <- objective(p$Z0, p$Z1, p$v, sc_weights_simplex(p$Z0, p$Z1, p$v, penalty = lambda), lambda)
    (f_got - f_ref) / max(1, f_ref)
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

test_that("the Proposition 99 weights are optimal", {
  # A first-order solve needed >12,000 iterations here; a 500-iteration cap
  # stopped 2.7% above this optimum (0.05638 vs 0.05790).
  p <- prop99_predictors(); v <- rep(1, ncol(p$Z0))
  w <- sc_weights_simplex(p$Z0, p$Z1, v)
  expect_lt(kkt_violation(p$Z0, p$Z1, v, w), 1e-6)
  expect_equal(objective(p$Z0, p$Z1, v, w), 0.0563814, tolerance = 1e-6)
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
