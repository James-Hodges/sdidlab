# Extracted from test-sc_weights_qp.R:48

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
withr::local_seed(7)
worst <- max(vapply(1:100, function(rep) {
    p <- random_problem(rep); lambda <- 10^stats::runif(1, -4, 1)
    kkt_violation(p$Z0, p$Z1, p$v, sc_weights_simplex(p$Z0, p$Z1, p$v, penalty = lambda), lambda)
  }, numeric(1)))
expect_lt(worst, 1e-5)
