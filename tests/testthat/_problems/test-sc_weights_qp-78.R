# Extracted from test-sc_weights_qp.R:78

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

# test -------------------------------------------------------------------------
withr::local_seed(31)
excess <- vapply(1:40, function(rep) {
    p <- random_problem(rep); lambda <- if (rep %% 2) 0 else 10^stats::runif(1, -4, 0)
    f_ref <- objective(p$Z0, p$Z1, p$v, reference_fw(p$Z0, p$Z1, p$v, lambda), lambda)
    f_got <- objective(p$Z0, p$Z1, p$v, sc_weights_simplex(p$Z0, p$Z1, p$v, penalty = lambda), lambda)
    (f_got - f_ref) / max(1, f_ref)
  }, numeric(1))
expect_lt(max(excess), 1e-6)
