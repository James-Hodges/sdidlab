covariate_panel <- function(n_ctrl = 15, T_pre = 8, T_post = 3, K = 2, seed = 5) {
  withr::with_seed(seed, {
    n <- n_ctrl + 1; Tt <- T_pre + T_post
    X <- array(stats::rnorm(n * Tt * K), dim = c(n, Tt, K),
               dimnames = list(NULL, NULL, paste0("cov", seq_len(K))))
    beta <- c(2, -1.5)[seq_len(K)]
    Y <- outer(seq(10, by = 1, length.out = n), seq_len(Tt), "+") +
      matrix(stats::rnorm(n * Tt, sd = 0.3), n, Tt)
    for (k in seq_len(K)) Y <- Y + beta[k] * X[, , k]
    Y[n, (T_pre + 1):Tt] <- Y[n, (T_pre + 1):Tt] + 4   # the effect
    list(Y = Y, X = X, N0 = n_ctrl, T0 = T_pre, beta = beta)
  })
}

test_that("it runs on a 3-D covariate array and recovers the coefficients", {
  d <- covariate_panel()
  Yres <- project_out_covariates(d$Y, d$X, d$N0, d$T0)
  expect_equal(dim(Yres), dim(d$Y))
  expect_equal(unname(attr(Yres, "beta")), d$beta, tolerance = 0.1)
  expect_named(attr(Yres, "beta"), c("cov1", "cov2"))
})

test_that("coefficients come from untreated cells only", {
  # Two-way fixed effects OLS on the untreated cells, done independently.
  d <- covariate_panel()
  n <- nrow(d$Y); Tt <- ncol(d$Y)
  long <- data.frame(unit = factor(rep(seq_len(n), Tt)), time = factor(rep(seq_len(Tt), each = n)),
                     y = as.vector(d$Y), x1 = as.vector(d$X[, , 1]), x2 = as.vector(d$X[, , 2]))
  untreated <- as.integer(long$unit) <= d$N0 | as.integer(long$time) <= d$T0
  b <- stats::coef(stats::lm(y ~ x1 + x2 + unit + time, data = long[untreated, ]))[c("x1", "x2")]
  Yres <- project_out_covariates(d$Y, d$X, d$N0, d$T0)
  expect_equal(unname(attr(Yres, "beta")), unname(b), tolerance = 1e-10)
  expect_equal(Yres, d$Y - b[[1]] * d$X[, , 1] - b[[2]] * d$X[, , 2],
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("it agrees with fixest", {
  skip_if_not_installed("fixest")
  d <- covariate_panel(K = 1)
  long <- data.frame(unit = rep(seq_len(nrow(d$Y)), ncol(d$Y)),
                     time = rep(seq_len(ncol(d$Y)), each = nrow(d$Y)),
                     y = as.vector(d$Y), x = as.vector(d$X[, , 1]))
  long <- long[long$unit <= d$N0 | long$time <= d$T0, ]
  fe <- fixest::feols(y ~ x | unit + time, data = long)
  Yres <- project_out_covariates(d$Y, d$X, d$N0, d$T0)
  expect_equal(unname(attr(Yres, "beta")), unname(stats::coef(fe)[["x"]]), tolerance = 1e-8)
})

test_that("a single covariate may be given as a matrix", {
  d <- covariate_panel(K = 1)
  expect_equal(project_out_covariates(d$Y, d$X[, , 1], d$N0, d$T0),
               project_out_covariates(d$Y, d$X, d$N0, d$T0), ignore_attr = TRUE)
})

test_that("the treated unit's post-period does not leak into the coefficients", {
  d <- covariate_panel(K = 1)
  shocked <- d$Y
  shocked[nrow(shocked), (d$T0 + 1):ncol(shocked)] <- shocked[nrow(shocked), (d$T0 + 1):ncol(shocked)] + 1000
  expect_equal(attr(project_out_covariates(shocked, d$X, d$N0, d$T0), "beta"),
               attr(project_out_covariates(d$Y, d$X, d$N0, d$T0), "beta"))
})

test_that("mismatched dimensions are rejected", {
  d <- covariate_panel(K = 1)
  expect_error(project_out_covariates(d$Y, d$X[-1, , , drop = FALSE], d$N0, d$T0), "per covariate")
})
