# The Abadie-L'Hour (2021) penalised synthetic control:
#   min_w (Z1 - Z0'w)' V (Z1 - Z0'w) + lambda * sum_j w_j (Z1 - Z0j)' V (Z1 - Z0j)
# over the simplex. Optimality itself is checked against the KKT conditions
# in test-sc_weights_qp.R; these tests pin down what the penalty is for.

pen_objective <- function(Z0, Z1, v, w, lambda) {
  d <- as.numeric((sweep(Z0, 2, Z1, "-")^2) %*% v)
  sum(v * (Z1 - as.numeric(crossprod(Z0, w)))^2) + lambda * sum(w * d)
}

# A target that several donor mixtures hit exactly: it is the midpoint of
# d1 and d2, and also of d3 and d4, and a combination of all four.
exact_problem <- function() {
  Z0 <- rbind(d1 = c(0, 0), d2 = c(2, 2), d3 = c(0, 2), d4 = c(2, 0),
              d5 = c(1.1, 0.9), d6 = c(9, 9))
  list(Z0 = Z0, Z1 = c(1, 1))
}

test_that("without a penalty an exact fit is not unique", {
  # Two different donor mixtures reproduce the target exactly, so the
  # unpenalised objective cannot choose between them -- yet they would imply
  # different counterfactuals. This is the non-identification the penalty
  # exists to fix.
  p <- exact_problem(); v <- rep(1, 2)
  a <- c(0.5, 0.5, 0, 0, 0, 0)
  b <- c(0, 0, 0.5, 0.5, 0, 0)
  expect_equal(pen_objective(p$Z0, p$Z1, v, a, 0), 0)
  expect_equal(pen_objective(p$Z0, p$Z1, v, b, 0), 0)
  # With a penalty they are no longer tied, and neither is the optimum.
  expect_gt(pen_objective(p$Z0, p$Z1, v, a, 1e-3),
            pen_objective(p$Z0, p$Z1, v, sc_weights_simplex(p$Z0, p$Z1, penalty = 1e-3), 1e-3))
})

test_that("a small penalty picks one solution regardless of donor order", {
  p <- exact_problem()
  a <- sc_weights_simplex(p$Z0, p$Z1, penalty = 1e-3)
  for (ord in list(c(3, 4, 1, 2, 6, 5), c(6, 5, 4, 3, 2, 1), c(2, 6, 1, 5, 3, 4))) {
    b <- sc_weights_simplex(p$Z0[ord, ], p$Z1, penalty = 1e-3)[order(ord)]
    expect_equal(b, a, tolerance = 1e-6)
  }
  # And it is the local one: d5 sits next to the target.
  expect_gt(a[5], 0.5)
})

test_that("a large penalty collapses to the nearest neighbour", {
  p <- exact_problem()
  w <- sc_weights_simplex(p$Z0, p$Z1, penalty = 1e3)
  expect_equal(unname(which.max(w)), 5L)
  expect_equal(max(w), 1, tolerance = 1e-8)
})

test_that("zero penalty is the classic synthetic control", {
  withr::local_seed(4)
  Z0 <- matrix(stats::rnorm(60), 20, 3); Z1 <- stats::rnorm(3)
  expect_identical(sc_weights_simplex(Z0, Z1, penalty = 0), sc_weights_simplex(Z0, Z1))
})

test_that("fit_scm threads the penalty through and stops warning on exact fits", {
  Y <- toy()
  Y["treated", ] <- 0.5 * Y["d1", ] + 0.5 * Y["d3", ]
  Y["treated", 7:8] <- Y["treated", 7:8] + 5
  expect_warning(fit_scm(Y, N0 = 5, T0 = 6, bias_correct = FALSE), class = "sdidlab_exact_fit")
  expect_no_warning(fit <- fit_scm(Y, N0 = 5, T0 = 6, bias_correct = FALSE, penalty = 1e-4))
  expect_equal(sum(fit$omega), 1)
})

test_that("an invalid penalty is rejected", {
  expect_error(sc_weights_simplex(diag(2), c(0.5, 0.5), penalty = -1), "non-negative")
  expect_error(sc_weights_simplex(diag(2), c(0.5, 0.5), penalty = NA_real_), "non-negative")
})
