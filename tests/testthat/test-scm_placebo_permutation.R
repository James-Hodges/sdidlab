# Noise is deliberate -- see the note on the fixture in test-fit_scm.R.
perm_toy <- function() {
  periods <- 2010:2017
  Y <- outer(c(10, 12, 14, 16, 18, 13), rep(1, length(periods))) +
    outer(rep(1, 6), seq_along(periods)) +
    withr::with_seed(42, matrix(stats::rnorm(6 * length(periods), sd = 0.4),
                                nrow = 6))
  Y[6, 7:8] <- Y[6, 7:8] + 5
  dimnames(Y) <- list(c(paste0("d", 1:5), "treated"), periods)
  Y
}

test_that("every unit gets a turn as the treated one", {
  out <- scm_placebo_permutation(perm_toy(), N0 = 5, T0 = 6,
                                 bias_correct = FALSE)
  expect_equal(nrow(out$results), 6L)
  expect_equal(sum(out$results$treated), 1L)
  expect_equal(out$results$unit[out$results$treated], "treated")
  expect_equal(dim(out$gaps), c(6L, 8L))
  expect_equal(rownames(out$gaps), rownames(perm_toy()))
})

test_that("the treated unit is kept out of every placebo's donor pool", {
  # If it were left in, a placebo could load onto the treated unit's
  # post-period jump and inflate the null distribution.
  Y <- perm_toy()
  Y[6, 7:8] <- Y[6, 7:8] + 500          # make contamination impossible to miss
  out <- scm_placebo_permutation(Y, N0 = 5, T0 = 6, bias_correct = FALSE)
  placebo_gaps <- out$gaps[!out$results$treated, 7:8]
  expect_true(all(abs(placebo_gaps) < 50))
})

test_that("the real treated unit stands out on the RMSPE ratio", {
  out <- scm_placebo_permutation(perm_toy(), N0 = 5, T0 = 6,
                                 bias_correct = FALSE)
  p <- rmspe_pvalue(stats::setNames(out$results$rmspe_ratio, out$results$unit),
                    out$results$att, treated = "treated")
  expect_equal(p$rank, 1L)
})

test_that("results are deterministic across calls", {
  # Some placebo fits on a five-donor toy land on a single donor, which
  # legitimately warns that the bias-correction regression was skipped.
  # The point here is reproducibility, so the warning is not the subject.
  a <- suppressWarnings(scm_placebo_permutation(perm_toy(), N0 = 5, T0 = 6))
  b <- suppressWarnings(scm_placebo_permutation(perm_toy(), N0 = 5, T0 = 6))
  expect_equal(a$results, b$results)
  expect_equal(a$gaps, b$gaps)
})

test_that("a predictor matrix is permuted in step with the outcome rows", {
  Y <- perm_toy()
  X <- cbind(size = c(100, 200, 300, 400, 500, 250))
  out <- scm_placebo_permutation(Y, N0 = 5, T0 = 6, X = X,
                                 bias_correct = FALSE)
  expect_equal(nrow(out$results), 6L)
  expect_true(all(is.finite(out$results$att)))
})
