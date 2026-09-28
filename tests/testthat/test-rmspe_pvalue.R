ratio <- c(wm = 6.2, a = 1.1, b = 2.0, c = 0.7)
att   <- c(28.9, -8.0, 4.0, 1.0)

test_that("a two-sided p-value is the rank of the unsigned ratio", {
  out <- rmspe_pvalue(ratio, att, treated = "wm")
  expect_equal(out$rank, 1L)
  expect_equal(out$p_value, 0.25)   # 1 of 4, the floor
  expect_equal(out$floor, 0.25)
  expect_equal(out$n_units, 4L)
})

test_that("a one-sided test demotes placebos that moved the other way", {
  # b's ratio (2.0) beats the treated unit's on a two-sided test only when
  # its sign is ignored; here everything is positive, so nothing changes.
  r <- c(wm = 2.0, a = 3.0, b = 1.0, c = 0.5)
  a <- c(10, -20, 5, 1)
  two <- rmspe_pvalue(r, a, treated = "wm", alternative = "two.sided")
  one <- rmspe_pvalue(r, a, treated = "wm", alternative = "greater")
  expect_equal(two$rank, 2L)        # a's 3.0 outranks it
  expect_equal(one$rank, 1L)        # a moved down, so it no longer does
  expect_lt(one$p_value, two$p_value)
})

test_that("alternative = 'less' mirrors 'greater'", {
  r <- c(wm = 2.0, a = 3.0, b = 1.0)
  a <- c(-10, 20, -5)
  expect_equal(rmspe_pvalue(r, a, treated = "wm", alternative = "less")$rank, 1L)
  expect_equal(rmspe_pvalue(r, a, treated = "wm", alternative = "greater")$rank, 3L)
})

test_that("ties count against the treated unit", {
  out <- rmspe_pvalue(c(wm = 2, a = 2, b = 1), c(1, 1, 1), treated = "wm")
  expect_equal(out$rank, 2L)
})

test_that("the p-value can never reach zero", {
  out <- rmspe_pvalue(c(wm = 1e6, a = 1, b = 1), c(1, 1, 1), treated = "wm")
  expect_equal(out$p_value, 1 / 3)
  expect_gt(out$p_value, 0)
})

test_that("a non-finite ratio is an error, not a silent drop", {
  expect_error(
    rmspe_pvalue(c(wm = 2, a = Inf, b = 1), c(1, 1, 1), treated = "wm"),
    "non-finite RMSPE ratio"
  )
})

test_that("the treated unit must exist and the vectors must align", {
  expect_error(rmspe_pvalue(ratio, att, treated = "nope"), "does not identify")
  expect_error(rmspe_pvalue(ratio, att[1:2], treated = "wm"), "same length")
})

test_that("the paper's reported p-values are attainable ranks out of 40", {
  # Table 1 reports 0.050 and 0.075 with N = 40: ranks 2 and 3.
  r <- c(wm = 10, stats::setNames(c(20, 15, seq(9, 1, length.out = 37)),
                                  paste0("d", 1:39)))
  a <- rep(1, 40)
  out <- rmspe_pvalue(r, a, treated = "wm")
  expect_equal(out$rank, 3L)
  expect_equal(out$p_value, 0.075)
})
