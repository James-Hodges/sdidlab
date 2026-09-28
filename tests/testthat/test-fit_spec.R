fast <- function(...) {
  args <- list(...)
  args$inference <- utils::modifyList(list(replications = 20), args$inference %||% list())
  do.call(sdid_spec, args)
}

test_that("fit_spec reproduces synthdid on Prop 99", {
  data("california_prop99", package = "synthdid", envir = environment())
  s <- synthdid::panel.matrices(california_prop99)
  ref <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0)
  set.seed(8); theirs <- sqrt(stats::vcov(ref, method = "placebo", replications = 20))
  run <- fit_spec(fast(inference = list(replications = 20, trim_k = Inf)), prop99_panel(), seed = 8)
  expect_equal(run$fit$att, as.numeric(ref), tolerance = 1e-12)
  expect_equal(run$fit$se, as.numeric(theirs), tolerance = 1e-10)
})

test_that("a seeded fit is reproducible and leaves the caller's RNG alone", {
  set.seed(1); before <- runif(1); set.seed(1)
  a <- fit_spec(fast(), prop99_panel(), seed = 5)
  expect_equal(runif(1), before)
  b <- fit_spec(fast(), prop99_panel(), seed = 5)
  expect_equal(a$fit$se, b$fit$se)
})

test_that("seed_for depends only on the spec's identity", {
  expect_identical(seed_for(fast(label = "a")), seed_for(fast(label = "b")))
  expect_false(identical(seed_for(fast()), seed_for(fast(estimator = "sc"))))
})

test_that("every estimator runs through fit_spec", {
  p <- prop99_panel()
  for (est in c("sdid", "sc", "did")) {
    run <- fit_spec(fast(estimator = est), p, seed = 1)
    expect_true(is.finite(run$fit$att) && is.finite(run$fit$se), label = est)
  }
  scm <- suppressWarnings(fit_spec(sdid_spec(estimator = "scm", inference = list(se_method = "none")), p))
  expect_true(is.finite(scm$fit$att))
  expect_equal(scm$fit$estimator, "scm")
})

test_that("both covariate routes run and the projected one uses no X", {
  p <- prop99_panel()
  p$data$x <- withr::with_seed(4, stats::rnorm(nrow(p$data))) + 0.02 * p$data$PacksPerCapita
  opt <- fit_spec(fast(covariates = list(vars = "x", mode = "optimized")), p, seed = 1)
  prj <- fit_spec(fast(covariates = list(vars = "x", mode = "projected")), p, seed = 1)
  expect_true(is.finite(opt$fit$att) && is.finite(prj$fit$att))
  expect_false(isTRUE(all.equal(opt$fit$att, prj$fit$att)))
  expect_equal(dim(attr(opt$fit$est, "setup")$X)[3], 1)
})

test_that("run_specs skips a failing specification and carries on", {
  specs <- list(fast(label = "ok"), fast(label = "bad", pool = list(units = c("Utah", "Nevada"), min_donors = 3)))
  runs <- suppressMessages(run_specs(specs, prop99_panel(), progress = FALSE))
  expect_s3_class(runs[[1]], "sdid_run")
  expect_null(runs[[2]])
  tab <- sdid_results_table(runs)
  expect_equal(nrow(tab), 1)
  expect_equal(tab$label, "ok")
})

test_that("keep = 'none' drops the fitted object but keeps the summary", {
  run <- fit_spec(fast(), prop99_panel(), keep = "none", seed = 1)
  expect_null(run$fit$est)
  expect_true(is.finite(run$fit$att))
  expect_error(plot_gap_path(run), "keep")
})

test_that("the run store saves, reloads, indexes and prunes", {
  store <- withr::local_tempdir()
  p <- prop99_panel()
  runs <- lapply(1:3, function(i) fit_spec(fast(), p, store = store, seed = i))
  other <- fit_spec(fast(estimator = "sc"), p, store = store, seed = 1)
  idx <- refresh_runs_index(store)
  expect_equal(nrow(idx), 4)
  expect_true(file.exists(file.path(store, "index.csv")))
  back <- load_run(runs[[1]]$path)
  expect_equal(back$fit$att, runs[[1]]$fit$att)
  expect_equal(back$fit$placebo_draws, runs[[1]]$fit$placebo_draws)
  expect_identical(spec_hash(back$spec), runs[[1]]$spec_hash)
  expect_s3_class(plot_gap_path(back), "ggplot")
  pruned <- suppressMessages(prune_runs(store, keep_latest_per_spec = 1))
  expect_equal(nrow(pruned), 2)
  expect_equal(nrow(refresh_runs_index(store)), 2)
})

test_that("an incomplete run directory is not loadable", {
  expect_error(load_run(withr::local_tempdir()), "No completed run")
})
