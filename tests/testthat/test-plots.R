run_sdid <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- fit_spec(sdid_spec(label = "plots", inference = list(replications = 30)),
                                           prop99_panel(), seed = 1)
    cache
  }
})

builds <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(p)))
}

test_that("every SDID plot builds", {
  run <- run_sdid()
  builds(plot_sdid_trajectory(run))
  builds(plot_gap_path(run))
  builds(plot_event_study(run))
  builds(plot_donor_weights(run))
  builds(plot_lambda_weights(run))
  builds(plot_sdid_units(run))
  builds(plot_placebo_distribution(run))
  builds(plot_placebo_trim(run))
  builds(plot_donor_influence(run))
})

test_that("the gap path's post-period points average to the ATT", {
  run <- run_sdid()
  d <- plot_gap_path(run)$data
  post <- d$period >= 1989
  expect_equal(mean(d$gap[post]), run$fit$att, tolerance = 1e-10)
})

test_that("the event study is synthdid's effect curve", {
  run <- run_sdid()
  d <- plot_event_study(run)$data
  expect_equal(d$effect, as.numeric(synthdid::synthdid_effect_curve(run$fit$est)))
  expect_equal(mean(d$effect), run$fit$att, tolerance = 1e-10)
})

test_that("plots use display labels from the panel", {
  p <- prop99_panel()
  p$data$name <- paste("State of", p$data$State)
  p$label <- "name"
  run <- fit_spec(sdid_spec(inference = list(se_method = "none")), p)
  expect_true(all(grepl("^State of", donor_weight_table(run)$label)))
})

test_that("SCM runs have their own plots", {
  scm <- suppressWarnings(fit_spec(
    sdid_spec(estimator = "scm", pool = list(units = c("Nevada", "Utah", "Montana", "Idaho", "Kansas",
                                                       "Colorado", "New Mexico", "Wyoming"))),
    prop99_panel()))
  builds(plot_scm_gaps(scm))
  builds(plot_gap_path(scm))
  builds(plot_event_study(scm))
  builds(plot_donor_weights(scm))
  expect_error(plot_placebo_trim(scm), "does not trim")
  expect_error(plot_scm_gaps(run_sdid()), "scm")
})

test_that("the donor influence shares each sum to one", {
  d <- donor_influence(run_sdid())
  expect_equal(sum(d$weight_share), 1)
  expect_equal(sum(d$se_share), 1, tolerance = 1e-10)
})

test_that("the forest plot checks its columns", {
  builds(plot_att_forest(data.frame(label = c("a", "b"), att = c(1, 2), ci_lo = c(0, 1), ci_hi = c(2, 3))))
  expect_error(plot_att_forest(data.frame(label = "a", att = 1)), "missing")
})
