# Extracted from test-plots.R:54

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
scm <- suppressWarnings(fit_spec(
    sdid_spec(estimator = "scm", pool = list(units = c("Nevada", "Utah", "Oregon", "Montana", "Idaho",
                                                       "Colorado", "New Mexico", "Wyoming"))),
    prop99_panel()))
