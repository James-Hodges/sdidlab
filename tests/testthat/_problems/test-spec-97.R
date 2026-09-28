# Extracted from test-spec.R:97

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
specs <- list(
    sdid_spec(label = "A", design = list(treat_time = 1989), inference = list(replications = 100)),
    sdid_spec(label = "B", design = list(treat_time = 1989), inference = list(replications = 100), estimator = "sc"),
    sdid_spec(label = "C", design = list(treat_time = 1989), inference = list(replications = 100), estimator = "did")
  )
text <- write_specs_yaml(specs, NULL)
raw <- yaml::yaml.load(text)
expect_equal(raw$defaults$design$treat_time, 1989)
expect_null(raw$specs[[1]]$design)
