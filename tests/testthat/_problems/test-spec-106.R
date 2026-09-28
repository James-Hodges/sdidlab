# Extracted from test-spec.R:106

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
f <- withr::local_tempfile(fileext = ".yml")
writeLines("specs:\n- label: same\n- label: same\n  estimator: sc\n", f)
expect_error(read_specs(f), "Duplicate label")
