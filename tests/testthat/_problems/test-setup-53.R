# Extracted from test-setup.R:53

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
st <- build_setup(sdid_spec(pool = list(units = c("Nevada", "Utah", "Oregon"))), prop99_panel())
