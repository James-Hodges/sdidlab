# Extracted from test-setup.R:74

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
register_pool("test_west", function(panel, spec) {
    data.frame(unit = c("Oregon", "Nevada", "Utah"), included = c(TRUE, TRUE, FALSE),
               reason = c(NA, NA, "not western enough"))
  })
register_outcome("test_log_packs", function(panel, spec) {
    d <- panel$data
    data.frame(unit = d$State, time = d$Year, y = log(d$PacksPerCapita))
  }, label = "Log packs")
st <- build_setup(sdid_spec(pool = list(base = "test_west"), data = list(outcome = "test_log_packs")),
                    prop99_panel())
expect_equal(st$N0, 2)
