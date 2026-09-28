# Extracted from test-fit_spec.R:8

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "sdidlab", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
fast <- function(...) sdid_spec(inference = list(replications = 20), ...)

# test -------------------------------------------------------------------------
data("california_prop99", package = "synthdid", envir = environment())
s <- synthdid::panel.matrices(california_prop99)
ref <- synthdid::synthdid_estimate(s$Y, s$N0, s$T0)
set.seed(8)
theirs <- sqrt(stats::vcov(ref, method = "placebo", replications = 20))
run <- fit_spec(fast(inference = list(replications = 20, trim_k = Inf)), prop99_panel(), seed = 8)
