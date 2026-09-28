test_that("the Prop 99 setup is exactly synthdid's panel.matrices", {
  data("california_prop99", package = "synthdid", envir = environment())
  ref <- synthdid::panel.matrices(california_prop99)
  st <- build_setup(sdid_spec(), prop99_panel())
  expect_equal(st$Y, ref$Y)
  expect_equal(st$N0, ref$N0)
  expect_equal(st$T0, ref$T0)
  expect_equal(st$treated_unit, "California")
  expect_equal(st$treat_time, 1989)
})

test_that("explicit timing overrides the treatment column", {
  st <- build_setup(sdid_spec(design = list(treated_unit = "California", treat_time = 1985)), prop99_panel())
  expect_equal(st$T0, 15)
})

test_that("the window and dropped periods select the columns", {
  st <- build_setup(sdid_spec(design = list(window = list(start = 1975, end = 1995),
                                            periods_drop = c(1980, 1981))), prop99_panel())
  expect_equal(st$periods, setdiff(1975:1995, c(1980, 1981)))
  expect_equal(st$T0, sum(setdiff(1975:1988, c(1980, 1981)) < 1989))
})

test_that("an effect_time blacks out the phase-in periods", {
  st <- build_setup(sdid_spec(design = list(treat_time = 1989, effect_time = 1992)), prop99_panel())
  expect_equal(st$blackout, 1989:1991)
  expect_false(any(1989:1991 %in% st$periods))
  expect_equal(st$periods[st$T0 + 1], 1992)
})

test_that("anticipation moves the first treated period earlier", {
  st <- build_setup(sdid_spec(design = list(treat_time = 1989, anticipation = 2)), prop99_panel())
  expect_equal(st$periods[st$T0 + 1], 1987)
})

test_that("a panel without a treatment column needs explicit timing", {
  p <- prop99_panel(); p$treatment <- NULL
  expect_error(build_setup(sdid_spec(), p), "treated_unit")
  expect_no_error(build_setup(sdid_spec(design = list(treated_unit = "California", treat_time = 1989)), p))
})

test_that("pool edits are recorded, with reasons", {
  spec <- sdid_spec(pool = list(drop_units = c("Nevada", "Utah"),
                                exclusion_notes = list(Nevada = "tax change")))
  st <- build_setup(spec, prop99_panel())
  expect_equal(st$N0, 36)
  dropped <- st$pool[!st$pool$included, ]
  expect_setequal(dropped$unit, c("Nevada", "Utah"))
  expect_match(dropped$reason[dropped$unit == "Nevada"], "tax change")
})

test_that("an explicit unit list is the pool", {
  st <- build_setup(sdid_spec(pool = list(units = c("Nevada", "Utah", "Montana"))), prop99_panel())
  expect_equal(st$N0, 3)
  expect_setequal(rownames(st$Y)[1:3], c("Nevada", "Utah", "Montana"))
})

test_that("too few donors is an error that says why", {
  expect_error(build_setup(sdid_spec(pool = list(units = c("Nevada", "Utah"), min_donors = 3)),
                           prop99_panel()), "min_donors")
})

test_that("registered pools and outcomes are used", {
  register_pool("test_west", function(panel, spec) {
    data.frame(unit = c("Montana", "Nevada", "Utah"), included = c(TRUE, TRUE, FALSE),
               reason = c(NA, NA, "not western enough"))
  })
  register_outcome("test_log_packs", function(panel, spec) {
    d <- panel$data
    data.frame(unit = d$State, time = d$Year, y = log(d$PacksPerCapita))
  }, label = "Log packs")
  st <- build_setup(sdid_spec(pool = list(base = "test_west"), data = list(outcome = "test_log_packs")),
                    prop99_panel())
  expect_equal(st$N0, 2)
  expect_equal(st$pool$reason[st$pool$unit == "Utah"], "not western enough")
  expect_equal(st$Y["California", "1970"], log(123))
  expect_error(build_setup(sdid_spec(pool = list(base = "nope")), prop99_panel()), "not registered")
  expect_error(build_setup(sdid_spec(data = list(outcome = "nope")), prop99_panel()), "neither")
})

test_that("units missing a kept period are dropped for balance and recorded", {
  p <- prop99_panel()
  p$data$PacksPerCapita[p$data$State == "Utah" & p$data$Year == 1975] <- NA
  st <- build_setup(sdid_spec(), p)
  expect_false("Utah" %in% rownames(st$Y))
  expect_equal(st$balance$dropped$unit, "Utah")
  # ...but not when that period is dropped from the design.
  st2 <- build_setup(sdid_spec(design = list(periods_drop = 1975)), p)
  expect_true("Utah" %in% rownames(st2$Y))
})

test_that("a treated unit with a hole is an error, not a silent drop", {
  p <- prop99_panel()
  p$data$PacksPerCapita[p$data$State == "California" & p$data$Year == 1975] <- NA
  expect_error(build_setup(sdid_spec(), p), "treated unit")
})

test_that("covariates are aligned with Y and z-scored on the estimation sample", {
  p <- prop99_panel()
  p$data$x <- p$data$PacksPerCapita * 2 + 7
  st <- build_setup(sdid_spec(covariates = list(vars = "x", mode = "optimized")), p)
  expect_equal(dim(st$X), c(dim(st$Y), 1))
  z <- (st$Y * 2 + 7 - mean(st$Y * 2 + 7)) / stats::sd(as.vector(st$Y * 2 + 7))
  expect_equal(unname(st$X[, , 1]), unname(z), tolerance = 1e-10)
})

test_that("log_vars adds log1p columns", {
  p <- prop99_panel()
  p$data$x <- p$data$PacksPerCapita
  st <- build_setup(sdid_spec(covariates = list(vars = "log_x", log_vars = "x", mode = "optimized",
                                                z_score = FALSE)), p)
  expect_equal(unname(st$X[, , 1]), unname(log1p(st$Y)), tolerance = 1e-10)
})

test_that("a missing covariate cell is an error", {
  p <- prop99_panel()
  p$data$x <- p$data$PacksPerCapita
  p$data$x[5] <- NA
  expect_error(build_setup(sdid_spec(covariates = list(vars = "x", mode = "optimized")), p), "missing")
})

test_that("a pool naming units the outcome does not cover reports them", {
  register_pool("test_ghosts", function(panel, spec) c("Nevada", "Utah", "Atlantis"))
  st <- build_setup(sdid_spec(pool = list(base = "test_ghosts")), prop99_panel())
  expect_equal(st$N0, 2)
  expect_equal(st$balance$dropped$unit, "Atlantis")
})
