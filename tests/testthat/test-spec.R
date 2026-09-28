test_that("sdid_panel validates its input", {
  d <- data.frame(u = rep(c("a", "b"), each = 3), t = rep(1:3, 2), y = 1:6)
  expect_s3_class(sdid_panel(d, "u", "t", "y"), "sdid_panel")
  expect_error(sdid_panel(d, "u", "t", "nope"), "not found")
  expect_error(sdid_panel(rbind(d, d[1, ]), "u", "t", "y"), "more than one row")
  d2 <- d; d2$t <- as.character(d2$t)
  expect_error(sdid_panel(d2, "u", "t", "y"), "numeric")
  d3 <- d; d3$tr <- 2
  expect_error(sdid_panel(d3, "u", "t", "y", treatment = "tr"), "0/1")
})

test_that("unit identifiers are stored as character", {
  d <- data.frame(u = rep(c(10L, 20L), each = 3), t = rep(1:3, 2), y = 1:6)
  expect_type(sdid_panel(d, "u", "t", "y")$data$u, "character")
})

test_that("prop99_panel wraps synthdid's data", {
  p <- prop99_panel()
  expect_equal(length(unique(p$data$State)), 39)
  expect_equal(p$treatment, "treated")
})

test_that("the default specification is valid and prints", {
  expect_s3_class(sdid_spec(), "sdid_spec")
  expect_output(print(sdid_spec(label = "x")), NA)
  expect_no_error(capture.output(print(sdid_spec(label = "x"))))
})

test_that("overrides merge one level deep, two for the window", {
  s <- sdid_spec(design = list(window = list(start = 1975)))
  expect_equal(s$design$window$start, 1975)
  expect_true(is.na(s$design$window$end))
  expect_equal(s$design$anticipation, 0)
})

test_that("unknown and duplicated fields are errors, not silently dropped", {
  expect_error(sdid_spec(desgin = list()), "Unknown spec field")
  expect_error(sdid_spec(design = list(treat_yr = 1)), "Unknown field")
  expect_error(do.call(sdid_spec, list(design = list(), design = list())), "more than once")
})

test_that("the hash ignores metadata, type and construction order", {
  a <- sdid_spec(design = list(treat_time = 1989L, window = list(start = 1970)), label = "a")
  b <- sdid_spec(label = "b", note = "different", design = list(window = list(start = 1970L), treat_time = 1989))
  expect_identical(spec_hash(a), spec_hash(b))
  expect_false(identical(spec_hash(a), spec_hash(sdid_spec(design = list(treat_time = 1990)))))
})

test_that("unit lists hash the same in any order", {
  a <- sdid_spec(pool = list(drop_units = c("x", "y")))
  b <- sdid_spec(pool = list(drop_units = c("y", "x", "x")))
  expect_identical(spec_hash(a), spec_hash(b))
})

test_that("the default specification hashes to its defaults alone", {
  # Adding a field with a default must not move existing hashes; the empty
  # payload is what guarantees that.
  expect_identical(spec_hash(sdid_spec()), spec_hash(sdid_spec(label = "anything")))
})

test_that("validation rejects impossible designs", {
  expect_error(sdid_spec(estimator = "ols"), "estimator")
  expect_error(sdid_spec(covariates = list(mode = "projected")), "no")
  expect_error(sdid_spec(design = list(window = list(start = 2000, end = 1990))), "before its end")
  expect_error(sdid_spec(design = list(treat_time = 1970, window = list(start = 1970))), "inside the window")
  expect_error(sdid_spec(design = list(anticipation = -1)), ">= 0")
  expect_error(sdid_spec(design = list(treat_time = 1989, effect_time = 1989)), "strictly later")
  expect_error(sdid_spec(pool = list(drop_units = "a", add_units = "a")), "both")
  expect_error(sdid_spec(pool = list(base = "west", units = "a")), "not both")
  expect_error(sdid_spec(design = list(treated_unit = "a"), pool = list(add_units = "a")), "own donor pool")
  expect_error(sdid_spec(inference = list(conf_level = 1)), "between 0 and 1")
  expect_error(sdid_spec(estimator = "scm", covariates = list(vars = "x", mode = "projected")), "not supported")
})

test_that("a single specification round-trips through YAML", {
  s <- sdid_spec(label = "rt", design = list(treat_time = 1989, periods_drop = c(1980, 1981)),
                 pool = list(drop_units = c("Nevada", "Utah"), exclusion_notes = list(Nevada = "tax change")),
                 covariates = list(vars = c("a", "log_b"), log_vars = "b", mode = "projected"),
                 inference = list(trim_k = Inf), options = list(region = "west"))
  f <- withr::local_tempfile(fileext = ".yml")
  write_spec_yaml(s, f)
  back <- read_spec_yaml(f)
  expect_identical(spec_hash(back), spec_hash(s))
  expect_equal(back$pool$exclusion_notes$Nevada, "tax change")
  expect_equal(back$options$region, "west")
})

test_that("a sweep hoists shared fields and round-trips", {
  specs <- list(
    sdid_spec(label = "A", design = list(treat_time = 1989), inference = list(replications = 100)),
    sdid_spec(label = "B", design = list(treat_time = 1989), inference = list(replications = 100), estimator = "sc"),
    sdid_spec(label = "C", design = list(treat_time = 1989), inference = list(replications = 100), estimator = "did")
  )
  text <- write_specs_yaml(specs, NULL)
  raw <- yaml::yaml.load(text)
  expect_equal(raw$defaults$design$treat_time, 1989)
  expect_null(raw$specs[[1]]$design)
  f <- withr::local_tempfile(fileext = ".yml")
  write_specs_yaml(specs, f)
  expect_identical(vapply(read_specs(f), spec_hash, ""), vapply(specs, spec_hash, ""))
})

test_that("duplicate labels in a sweep are rejected", {
  f <- withr::local_tempfile(fileext = ".yml")
  writeLines("specs:\n- label: same\n- label: same\n  estimator: sc\n", f)
  expect_error(read_specs(f), "Duplicate label")
})

test_that("an invalid spec in a sweep is reported with its position", {
  f <- withr::local_tempfile(fileext = ".yml")
  writeLines("specs:\n- label: ok\n- label: bad\n  estimator: ols\n", f)
  expect_error(read_specs(f), "Spec 2")
})
