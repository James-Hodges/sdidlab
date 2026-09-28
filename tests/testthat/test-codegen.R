test_that("a generated script parses and reproduces the fits exactly", {
  specs <- list(
    sdid_spec(label = "SDID", inference = list(replications = 20)),
    sdid_spec(label = "SC", estimator = "sc", inference = list(replications = 20),
              pool = list(drop_units = "Nevada", exclusion_notes = list(Nevada = "note with \"quotes\"")))
  )
  code <- write_spec_script(specs, NULL, data_code = "panel <- prop99_panel()", title = "Test")
  expect_no_error(parse(text = code))

  # Run it without the library() call, against the package under test.
  env <- new.env(parent = asNamespace("sdidlab"))
  script <- sub("library(sdidlab)", "", code, fixed = TRUE)
  suppressMessages(capture.output(eval(parse(text = script), envir = env)))
  direct <- run_specs(specs, prop99_panel(), progress = FALSE)
  expect_equal(vapply(env$runs, function(r) r$fit$att, 0), vapply(direct, function(r) r$fit$att, 0))
  expect_equal(vapply(env$runs, function(r) r$fit$se, 0), vapply(direct, function(r) r$fit$se, 0))
  expect_identical(vapply(env$runs, `[[`, "", "spec_hash"), vapply(specs, spec_hash, ""))
})

test_that("a script can refer to an external spec file", {
  code <- write_spec_script(sdid_spec(), NULL, data_code = "panel <- prop99_panel()",
                            spec_file = "specs/prop99.yml")
  expect_match(code, 'read_specs("specs/prop99.yml")', fixed = TRUE)
  expect_false(grepl("spec_yaml", code))
})

test_that("bad data_code is rejected", {
  expect_error(write_spec_script(sdid_spec(), NULL, data_code = c("a", "b")), "single string")
})
