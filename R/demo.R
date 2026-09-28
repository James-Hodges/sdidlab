#' Launch the interactive demo
#'
#' A Shiny app for building a specification, fitting it, inspecting every
#' diagnostic plot, comparing specifications, and downloading them as YAML
#' plus an R script that reproduces the fits. It runs on the built-in
#' Proposition 99 panel or an uploaded CSV.
#'
#' @param ... Passed to [shiny::runApp()], e.g. `port` or `launch.browser`.
#'
#' @return Called for its side effect of running the app.
#'
#' @examples
#' if (interactive()) run_demo()
#'
#' @export
run_demo <- function(...) {
  for (pkg in c("shiny", "bslib")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort("The demo needs {.pkg {pkg}}: {.code install.packages(\"{pkg}\")}.")
    }
  }
  shiny::runApp(system.file("shiny", "demo", package = "sdidlab"), ...)
}
