#' Declare a long panel for sdidlab
#'
#' Wraps a long data frame -- one row per unit and period -- with the names
#' of its unit, time and outcome columns, so that specifications can refer to
#' columns by name. Any other columns (further outcomes, covariates, a
#' display label) stay available: a specification's `data$outcome` and
#' `covariates$vars` name them.
#'
#' Units are stored as character, so numeric identifiers and names behave
#' the same way everywhere, and times must be numeric (years, or any ordered
#' numeric period index).
#'
#' @param data A data frame in long format.
#' @param unit,time Names of the unit and time columns.
#' @param outcome Name of the default outcome column, used when a
#'   specification does not name one. Optional.
#' @param treatment Name of a 0/1 treatment column. Optional. When given, a
#'   specification that leaves `design$treated_unit` or `design$treat_time`
#'   unset takes them from it: the one unit that is ever treated, and its
#'   first treated period.
#' @param label Name of a column holding a display name for each unit, used
#'   in plots and tables. Optional.
#'
#' @return An `sdid_panel` object: a list with `data` (the data frame, with
#'   the unit column coerced to character) and the column names.
#'
#' @examples
#' data("california_prop99", package = "synthdid")
#' panel <- sdid_panel(california_prop99, unit = "State", time = "Year",
#'                     outcome = "PacksPerCapita", treatment = "treated")
#' panel
#'
#' @export
sdid_panel <- function(data, unit, time, outcome = NULL, treatment = NULL,
                       label = NULL) {
  if (!is.data.frame(data)) cli::cli_abort("{.arg data} must be a data frame.")
  cols <- c(unit = unit, time = time, outcome = outcome, treatment = treatment, label = label)
  missing <- cols[!cols %in% names(data)]
  if (length(missing)) {
    cli::cli_abort("Column{?s} {.field {missing}} not found in {.arg data}.")
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data[[unit]] <- as.character(data[[unit]])
  if (!is.numeric(data[[time]])) {
    cli::cli_abort(c(
      "The time column {.field {time}} must be numeric.",
      i = "Convert dates to a numeric period (a year, or a period index) first."
    ))
  }
  key <- paste(data[[unit]], data[[time]], sep = "\r")
  if (anyDuplicated(key)) {
    cli::cli_abort("{.arg data} has more than one row for some unit and period.")
  }
  if (!is.null(outcome) && !is.numeric(data[[outcome]])) {
    cli::cli_abort("The outcome column {.field {outcome}} must be numeric.")
  }
  if (!is.null(treatment) && !all(data[[treatment]] %in% c(0, 1, NA))) {
    cli::cli_abort("The treatment column {.field {treatment}} must be 0/1.")
  }
  structure(
    list(data = data, unit = unit, time = time, outcome = outcome,
         treatment = treatment, label = label),
    class = "sdid_panel"
  )
}

#' @export
print.sdid_panel <- function(x, ...) {
  d <- x$data
  cli::cli_h3("sdid_panel")
  cli::cli_bullets(c(
    "*" = "{length(unique(d[[x$unit]]))} unit{?s} ({.field {x$unit}}) x {length(unique(d[[x$time]]))} period{?s} ({.field {x$time}}, {min(d[[x$time]])}-{max(d[[x$time]])})",
    "*" = "outcome   {x$outcome %||% 'none by default'}",
    "*" = "treatment {if (is.null(x$treatment)) 'none' else x$treatment}",
    "*" = "columns   {.field {setdiff(names(d), c(x$unit, x$time))}}"
  ))
  invisible(x)
}

#' California's Proposition 99 as an sdid_panel
#'
#' The example panel from Abadie, Diamond and Hainmueller (2010) and
#' Arkhangelsky et al. (2021): per-capita cigarette sales in 39 US states,
#' 1970-2000, with California treated from 1989. The data ship with
#' `synthdid`; this wraps them for use with [sdid_spec()] and [fit_spec()].
#'
#' @return An [sdid_panel()].
#'
#' @examples
#' prop99_panel()
#'
#' @export
prop99_panel <- function() {
  env <- new.env()
  utils::data("california_prop99", package = "synthdid", envir = env)
  sdid_panel(env$california_prop99, unit = "State", time = "Year",
             outcome = "PacksPerCapita", treatment = "treated")
}

# Keep only units observed in every period, and report the rest. The
# treated unit must survive: dropping it would silently fit a different
# question.
balance_panel <- function(df, periods, treated, candidates = unique(df$unit)) {
  by_unit <- split(df$time[!is.na(df$y)], df$unit[!is.na(df$y)])
  # Candidates with no rows at all (a pool naming a unit the outcome does not
  # cover) are reported like any other gap, not silently lost.
  units <- unique(c(treated, candidates))
  n_obs <- vapply(units, function(u) length(intersect(by_unit[[u]] %||% numeric(0), periods)), integer(1))
  complete <- units[n_obs == length(periods)]
  if (!treated %in% complete) {
    cli::cli_abort(c(
      "The treated unit {.val {treated}} is missing the outcome in {length(periods) - n_obs[[treated]]} of {length(periods)} period{?s}.",
      i = "Narrow the window or drop those periods with {.field design$periods_drop}."
    ))
  }
  dropped <- units[n_obs < length(periods)]
  list(
    data = df[df$unit %in% complete & df$time %in% periods, , drop = FALSE],
    dropped = tibble::tibble(
      unit = dropped,
      reason = sprintf("missing %d of %d periods", length(periods) - n_obs[dropped], length(periods))
    )
  )
}
