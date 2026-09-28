#' Build the estimator inputs for a specification
#'
#' Resolves a specification against a panel: the outcome series, the treated
#' unit and its timing, the donor pool, the periods (window, dropped
#' periods, blackout), a balanced panel, and the covariate array. The result
#' is what [fit_spec()] hands to the estimator, and is useful on its own for
#' checking what a specification actually selects.
#'
#' Every exclusion is recorded rather than silent: `pool` lists each
#' candidate donor with whether it was included and why not, and `balance`
#' lists donors dropped for missing the outcome in a kept period.
#'
#' @param spec An [sdid_spec()].
#' @param panel An [sdid_panel()].
#'
#' @return An `sdid_setup`: `Y`, `N0`, `T0` (as `synthdid::panel.matrices()`
#'   returns them), `X` (covariate array or `NULL`), `pool`, `balance`,
#'   `periods`, `blackout`, `treated_unit`, `treated_label`, `treat_time` and
#'   `spec`.
#'
#' @examples
#' setup <- build_setup(sdid_spec(), prop99_panel())
#' dim(setup$Y)
#' c(N0 = setup$N0, T0 = setup$T0)
#'
#' @export
build_setup <- function(spec, panel) {
  if (!inherits(panel, "sdid_panel")) cli::cli_abort("{.arg panel} must be an {.fn sdid_panel}.")
  timing <- resolve_timing(spec, panel)
  treated <- timing$treated_unit
  series <- outcome_series(spec, panel)

  pool <- resolve_pool(spec, panel, treated)
  donors <- pool$unit[pool$included]
  if (length(donors) < spec$pool$min_donors) {
    cli::cli_abort(c(
      "Only {length(donors)} donor{?s} in the pool; {.field pool$min_donors} is {spec$pool$min_donors}.",
      i = "{sum(!pool$included)} candidate{?s} {?was/were} excluded; see {.code build_setup(...)$pool}."
    ))
  }

  times <- sort(unique(series$time))
  w <- spec$design$window
  in_window <- times[times >= (w$start %|NA|% -Inf) & times <= (w$end %|NA|% Inf)]
  blackout <- blackout_periods(spec_with_timing(spec, timing), in_window)
  periods <- setdiff(in_window, c(spec$design$periods_drop, blackout))
  first_treated <- if (is.na(spec$design$effect_time)) {
    timing$treat_time - spec$design$anticipation
  } else {
    spec$design$effect_time
  }
  n_pre <- sum(periods < first_treated)
  n_post <- sum(periods >= first_treated)
  if (n_pre < 2L) {
    cli::cli_abort(c(
      "The design leaves {n_pre} pre-treatment period{?s}; at least 2 are needed.",
      i = if (length(blackout)) "Blacked out: {.val {blackout}}." else NULL
    ))
  }
  if (n_post < 1L) cli::cli_abort("The design leaves no post-treatment period.")

  df <- series[series$unit %in% c(treated, donors) & series$time %in% periods, , drop = FALSE]
  bal <- balance_panel(df, periods, treated, candidates = donors)
  df <- bal$data
  df$treated <- as.integer(df$unit == treated & df$time >= first_treated)

  m <- synthdid::panel.matrices(df, unit = "unit", time = "time",
                                outcome = "y", treatment = "treated")

  structure(
    list(
      Y = m$Y, N0 = m$N0, T0 = m$T0,
      X = build_covariate_array(spec, panel, m$Y),
      pool = pool, balance = list(dropped = bal$dropped),
      periods = as.numeric(colnames(m$Y)), blackout = blackout,
      treated_unit = treated, treated_label = unit_labels(panel, treated),
      treat_time = timing$treat_time, spec = spec
    ),
    class = "sdid_setup"
  )
}

# The treated unit and first treated period: from the spec, or else from the
# panel's treatment column.
resolve_timing <- function(spec, panel) {
  treated <- spec$design$treated_unit
  treat_time <- spec$design$treat_time
  if (is.na(treated) || is.na(treat_time)) {
    if (is.null(panel$treatment)) {
      cli::cli_abort(c(
        "The specification does not set {.field design$treated_unit} and {.field design$treat_time}.",
        i = "Set them, or declare a {.arg treatment} column in {.fn sdid_panel}."
      ))
    }
    d <- panel$data
    on <- d[[panel$treatment]] %in% 1
    ever <- unique(d[[panel$unit]][on])
    if (length(ever) != 1L) {
      cli::cli_abort("The panel's treatment column marks {length(ever)} treated unit{?s}; sdidlab needs exactly one.")
    }
    if (is.na(treated)) treated <- ever
    if (is.na(treat_time)) treat_time <- min(d[[panel$time]][on & d[[panel$unit]] == treated])
  }
  if (!treated %in% panel$data[[panel$unit]]) {
    cli::cli_abort("Treated unit {.val {treated}} is not in the panel.")
  }
  list(treated_unit = treated, treat_time = treat_time)
}

spec_with_timing <- function(spec, timing) {
  spec$design$treated_unit <- timing$treated_unit
  spec$design$treat_time <- timing$treat_time
  spec
}

# unit / time / y for the spec's outcome: a panel column or a registered
# outcome.
outcome_series <- function(spec, panel) {
  name <- spec$data$outcome %|NA|% panel$outcome
  if (is.null(name) || is.na(name)) {
    cli::cli_abort("No outcome: set {.field data$outcome} or give {.fn sdid_panel} a default {.arg outcome}.")
  }
  d <- panel$data
  if (name %in% names(d)) {
    if (!is.numeric(d[[name]])) cli::cli_abort("Outcome column {.field {name}} is not numeric.")
    return(data.frame(unit = d[[panel$unit]], time = d[[panel$time]], y = d[[name]],
                      stringsAsFactors = FALSE))
  }
  reg <- .registry$outcomes[[name]]
  if (is.null(reg)) {
    cli::cli_abort(c(
      "Outcome {.val {name}} is neither a panel column nor a registered outcome.",
      i = "Registered outcomes: {.val {registered_outcomes()}}."
    ))
  }
  out <- as.data.frame(reg$fn(panel, spec), stringsAsFactors = FALSE)
  if (!all(c("unit", "time", "y") %in% names(out))) {
    cli::cli_abort("Registered outcome {.val {name}} must return columns {.field unit}, {.field time} and {.field y}.")
  }
  out$unit <- as.character(out$unit)
  out[, c("unit", "time", "y")]
}

#' Resolve a specification's donor pool
#'
#' Every candidate donor with whether it is included and, if not, why.
#' The candidates are the panel's other units (`pool$base = "all"`), the
#' units listed in `pool$units`, or the result of a registered pool; then
#' `drop_units` are excluded and `add_units` forced in.
#'
#' @inheritParams build_setup
#' @param treated The treated unit, which is never a candidate.
#'
#' @return A tibble with `unit`, `label`, `included` and `reason`.
#'
#' @examples
#' pool <- resolve_pool(sdid_spec(pool = list(drop_units = "Nevada")),
#'                      prop99_panel(), treated = "California")
#' pool[!pool$included, ]
#'
#' @export
resolve_pool <- function(spec, panel, treated) {
  p <- spec$pool
  all_units <- sort(unique(panel$data[[panel$unit]]))
  base <- if (length(p$units)) {
    missing <- setdiff(p$units, all_units)
    if (length(missing)) cli::cli_abort("{.field pool$units} names unit{?s} not in the panel: {.val {missing}}.")
    tibble::tibble(unit = p$units, included = TRUE, reason = NA_character_)
  } else if (identical(p$base, "all")) {
    tibble::tibble(unit = all_units, included = TRUE, reason = NA_character_)
  } else {
    fn <- .registry$pools[[p$base]]
    if (is.null(fn)) {
      cli::cli_abort(c("Pool {.val {p$base}} is not registered.",
                       i = "Registered pools: {.val {registered_pools()}}."))
    }
    res <- fn(panel, spec)
    if (is.data.frame(res)) {
      res <- tibble::as_tibble(res)
      res$unit <- as.character(res$unit)
      if (is.null(res$reason)) res$reason <- ifelse(res$included, NA_character_, "excluded by pool")
      res[, c("unit", "included", "reason")]
    } else {
      tibble::tibble(unit = as.character(res), included = TRUE, reason = NA_character_)
    }
  }
  base <- base[base$unit != treated, , drop = FALSE]

  drop_hit <- base$unit %in% p$drop_units & base$included
  base$included[drop_hit] <- FALSE
  base$reason[drop_hit] <- vapply(base$unit[drop_hit], function(u) {
    note <- p$exclusion_notes[[u]]
    if (is.null(note)) "dropped by hand" else paste("dropped by hand:", note)
  }, character(1))

  add <- setdiff(p$add_units, base$unit[base$included])
  if (length(add)) {
    missing <- setdiff(add, all_units)
    if (length(missing)) cli::cli_abort("{.field pool$add_units} names unit{?s} not in the panel: {.val {missing}}.")
    base <- base[!base$unit %in% add, , drop = FALSE]
    base <- rbind(base, tibble::tibble(unit = add, included = TRUE, reason = NA_character_))
  }
  base$label <- unit_labels(panel, base$unit)
  base[, c("unit", "label", "included", "reason")]
}

unit_labels <- function(panel, units) {
  if (is.null(panel$label)) return(units)
  d <- panel$data
  lookup <- stats::setNames(as.character(d[[panel$label]]), d[[panel$unit]])
  out <- unname(lookup[units])
  ifelse(is.na(out), units, out)
}

# Covariates as a units x periods x covariates array aligned with Y. Log
# transforms are applied first; z-scoring uses the estimation sample itself
# (the units and periods in Y), not the whole panel.
build_covariate_array <- function(spec, panel, Y) {
  cv <- spec$covariates
  if (!length(cv$vars) || identical(cv$mode, "none")) return(NULL)
  d <- panel$data
  for (v in cv$log_vars) {
    if (!v %in% names(d)) cli::cli_abort("No column {.field {v}} to log.")
    d[[paste0("log_", v)]] <- log1p(pmax(as.numeric(d[[v]]), 0))
  }
  missing <- setdiff(cv$vars, names(d))
  if (length(missing)) cli::cli_abort("Unknown covariate{?s}: {.field {missing}}.")

  units <- rownames(Y); times <- colnames(Y)
  key <- paste(d[[panel$unit]], d[[panel$time]], sep = "\r")
  cells <- outer(units, times, paste, sep = "\r")
  arr <- array(NA_real_, dim = c(length(units), length(times), length(cv$vars)),
               dimnames = list(units, times, cv$vars))
  for (k in seq_along(cv$vars)) {
    vals <- as.numeric(d[[cv$vars[k]]])[match(cells, key)]
    if (isTRUE(cv$z_score)) {
      s <- stats::sd(vals, na.rm = TRUE)
      vals <- if (is.finite(s) && s > 0) (vals - mean(vals, na.rm = TRUE)) / s else vals * 0
    }
    arr[, , k] <- matrix(vals, nrow = length(units))
  }
  if (anyNA(arr)) {
    cli::cli_abort(c(
      "Covariate array has missing values after balancing.",
      i = "Every covariate must be observed for every retained unit-period.",
      i = "Missing in: {.field {cv$vars[apply(arr, 3, anyNA)]}}."
    ))
  }
  arr
}
