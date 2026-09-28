#' Register derived outcomes and named donor pools
#'
#' Extension points for data that a plain column cannot express. A
#' specification names an outcome in `data$outcome` and a pool in
#' `pool$base`; when the name is not a column (or `"all"`), sdidlab looks
#' it up here.
#'
#' * An **outcome** function takes `(panel, spec)` and returns a data frame
#'   with columns `unit`, `time` and `y` -- for example a share computed
#'   from two columns, or a series built from another table.
#' * A **pool** function takes `(panel, spec)` and returns either a
#'   character vector of donor units, or a data frame with columns `unit`,
#'   `included` (logical) and optionally `reason`, which is how a pool
#'   records *why* a candidate was left out.
#'
#' Both receive the whole specification, so a spec's free-form `options`
#' section can parameterise them. Registrations last for the R session;
#' put them in a package's `.onLoad()` or at the top of a script.
#'
#' @param name The name specifications will use.
#' @param fn The function, as described above.
#' @param label,units For outcomes: how axes and tables should describe it.
#'
#' @return The name, invisibly.
#'
#' @examples
#' register_outcome("packs_log", function(panel, spec) {
#'   d <- panel$data
#'   data.frame(unit = d[[panel$unit]], time = d[[panel$time]],
#'              y = log(d$PacksPerCapita))
#' }, label = "Log packs per capita")
#'
#' register_pool("west", function(panel, spec) {
#'   c("Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico",
#'     "Oregon", "Utah", "Washington", "Wyoming")
#' })
#'
#' registered_outcomes()
#' registered_pools()
#' @name registry
NULL

.registry <- new.env(parent = emptyenv())
.registry$outcomes <- list()
.registry$pools <- list()

#' @rdname registry
#' @export
register_outcome <- function(name, fn, label = name, units = NA_character_) {
  check_registration(name, fn)
  .registry$outcomes[[name]] <- list(fn = fn, label = label, units = units)
  invisible(name)
}

#' @rdname registry
#' @export
register_pool <- function(name, fn) {
  check_registration(name, fn)
  if (identical(name, "all")) cli::cli_abort("{.val all} is reserved for the whole panel.")
  .registry$pools[[name]] <- fn
  invisible(name)
}

#' @rdname registry
#' @export
registered_outcomes <- function() names(.registry$outcomes)

#' @rdname registry
#' @export
registered_pools <- function() names(.registry$pools)

check_registration <- function(name, fn) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    cli::cli_abort("{.arg name} must be a single non-empty string.")
  }
  if (!is.function(fn)) cli::cli_abort("{.arg fn} must be a function.")
}
