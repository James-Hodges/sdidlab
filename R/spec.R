#' Default specification
#'
#' Every field of a specification, at its default. [sdid_spec()] merges
#' overrides onto this, and [spec_hash()] drops fields still at their
#' default before hashing, so adding a new field with a default never
#' changes an existing specification's identity.
#'
#' @return A named list.
#' @keywords internal
#' @export
spec_defaults <- function() {
  list(
    # Metadata. Never hashed -- renaming a spec must not orphan its runs.
    schema_version = 1L,
    id             = NA_character_,
    label          = NA_character_,
    note           = "",
    tags           = character(0),

    data = list(
      outcome       = NA_character_,   # a panel column or a registered outcome
      outcome_label = NA_character_    # axis label; defaults to the outcome name
    ),

    design = list(
      treated_unit = NA_character_,    # NA: take it from the panel's treatment column
      window       = list(start = NA_real_, end = NA_real_),  # NA: the panel's range
      periods_drop = numeric(0),       # periods removed from both match and effect
      treat_time   = NA_real_,         # first treated period; NA: from the panel
      # First period in which the outcome could reflect the treatment, for
      # outcomes that phase in. Periods from treat_time - anticipation up to
      # (not including) effect_time are dropped: contaminated but not yet
      # saturated, so they are uninterpretable from either side.
      effect_time  = NA_real_,
      anticipation = 0                 # moves the first treated period earlier
    ),

    pool = list(
      base            = "all",         # "all" or a registered pool
      units           = character(0),  # an explicit donor list, instead of base
      drop_units      = character(0),
      add_units       = character(0),
      min_donors      = 2L,
      exclusion_notes = list()         # why, keyed by unit; not hashed
    ),

    covariates = list(
      vars     = character(0),
      mode     = "none",               # none | optimized | projected
      z_score  = TRUE,
      log_vars = character(0)          # log1p-transformed, used as log_<name>
    ),

    estimator = "sdid",                # sdid | sc | did | scm

    # Synthetic control, used only when estimator == "scm".
    scm = list(
      normalize    = FALSE,            # rescale every unit to 100 at the last pre-period
      bias_correct = TRUE,
      v_method     = "equal",          # equal | mspe | synth
      v_starts     = 4L,
      v_maxit      = 300L,
      bc_donors    = "auto",           # auto | positive | all
      match_times  = numeric(0),       # pre-periods matched on; empty means all
      penalty      = 0                 # Abadie-L'Hour lambda
    ),

    inference = list(
      se_method    = "placebo",        # placebo | none
      replications = 200L,
      trim_k       = 5,
      conf_level   = 0.95
    ),

    # Free-form settings for registered outcomes and pools. Hashed.
    options = list()
  )
}

#' Build a specification
#'
#' A specification is a complete, declarative description of one model: the
#' outcome, the treated unit and timing, the donor pool, covariates, the
#' estimator and the inference. Specifications hash to a stable identifier
#' ([spec_hash()]), round-trip through YAML ([write_spec_yaml()]), and are
#' fitted by [fit_spec()].
#'
#' Pass any subset of sections; each merges onto [spec_defaults()], one level
#' deep (two for `design$window`), so `design = list(treat_time = 1989)`
#' leaves the other design fields at their defaults. Unknown fields are an
#' error rather than silently ignored.
#'
#' @param ... Named sections to override: `data`, `design`, `pool`,
#'   `covariates`, `estimator`, `scm`, `inference`, `options`, and the
#'   metadata fields `id`, `label`, `note`, `tags`.
#'
#' @return An `sdid_spec`.
#'
#' @examples
#' spec <- sdid_spec(
#'   label = "Prop 99, SDID",
#'   design = list(treat_time = 1989, window = list(start = 1970)),
#'   inference = list(replications = 100)
#' )
#' spec
#' spec_hash(spec)
#'
#' @export
sdid_spec <- function(...) {
  overrides <- list(...)
  spec <- merge_spec(spec_defaults(), overrides)
  spec <- coerce_spec(spec)
  validate_spec(spec)
  structure(spec, class = "sdid_spec", spec_hash = spec_hash_impl(spec))
}

#' @export
print.sdid_spec <- function(x, ...) {
  w <- x$design$window
  inf <- x$inference
  cli::cli_h3("sdid_spec {.emph {if (is.na(x$label)) '(unlabelled)' else x$label}}")
  cli::cli_bullets(c(
    "*" = "hash      {.strong {spec_hash(x)}}",
    "*" = "outcome   {if (is.na(x$data$outcome)) '(panel default)' else x$data$outcome}",
    "*" = "treated   {if (is.na(x$design$treated_unit)) 'from the panel' else x$design$treated_unit}, {spec_timing_text(x)}",
    "*" = "window    {fmt_na(w$start, 'start')}-{fmt_na(w$end, 'end')}",
    "*" = "pool      {pool_text(x$pool)}",
    "*" = "estimator {x$estimator}",
    "*" = "inference {inf$se_method}{if (inf$se_method == 'placebo') paste0(', ', inf$replications, ' reps, trim_k ', inf$trim_k) else ''}"
  ))
  invisible(x)
}

# Merge a partial override list over a base spec. One level deep for
# sections, two for design$window, so `design = list(window = list(end =
# 2016))` does not silently discard `start`.
merge_spec <- function(base, overrides) {
  if (!length(overrides)) return(base)
  # `list(design = a, design = b)` is legal R, and `overrides[["design"]]`
  # silently returns the first -- the fit would run and answer a question
  # nobody asked.
  dupes <- unique(names(overrides)[duplicated(names(overrides))])
  if (length(dupes)) {
    cli::cli_abort(c(
      "Section{?s} {.field {dupes}} given more than once.",
      i = "Merge them into a single argument: {.code {dupes[1]} = list(...)}."
    ))
  }
  unknown <- setdiff(names(overrides), names(base))
  if (length(unknown)) {
    cli::cli_abort(c(
      "Unknown spec field{?s}: {.field {unknown}}.",
      i = "Valid fields: {.field {names(base)}}."
    ))
  }
  free_form <- c("options")
  for (nm in names(overrides)) {
    val <- overrides[[nm]]
    if (nm %in% free_form) {
      base[[nm]] <- utils::modifyList(base[[nm]], as.list(val))
    } else if (is.list(base[[nm]]) && is.list(val) && !is.null(names(val))) {
      bad <- setdiff(names(val), names(base[[nm]]))
      if (length(bad)) {
        cli::cli_abort("Unknown {cli::qty(bad)}field{?s} in {.field {nm}}: {.field {bad}}.")
      }
      if (identical(nm, "design") && !is.null(val$window)) {
        val$window <- utils::modifyList(base[[nm]]$window, as.list(val$window))
      }
      # exclusion_notes is a map, so it is replaced rather than merged.
      notes <- if (identical(nm, "pool")) val$exclusion_notes else NULL
      base[[nm]] <- utils::modifyList(base[[nm]], val)
      if (!is.null(notes)) base[[nm]]$exclusion_notes <- as.list(notes)
    } else {
      base[[nm]] <- val
    }
  }
  base
}

# Canonical types. Without this, sdid_spec(design = list(treat_time = 1989))
# and treat_time = 1989L would hash differently, because default-dropping
# compares with identical().
coerce_spec <- function(spec) {
  chr1 <- function(x) if (is.null(x) || (length(x) == 1L && is.na(x))) NA_character_ else as.character(x)
  num1 <- function(x) if (is.null(x) || (length(x) == 1L && is.na(x))) NA_real_ else as.numeric(x)
  chrv <- function(x) sort(unique(as.character(unlist(x))))
  numv <- function(x) sort(unique(as.numeric(unlist(x))))

  spec$schema_version <- as.integer(spec$schema_version)
  spec$id    <- chr1(spec$id)
  spec$label <- chr1(spec$label)
  spec$note  <- as.character(spec$note %||% "")
  spec$tags  <- as.character(unlist(spec$tags))

  spec$data$outcome       <- chr1(spec$data$outcome)
  spec$data$outcome_label <- chr1(spec$data$outcome_label)

  d <- spec$design
  d$treated_unit <- chr1(d$treated_unit)
  d$window$start <- num1(d$window$start)
  d$window$end   <- num1(d$window$end)
  d$periods_drop <- numv(d$periods_drop)
  d$treat_time   <- num1(d$treat_time)
  # The NULL guard matters: `effect_time: ~` in YAML makes modifyList()
  # delete the element rather than set it to NA.
  d$effect_time  <- num1(d$effect_time)
  d$anticipation <- as.numeric(d$anticipation %||% 0)
  spec$design <- d

  p <- spec$pool
  p$base       <- as.character(p$base)
  p$units      <- chrv(p$units)
  p$drop_units <- chrv(p$drop_units)
  p$add_units  <- chrv(p$add_units)
  p$min_donors <- as.integer(p$min_donors)
  p$exclusion_notes <- as.list(p$exclusion_notes)
  spec$pool <- p

  cv <- spec$covariates
  cv$vars     <- as.character(unlist(cv$vars))
  cv$mode     <- as.character(cv$mode)
  cv$z_score  <- as.logical(cv$z_score)
  cv$log_vars <- as.character(unlist(cv$log_vars))
  spec$covariates <- cv

  spec$estimator <- as.character(spec$estimator)

  s <- spec$scm
  s$normalize    <- as.logical(s$normalize)
  s$bias_correct <- as.logical(s$bias_correct)
  s$v_method     <- as.character(s$v_method)
  s$v_starts     <- as.integer(s$v_starts)
  s$v_maxit      <- as.integer(s$v_maxit)
  s$bc_donors    <- as.character(s$bc_donors)
  s$match_times  <- numv(s$match_times)
  s$penalty      <- as.numeric(s$penalty)
  spec$scm <- s

  i <- spec$inference
  i$se_method    <- as.character(i$se_method)
  i$replications <- as.integer(i$replications)
  i$trim_k       <- as.numeric(i$trim_k)
  i$conf_level   <- as.numeric(i$conf_level)
  spec$inference <- i

  spec$options <- as.list(spec$options)
  spec
}

validate_spec <- function(spec) {
  one_of <- function(x, allowed, field) {
    if (length(x) != 1L || !x %in% allowed) {
      cli::cli_abort("{.field {field}} must be one of {.val {allowed}}; got {.val {x}}.")
    }
  }
  one_of(spec$estimator, c("sdid", "sc", "did", "scm"), "estimator")
  one_of(spec$covariates$mode, c("none", "optimized", "projected"), "covariates$mode")
  one_of(spec$inference$se_method, c("placebo", "none"), "inference$se_method")
  one_of(spec$scm$v_method, c("equal", "mspe", "synth"), "scm$v_method")
  one_of(spec$scm$bc_donors, c("auto", "positive", "all"), "scm$bc_donors")

  if (!identical(spec$covariates$mode, "none") && !length(spec$covariates$vars)) {
    cli::cli_abort("{.field covariates$mode} is {.val {spec$covariates$mode}} but no {.field covariates$vars} are named.")
  }
  if (identical(spec$estimator, "scm") && !identical(spec$covariates$mode, "none")) {
    cli::cli_abort("Time-varying covariates are not supported with {.code estimator = \"scm\"}.")
  }
  if (length(spec$covariates$log_vars) &&
      !all(paste0("log_", spec$covariates$log_vars) %in% spec$covariates$vars)) {
    cli::cli_abort("Each {.field covariates$log_vars} entry {.val x} must appear in {.field covariates$vars} as {.val log_x}.")
  }
  w <- spec$design$window
  if (!is.na(w$start) && !is.na(w$end) && w$start >= w$end) {
    cli::cli_abort("{.field design$window} start ({w$start}) must be before its end ({w$end}).")
  }
  tt <- spec$design$treat_time
  if (!is.na(tt) && ((!is.na(w$start) && tt <= w$start) || (!is.na(w$end) && tt > w$end))) {
    cli::cli_abort(c(
      "{.field design$treat_time} ({tt}) must fall inside the window, after its start.",
      i = "There must be at least one pre-period."
    ))
  }
  if (spec$design$anticipation < 0) {
    cli::cli_abort(c(
      "{.field design$anticipation} must be >= 0; got {spec$design$anticipation}.",
      i = "It moves the first treated period earlier, not later. To delay when the effect is measured, set {.field design$effect_time}."
    ))
  }
  et <- spec$design$effect_time
  if (!is.na(et) && !is.na(tt) && et <= tt) {
    cli::cli_abort("{.field design$effect_time} ({et}) must be strictly later than {.field design$treat_time} ({tt}).")
  }
  p <- spec$pool
  if (length(p$units) && !identical(p$base, "all")) {
    cli::cli_abort("Give either {.field pool$units} or a named {.field pool$base}, not both.")
  }
  both <- intersect(p$drop_units, p$add_units)
  if (length(both)) cli::cli_abort("Unit{?s} {.val {both}} appear{?s/} in both {.field drop_units} and {.field add_units}.")
  if (!is.na(spec$design$treated_unit) &&
      spec$design$treated_unit %in% c(p$add_units, p$units)) {
    cli::cli_abort("The treated unit cannot be in its own donor pool.")
  }
  if (p$min_donors < 2L) cli::cli_abort("{.field pool$min_donors} must be at least 2.")

  inf <- spec$inference
  if (inf$conf_level <= 0 || inf$conf_level >= 1) cli::cli_abort("{.field inference$conf_level} must be strictly between 0 and 1.")
  if (inf$replications < 2L) cli::cli_abort("{.field inference$replications} must be at least 2.")
  if (inf$trim_k <= 0) cli::cli_abort("{.field inference$trim_k} must be positive (use Inf for no trimming).")
  if (spec$scm$penalty < 0) cli::cli_abort("{.field scm$penalty} must be non-negative.")
  invisible(spec)
}

# Periods dropped between the first treated period and the first period in
# which the outcome could reflect it.
blackout_periods <- function(spec, times) {
  d <- spec$design
  if (is.na(d$effect_time) || is.na(d$treat_time)) return(times[0])
  times[times >= d$treat_time - d$anticipation & times < d$effect_time]
}

fmt_na <- function(x, what) if (is.na(x)) paste0("(panel ", what, ")") else format(x)

spec_timing_text <- function(spec) {
  d <- spec$design
  out <- if (is.na(d$treat_time)) "timing from the panel" else paste0("treated from ", format(d$treat_time))
  if (d$anticipation > 0) out <- paste0(out, ", anticipation ", format(d$anticipation))
  if (!is.na(d$effect_time)) out <- paste0(out, ", effect measured from ", format(d$effect_time))
  out
}

pool_text <- function(pool) {
  base <- if (length(pool$units)) paste(length(pool$units), "listed units") else pool$base
  mods <- c(
    if (length(pool$drop_units)) paste0("-", length(pool$drop_units), " by hand"),
    if (length(pool$add_units)) paste0("+", length(pool$add_units), " by hand")
  )
  if (!length(mods)) base else paste0(base, " (", paste(mods, collapse = ", "), ")")
}

spec_subtitle <- function(spec) {
  paste0(spec$label %|NA|% spec$estimator, " | pool: ", pool_text(spec$pool))
}

`%|NA|%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
