#' Fit one specification end to end
#'
#' Specification in, run out: builds the estimator inputs with
#' [build_setup()], fits, runs inference, and -- given a store -- saves the
#' whole run so it can be reloaded and re-plotted without refitting.
#'
#' # Reproducibility of the standard error
#'
#' Placebo standard errors are simulated, so they depend on the state of the
#' RNG stream when the fit begins. `seed = NULL` inherits the caller's
#' stream: reproducible for a whole sweep run after one `set.seed()`, but a
#' spec refitted on its own gets a different SE. Passing
#' `seed = seed_for(spec)` makes a spec reproducible in isolation. The
#' caller's own stream is left untouched either way.
#'
#' @param spec An [sdid_spec()].
#' @param panel An [sdid_panel()].
#' @param store Directory of a run store to save into, or `NULL`.
#' @param keep How much of the fit to retain: `"all"` (default) keeps the
#'   fitted object and every placebo draw, which every plot needs; `"est"`
#'   keeps only the fitted object; `"none"` keeps the summary alone.
#' @param seed `NULL` to inherit the caller's RNG stream, or an integer seed
#'   applied to this fit only.
#'
#' @return An `sdid_run`: `spec`, `spec_hash`, `run_id`, `fit`, `pool`,
#'   `balance`, `periods`, `timings` and `meta`.
#'
#' @examples
#' \donttest{
#' spec <- sdid_spec(label = "Prop 99", inference = list(replications = 50))
#' run <- fit_spec(spec, prop99_panel(), seed = 1)
#' run
#' }
#'
#' @export
fit_spec <- function(spec, panel, store = NULL, keep = c("all", "est", "none"),
                     seed = NULL) {
  keep <- match.arg(keep)
  started <- Sys.time()
  t_setup <- system.time(setup <- build_setup(spec, panel))[["elapsed"]]

  Y <- setup$Y
  X <- setup$X
  if (identical(spec$covariates$mode, "projected")) {
    Y <- project_out_covariates(Y, X, setup$N0, setup$T0)
    attr(Y, "beta") <- NULL
    X <- NULL
  }

  fit_once <- function() {
    if (identical(spec$estimator, "scm")) {
      return(fit_scm_spec(Y, setup$N0, setup$T0, setup$periods, spec, keep = keep))
    }
    fit_sdid(Y, setup$N0, setup$T0, X = X,
             replications = spec$inference$replications,
             conf_level   = spec$inference$conf_level,
             trim_k       = spec$inference$trim_k,
             se_method    = spec$inference$se_method,
             estimator    = spec$estimator,
             keep         = keep)
  }
  t_fit <- system.time(
    fit <- if (is.null(seed)) fit_once() else with_local_seed(seed, fit_once())
  )[["elapsed"]]
  fit$keep <- keep

  run <- structure(
    list(
      spec = spec, spec_hash = spec_hash(spec), run_id = new_run_id(spec),
      fit = fit, pool = setup$pool, balance = setup$balance,
      periods = setup$periods, blackout = setup$blackout,
      outcome = setup$outcome,
      treated_unit = setup$treated_unit, treated_label = setup$treated_label,
      treat_time = setup$treat_time,
      timings = list(setup_secs = t_setup, fit_secs = t_fit,
                     total_secs = as.numeric(difftime(Sys.time(), started, units = "secs"))),
      meta = run_meta(started, seed)
    ),
    class = "sdid_run"
  )
  if (!is.null(store)) run <- save_run(run, store)
  run
}

#' @export
print.sdid_run <- function(x, ...) {
  f <- x$fit
  ci <- if (is.na(f$se)) "no interval (se_method = none)" else sprintf("[%.3g, %.3g]", f$ci_lo, f$ci_hi)
  cli::cli_h3("sdid_run {.emph {x$spec$label %|NA|% x$spec_hash}}")
  cli::cli_bullets(c(
    "*" = "ATT       {signif(f$att, 4)}   {round(100 * x$spec$inference$conf_level)}% CI {ci}",
    "*" = "SE        {if (is.na(f$se)) 'NA' else sprintf('%.3g  (t = %.2f)', f$se, f$t)}",
    "*" = "fit       {x$spec$estimator}, N0 = {f$N0}, T0 = {f$T0}, pre-period RMSE {signif(f$trend_rmse, 3)}",
    "*" = "placebos  {if (is.na(f$placebos_kept)) 'none' else paste0(f$placebos_kept, '/', f$placebos_total, ' kept at trim_k ', f$trim_k)}",
    "*" = "treated   {x$treated_unit} from {x$treat_time}",
    "*" = "spec      {x$spec_hash}   run {x$run_id}"
  ))
  invisible(x)
}

#' Fit a list of specifications
#'
#' Fits each specification with [fit_spec()], in order. A specification that
#' cannot be fitted (too few donors, a too-short window) is reported and
#' skipped, and the rest carry on.
#'
#' @param specs A list of [sdid_spec()]s, for example from [read_specs()].
#' @inheritParams fit_spec
#' @param seed_each Seed each fit with [seed_for()], so every specification's
#'   standard error is reproducible on its own. `FALSE` inherits one stream
#'   across the whole sweep instead, so a result depends on its position.
#' @param progress Print one line per fitted specification.
#'
#' @return A list of `sdid_run`s, with `NULL` for specifications that failed.
#'
#' @examples
#' \donttest{
#' specs <- list(
#'   sdid_spec(label = "SDID", inference = list(replications = 30)),
#'   sdid_spec(label = "SC", estimator = "sc", inference = list(replications = 30))
#' )
#' runs <- run_specs(specs, prop99_panel())
#' sdid_results_table(runs)
#' }
#'
#' @export
run_specs <- function(specs, panel, store = NULL, keep = c("all", "est", "none"),
                      seed_each = TRUE, progress = TRUE) {
  keep <- match.arg(keep)
  total <- length(specs)
  failed <- character(0)
  runs <- lapply(seq_along(specs), function(i) {
    spec <- specs[[i]]
    label <- spec$label %|NA|% spec$id %|NA|% spec_hash(spec)
    tryCatch({
      run <- fit_spec(spec, panel, store = store, keep = keep,
                      seed = if (seed_each) seed_for(spec) else NULL)
      if (progress) {
        cli::cli_alert_success("{i}/{total} {label}: ATT {signif(run$fit$att, 4)}{if (is.na(run$fit$se)) '' else paste0(', SE ', signif(run$fit$se, 3))}")
      }
      run
    }, error = function(e) {
      failed[[length(failed) + 1L]] <<- label
      cli::cli_alert_danger("{i}/{total} {label}: {conditionMessage(e)}")
      NULL
    })
  })
  if (length(failed)) {
    cli::cli_alert_warning("{total - length(failed)}/{total} fitted; {length(failed)} failed: {.val {failed}}.")
  }
  runs
}

#' A specification's own seed
#'
#' Folds a specification's hash into a base seed, so refitting one
#' specification alone reproduces the standard error a sweep reported.
#'
#' @param spec An [sdid_spec()].
#' @param base Base seed to fold the hash into.
#'
#' @return A single integer.
#'
#' @examples
#' seed_for(sdid_spec())
#'
#' @export
seed_for <- function(spec, base = 20110901L) {
  digits <- strtoi(substr(spec_hash(spec), 1L, 7L), base = 16L)
  as.integer((base + digits) %% .Machine$integer.max)
}

# Seed this fit without disturbing the caller's stream.
with_local_seed <- function(seed, code) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  force(code)
}

# Must not touch the RNG: a sweep that inherits one stream would otherwise
# shift every later specification's placebo draws.
new_run_id <- function(spec) {
  stamp <- gsub("[.]", "", format(Sys.time(), "%Y%m%dT%H%M%OS6"))
  paste0(stamp, "-", substr(spec_hash(spec), 1L, 8L))
}

run_meta <- function(started, seed) {
  sd <- utils::packageDescription("synthdid")
  list(
    created_at = format(started, "%Y-%m-%dT%H:%M:%S"),
    seeded = !is.null(seed),
    seed = seed %||% NA_integer_,
    r_version = as.character(getRversion()),
    sdidlab_version = as.character(utils::packageVersion("sdidlab")),
    synthdid_version = sd$Version %||% NA_character_,
    synthdid_sha = sd$RemoteSha %||% NA_character_
  )
}

# The SCM path answers a different question and returns a different shape,
# so it goes through an adapter rather than a branch inside fit_sdid().
fit_scm_spec <- function(Y, N0, T0, periods, spec, keep = "all") {
  o <- spec$scm
  if (isTRUE(o$normalize)) {
    base <- Y[, T0]
    bad <- !is.finite(base) | base == 0
    if (bad[N0 + 1L]) {
      cli::cli_abort("The treated unit's last pre-period is {base[N0 + 1L]}; it cannot be normalised.")
    }
    if (any(bad)) {
      cli::cli_warn("Dropping {sum(bad)} unit{?s} with a zero or missing last pre-period before normalising.")
      Y <- Y[!bad, , drop = FALSE]
      N0 <- N0 - sum(bad[seq_len(N0)])
      base <- base[!bad]
    }
    Y <- 100 * Y / base
  }
  outcome_periods <- if (length(o$match_times)) {
    idx <- match(o$match_times, periods)
    idx <- idx[!is.na(idx) & idx <= T0]
    if (!length(idx)) {
      cli::cli_warn("No {.field scm$match_times} fall in the pre-period; matching on all of them.")
      NULL
    } else idx
  } else NULL

  args <- list(N0 = N0, T0 = T0, outcome_periods = outcome_periods,
               bias_correct = isTRUE(o$bias_correct), bc_donors = o$bc_donors,
               v_method = o$v_method, v_starts = o$v_starts, v_maxit = o$v_maxit,
               penalty = o$penalty)
  fit <- do.call(fit_scm, c(list(Y = Y), args))

  se <- p <- NA_real_; rank <- n_placebo <- NA_integer_; draws <- NULL; placebo_gaps <- NULL
  if (!identical(spec$inference$se_method, "none") && !identical(keep, "none")) {
    perm <- do.call(scm_placebo_permutation, c(list(Y = Y), args))
    r <- perm$results
    if (!any(!is.finite(r$rmspe_ratio_hybrid))) {
      pv <- rmspe_pvalue(stats::setNames(r$rmspe_ratio_hybrid, r$unit), r$att,
                         r$unit[r$treated], "two.sided")
      p <- pv$p_value; rank <- pv$rank
    } else {
      cli::cli_warn("Some units fit the pre-period exactly, so the RMSPE ranking is undefined; no p-value.")
    }
    n_placebo <- nrow(r)
    # The SD of the placebo effects: a scale for the interval, not a
    # sampling standard error.
    se <- stats::sd(r$att[!r$treated])
    draws <- tibble::tibble(draw = seq_len(nrow(r)), att = r$att, rmse = r$pre_rmspe_raw,
                            keep = TRUE, placebo_unit = r$unit, treated = r$treated)
    placebo_gaps <- list(corrected = perm$gaps, raw = perm$gaps_raw)
  }
  z <- stats::qnorm(1 - (1 - spec$inference$conf_level) / 2)
  out <- list(
    att = fit$att, se = se, t = fit$att / se,
    ci_lo = fit$att - z * se, ci_hi = fit$att + z * se,
    trend_rmse = fit$pre_rmspe_raw, level_gap = 0,
    N0 = N0, T0 = T0, n_periods = ncol(Y),
    omega = fit$omega, lambda = numeric(0), estimator = "scm",
    trim_k = NA_real_, placebos_kept = n_placebo, placebos_total = n_placebo,
    gap = fit$gap, gap_raw = fit$gap_raw, rmspe_p = p, rmspe_rank = rank,
    rmspe_ratio = fit$rmspe_ratio_hybrid, exact_fit = fit$exact_fit,
    exact_match = fit$exact_match,
    v = fit$v, n_predictors = fit$n_predictors
  )
  out["est"] <- list(NULL)
  out["placebo_draws"] <- list(draws)
  out["placebo_gaps"] <- list(if (identical(keep, "all")) placebo_gaps else NULL)
  out
}
