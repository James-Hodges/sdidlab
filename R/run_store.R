#' Save, load and index fitted runs
#'
#' A run store is a directory of fitted runs, one folder per run, grouped by
#' specification hash. Each folder holds the specification as YAML, the
#' fitted object and placebo draws, the pool and balance records, and a
#' one-row summary that [refresh_runs_index()] collects into an index. A run
#' reloaded with [load_run()] can be re-plotted without refitting.
#'
#' @param run An `sdid_run` from [fit_spec()].
#' @param store Directory of the run store.
#' @param path A single run's directory (the `path` column of the index).
#' @param write Cache the index as `index.csv` in the store.
#' @param keep_latest_per_spec How many runs to keep for each specification.
#' @param dry_run Report what would be deleted without deleting it.
#'
#' @return `save_run()` returns the run with its `path` set; `load_run()` an
#'   `sdid_run`; `refresh_runs_index()` a tibble with one row per run;
#'   `prune_runs()` the rows it deleted (or would delete).
#'
#' @examples
#' \donttest{
#' store <- tempfile("runs")
#' run <- fit_spec(sdid_spec(inference = list(replications = 20)), prop99_panel(),
#'                 store = store, seed = 1)
#' idx <- refresh_runs_index(store)
#' idx[, c("spec_hash", "att", "se")]
#' again <- load_run(idx$path[1])
#' all.equal(again$fit$att, run$fit$att)
#' }
#' @name run_store
NULL

#' @rdname run_store
#' @export
save_run <- function(run, store) {
  dir <- file.path(store, substr(run$spec_hash, 1L, 8L), run$run_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_spec_yaml(run$spec, file.path(dir, "spec.yml"))
  light <- run
  light$fit["est"] <- list(NULL)
  light$fit["placebo_draws"] <- list(NULL)
  light$fit["placebo_gaps"] <- list(NULL)
  light$path <- NULL
  saveRDS(light, file.path(dir, "run.rds"))
  if (!is.null(run$fit$est)) saveRDS(run$fit$est, file.path(dir, "est.rds"))
  if (!is.null(run$fit$placebo_draws)) saveRDS(run$fit$placebo_draws, file.path(dir, "placebo.rds"))
  if (!is.null(run$fit$placebo_gaps)) saveRDS(run$fit$placebo_gaps, file.path(dir, "placebo_gaps.rds"))
  saveRDS(run_index_row(run), file.path(dir, "row.rds"))
  run$path <- dir
  invisible(run)
}

#' @rdname run_store
#' @export
load_run <- function(path) {
  file <- file.path(path, "run.rds")
  if (!file.exists(file)) {
    cli::cli_abort(c("No completed run at {.file {path}}.",
                     i = "It was cancelled or crashed before it was saved."))
  }
  run <- readRDS(file)
  if (file.exists(file.path(path, "est.rds"))) run$fit$est <- readRDS(file.path(path, "est.rds"))
  if (file.exists(file.path(path, "placebo.rds"))) run$fit$placebo_draws <- readRDS(file.path(path, "placebo.rds"))
  if (file.exists(file.path(path, "placebo_gaps.rds"))) run$fit$placebo_gaps <- readRDS(file.path(path, "placebo_gaps.rds"))
  run$path <- path
  run
}

#' @rdname run_store
#' @export
refresh_runs_index <- function(store, write = TRUE) {
  rows <- list.files(store, pattern = "^row[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(rows)) return(empty_runs_index())
  index <- do.call(rbind, lapply(rows, function(f) {
    r <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r$path <- dirname(f)
    r
  }))
  index <- tibble::as_tibble(index[order(index$created_at, decreasing = TRUE), , drop = FALSE])
  if (write) utils::write.csv(index, file.path(store, "index.csv"), row.names = FALSE)
  index
}

#' @rdname run_store
#' @export
prune_runs <- function(store, keep_latest_per_spec = 3L, dry_run = FALSE) {
  index <- refresh_runs_index(store, write = FALSE)
  if (!nrow(index)) return(invisible(index))
  rank <- stats::ave(seq_len(nrow(index)), index$spec_hash, FUN = seq_along)
  doomed <- index[rank > keep_latest_per_spec, , drop = FALSE]
  if (!nrow(doomed)) {
    cli::cli_alert_info("Nothing to prune: no specification has more than {keep_latest_per_spec} run{?s}.")
  } else if (dry_run) {
    cli::cli_alert_info("Would delete {nrow(doomed)} run{?s}.")
  } else {
    unlink(doomed$path, recursive = TRUE)
    refresh_runs_index(store)
    cli::cli_alert_success("Deleted {nrow(doomed)} run{?s}.")
  }
  invisible(doomed)
}

run_index_row <- function(run) {
  s <- run$spec; f <- run$fit
  tibble::tibble(
    run_id = run$run_id, spec_hash = run$spec_hash,
    spec_id = s$id, label = s$label, note = s$note, tags = paste(s$tags, collapse = ","),
    created_at = run$meta$created_at, seeded = run$meta$seeded, seed = run$meta$seed,
    outcome = s$data$outcome, treated_unit = run$treated_unit, treat_time = run$treat_time,
    window_start = min(run$periods), window_end = max(run$periods),
    estimator = s$estimator, pool = pool_text(s$pool),
    covariate_mode = s$covariates$mode, n_covariates = length(s$covariates$vars),
    se_method = s$inference$se_method, replications = s$inference$replications,
    trim_k = s$inference$trim_k, conf_level = s$inference$conf_level,
    candidates = sum(run$pool$included),
    dropped_for_balance = nrow(run$balance$dropped %||% data.frame()),
    N0 = f$N0, T0 = f$T0, n_periods = f$n_periods,
    att = f$att, se = f$se, t = f$t, ci_lo = f$ci_lo, ci_hi = f$ci_hi,
    trend_rmse = f$trend_rmse, level_gap = f$level_gap,
    placebos_kept = f$placebos_kept, placebos_total = f$placebos_total,
    total_secs = run$timings$total_secs,
    sdidlab_version = run$meta$sdidlab_version, synthdid_sha = run$meta$synthdid_sha
  )
}

empty_runs_index <- function() {
  tibble::tibble(run_id = character(0), spec_hash = character(0), label = character(0),
                 created_at = character(0), att = numeric(0), se = numeric(0),
                 ci_lo = numeric(0), ci_hi = numeric(0), path = character(0))
}
