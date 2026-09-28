#' Summarise fitted runs as a table
#'
#' One row per run: what was estimated, the estimate, its interval, and the
#' fit and inference records. `NULL` entries (specifications that failed in
#' [run_specs()]) are skipped.
#'
#' @param runs A list of `sdid_run`s.
#' @param digits Rounding for the estimate columns.
#'
#' @return A data frame.
#'
#' @examples
#' \donttest{
#' runs <- run_specs(list(sdid_spec(label = "SDID", inference = list(replications = 20))),
#'                   prop99_panel(), progress = FALSE)
#' sdid_results_table(runs)
#' }
#'
#' @export
sdid_results_table <- function(runs, digits = 3) {
  runs <- Filter(Negate(is.null), runs)
  if (!length(runs)) return(NULL)
  rows <- lapply(runs, function(run) {
    s <- run$spec; f <- run$fit
    data.frame(
      label = s$label %|NA|% s$id %|NA|% run$spec_hash,
      estimator = s$estimator,
      outcome = s$data$outcome %|NA|% "(default)",
      treated = run$treated_unit, treat_time = run$treat_time,
      window = paste0(min(run$periods), "-", max(run$periods)),
      N0 = f$N0, T0 = f$T0,
      ATT = signif(f$att, digits + 1), SE = signif(f$se, digits),
      t = round(f$t, 2), ci_lo = signif(f$ci_lo, digits + 1), ci_hi = signif(f$ci_hi, digits + 1),
      trend_rmse = signif(f$trend_rmse, digits),
      placebos = if (is.na(f$placebos_total)) NA_character_ else paste0(f$placebos_kept, "/", f$placebos_total),
      spec_hash = run$spec_hash,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Donor weights, largest first
#'
#' @param run An `sdid_run`.
#'
#' @return A tibble of `unit`, `label`, `weight` and `cumulative` weight.
#'
#' @examples
#' \donttest{
#' run <- fit_spec(sdid_spec(inference = list(se_method = "none")), prop99_panel())
#' head(donor_weight_table(run))
#' }
#'
#' @export
donor_weight_table <- function(run) {
  omega <- run$fit$omega
  if (is.null(omega) || !length(omega)) cli::cli_abort("This run has no unit weights.")
  out <- tibble::tibble(unit = names(omega), weight = as.numeric(omega))
  out <- out[order(out$weight, decreasing = TRUE), , drop = FALSE]
  out$label <- run_labels(run, out$unit)
  out$cumulative <- cumsum(out$weight)
  out[, c("unit", "label", "weight", "cumulative")]
}

#' What each donor contributes to the estimate and to the standard error
#'
#' Sets each donor's share of the counterfactual (its unit weight) against
#' its share of the placebo variance behind the standard error -- the draws
#' in which it stood in for the treated unit. A donor far above the diagonal
#' supplies more uncertainty than estimate.
#'
#' Dropping a donor because it inflates the standard error is selection on
#' the result; this table is for understanding an interval, not for
#' choosing a pool.
#'
#' @param run An `sdid_run` fitted with `keep = "all"`.
#'
#' @return A tibble with one row per donor: `unit`, `label`, `weight`,
#'   `weight_share`, `n_placebo`, `n_trimmed`, `att_mean`, `abs_att_max`,
#'   `rmse_ratio` and `se_share`, largest `se_share` first.
#'
#' @export
donor_influence <- function(run) {
  require_run_part(run, "placebo_draws", "placebo draws")
  if (identical(run$spec$estimator, "scm")) {
    cli::cli_abort("Donor influence is defined for the simulated placebo SE, not SCM permutation inference.")
  }
  draws <- run$fit$placebo_draws
  w <- donor_weight_table(run)
  kept <- draws[draws$keep %in% TRUE, , drop = FALSE]
  centre <- mean(kept$att)
  ss_total <- sum((kept$att - centre)^2)
  rows <- lapply(w$unit, function(id) {
    mine <- draws[draws$placebo_unit %in% id, , drop = FALSE]
    mine_kept <- mine[mine$keep %in% TRUE, , drop = FALSE]
    tibble::tibble(
      unit = id,
      n_placebo = nrow(mine),
      n_trimmed = sum(!mine$keep %in% TRUE),
      att_mean = if (nrow(mine)) mean(mine$att, na.rm = TRUE) else NA_real_,
      abs_att_max = if (nrow(mine)) max(abs(mine$att), na.rm = TRUE) else NA_real_,
      rmse_ratio = if (nrow(mine)) mean(mine$rmse, na.rm = TRUE) / run$fit$treated_rmse else NA_real_,
      se_share = if (ss_total > 0 && nrow(mine_kept)) sum((mine_kept$att - centre)^2) / ss_total else 0
    )
  })
  out <- do.call(rbind, rows)
  out$label <- w$label
  out$weight <- w$weight
  out$weight_share <- w$weight / sum(w$weight)
  out <- out[order(out$se_share, decreasing = TRUE), , drop = FALSE]
  out[, c("unit", "label", "weight", "weight_share", "n_placebo", "n_trimmed",
          "att_mean", "abs_att_max", "rmse_ratio", "se_share")]
}

#' Donors carrying almost no weight
#'
#' Splits the donors at a weight threshold -- useful for reporting how
#' concentrated a counterfactual is, or for a robustness refit on the
#' material donors only (decided in advance, not after seeing the result).
#'
#' @param run An `sdid_run`.
#' @param threshold Weight below which a donor counts as negligible.
#'
#' @return A list: `drop` and `keep` (unit identifiers), their counts, the
#'   share of total weight below the threshold, and the full `table`.
#'
#' @export
low_weight_units <- function(run, threshold = 0.01) {
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    cli::cli_abort("{.arg threshold} must be a single number.")
  }
  w <- donor_weight_table(run)
  low <- w$weight < threshold
  list(drop = w$unit[low], keep = w$unit[!low], n_drop = sum(low), n_keep = sum(!low),
       weight_dropped = sum(w$weight[low]) / sum(w$weight), threshold = threshold, table = w)
}

run_labels <- function(run, units) {
  pool <- run$pool
  lookup <- if (!is.null(pool) && "label" %in% names(pool)) stats::setNames(pool$label, pool$unit) else character(0)
  if (!is.null(run$treated_label)) lookup[run$treated_unit] <- run$treated_label
  out <- unname(lookup[units])
  ifelse(is.na(out), units, out)
}
