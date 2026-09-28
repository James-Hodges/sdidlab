#' Plot a fitted run
#'
#' Every plot takes an `sdid_run` from [fit_spec()] (or [load_run()]) and
#' returns a `ggplot` object, so it can be themed, saved or combined further.
#'
#' * `plot_sdid_trajectory()` -- treated unit and synthetic control over time,
#'   with the time weights along the bottom (synthdid's own plot, restyled).
#' * `plot_gap_path()` -- treated minus synthetic in every period, measured
#'   against the baseline the estimator subtracts, so the post-period points
#'   average exactly to the ATT.
#' * `plot_event_study()` -- the effect in each post-treatment period.
#' * `plot_donor_weights()`, `plot_lambda_weights()` -- unit and time weights.
#' * `plot_sdid_units()` -- each donor's contribution (synthdid's units plot).
#' * `plot_placebo_distribution()`, `plot_placebo_trim()` -- the placebo
#'   reference distribution behind the standard error, and where trimming
#'   cut it.
#' * `plot_donor_influence()` -- each donor's share of the estimate against
#'   its share of the standard error (see [donor_influence()]).
#' * `plot_scm_gaps()` -- for `estimator = "scm"`: the treated unit's gap
#'   against every placebo's, the classic permutation figure.
#'
#' Plots that need the fitted object or the placebo draws require a run
#' fitted with `keep = "all"` (the default).
#'
#' @param run An `sdid_run`.
#' @param top_n Donors to show.
#' @param material Weight marked as the materiality threshold.
#' @param bins Histogram bins.
#' @param label_n Donors to label.
#' @param se_method Passed to `synthdid::synthdid_plot()`.
#' @param corrected For `plot_scm_gaps()`: plot bias-corrected gaps instead
#'   of the uncorrected ones.
#'
#' @return A `ggplot`.
#'
#' @examples
#' \donttest{
#' run <- fit_spec(sdid_spec(inference = list(replications = 50)), prop99_panel(), seed = 1)
#' plot_sdid_trajectory(run)
#' plot_gap_path(run)
#' plot_placebo_distribution(run)
#' }
#' @name plots
NULL

#' @rdname plots
#' @export
plot_sdid_trajectory <- function(run, se_method = "none") {
  require_run_part(run, "est", "fitted object")
  pal <- sdid_palette()
  p <- withCallingHandlers(
    synthdid::synthdid_plot(named_est(run), se.method = se_method,
                            treated.name = run$treated_label %||% "treated",
                            control.name = "synthetic control"),
    warning = function(w) if (grepl("deprecated", conditionMessage(w))) invokeRestart("muffleWarning")
  )
  p <- p +
    ggplot2::labs(
      title = "Treated unit and synthetic control",
      subtitle = spec_subtitle(run$spec), x = NULL, y = outcome_axis_label(run),
      caption = paste0(
        "The band along the bottom shows the time weights (lambda): the pre-periods the estimator leans on.\n",
        "SDID carries unit fixed effects, so a constant level gap between the two series is absorbed by design."
      )
    ) +
    ggplot2::geom_vline(xintercept = first_treated(run) - 0.5, colour = pal$synthetic,
                        linetype = "dashed", linewidth = 0.5) +
    theme_sdid()
  if (length(run$blackout)) {
    p <- p + ggplot2::annotate("rect", xmin = min(run$blackout) - 0.5, xmax = max(run$blackout) + 0.5,
                               ymin = -Inf, ymax = Inf, fill = pal$neutral, alpha = 0.25)
  }
  suppressMessages(p + ggplot2::scale_colour_manual(
    values = stats::setNames(c(pal$treated, pal$synthetic),
                             c(run$treated_label %||% "treated", "synthetic control"))
  ))
}

#' @rdname plots
#' @export
plot_gap_path <- function(run) {
  pal <- sdid_palette()
  periods <- run$periods
  T0 <- run$fit$T0
  raw <- NULL
  if (identical(run$spec$estimator, "scm")) {
    gap <- as.numeric(run$fit$gap)
    raw <- as.numeric(run$fit$gap_raw)
    baseline <- 0
    caption <- paste0("Solid: the bias-corrected gap, whose post-period average is the ATT (dotted). Periods used as predictors are pinned near zero by the correction.\n",
                      "Dashed grey: the uncorrected gap, which shows the synthetic control's actual pre-period fit.")
  } else {
    require_run_part(run, "est", "fitted object")
    est <- run$fit$est
    Y <- adjusted_outcome(est)
    N0 <- attr(est, "setup")$N0
    w <- attr(est, "weights")
    gap <- as.numeric(Y[N0 + 1L, ]) - as.numeric(crossprod(Y[seq_len(N0), , drop = FALSE], w$omega))
    # The baseline the estimator itself subtracts: the lambda-weighted
    # pre-period gap. Measured against it, the post-period points average
    # exactly to the ATT.
    baseline <- sum(w$lambda * gap[seq_len(T0)])
    caption <- paste0(
      "Measured against the lambda-weighted pre-period gap, the baseline ", run$spec$estimator,
      " subtracts, so the post-period points average to the ATT (dotted line).\n",
      "Wander in the shaded pre-period is the parallel-trends question; trend RMSE summarises it."
    )
  }
  df <- tibble::tibble(period = periods, gap = gap - baseline)
  missing <- setdiff(seq(min(periods), max(periods)), periods)
  if (length(missing) && all(periods == round(periods))) {
    df <- rbind(df, tibble::tibble(period = missing, gap = NA_real_))
    df <- df[order(df$period), , drop = FALSE]
  }
  last_pre <- periods[T0]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$gap)) +
    ggplot2::annotate("rect", xmin = min(periods) - 0.5, xmax = last_pre + 0.5,
                      ymin = -Inf, ymax = Inf, fill = pal$grid, alpha = 0.55) +
    ggplot2::geom_hline(yintercept = 0, colour = pal$zero, linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = run$fit$att, colour = pal$synthetic,
                        linetype = "dotted", linewidth = 0.6) +
    ggplot2::geom_line(colour = pal$treated, linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(colour = pal$treated, size = 2, na.rm = TRUE) +
    ggplot2::labs(
      title = "Treated minus synthetic",
      subtitle = sprintf("%s | ATT %s | pre-period trend RMSE %s", spec_subtitle(run$spec),
                         signif(run$fit$att, 4), signif(run$fit$trend_rmse, 3)),
      x = NULL, y = paste0("Gap: ", outcome_axis_label(run)), caption = caption
    ) +
    theme_sdid()
  if (!is.null(raw)) {
    p <- p + ggplot2::geom_line(data = tibble::tibble(period = periods, gap = raw),
                                colour = pal$neutral, linetype = "dashed", linewidth = 0.6)
  }
  p
}

#' @rdname plots
#' @export
plot_event_study <- function(run) {
  pal <- sdid_palette()
  post <- run$periods[-seq_len(run$fit$T0)]
  effect <- if (identical(run$spec$estimator, "scm")) {
    as.numeric(run$fit$gap)[-seq_len(run$fit$T0)]
  } else {
    require_run_part(run, "est", "fitted object")
    as.numeric(synthdid::synthdid_effect_curve(run$fit$est))
  }
  df <- tibble::tibble(period = post, effect = effect)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$effect)) +
    ggplot2::geom_hline(yintercept = 0, colour = pal$zero, linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = run$fit$att, colour = pal$synthetic,
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(colour = pal$treated, linewidth = 0.9) +
    ggplot2::geom_point(colour = pal$treated, size = 2.4) +
    ggplot2::labs(
      title = "Effect by post-treatment period",
      subtitle = sprintf("%s | average %s (dashed)", spec_subtitle(run$spec), signif(run$fit$att, 4)),
      x = NULL, y = outcome_axis_label(run),
      caption = "Per-period effects share one reference distribution; they are not separately tested."
    ) +
    theme_sdid()
}

#' @rdname plots
#' @export
plot_donor_weights <- function(run, top_n = 20L, material = 0.01) {
  pal <- sdid_palette()
  w <- donor_weight_table(run)
  shown <- utils::head(w, top_n)
  shown$label <- factor(shown$label, levels = rev(shown$label))
  ggplot2::ggplot(shown, ggplot2::aes(x = .data$weight, y = .data$label)) +
    ggplot2::geom_col(fill = pal$treated, width = 0.7) +
    ggplot2::geom_vline(xintercept = material, colour = pal$zero, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$weight)),
                       hjust = -0.15, size = 3, colour = pal$text_secondary) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = "Donor weights",
      subtitle = sprintf("%s | %d of %d donors above %s | top %d carry %.0f%%",
                         spec_subtitle(run$spec), sum(w$weight >= material), nrow(w),
                         format(material), nrow(shown), 100 * sum(shown$weight)),
      x = "Unit weight (omega)", y = NULL,
      caption = "Dashed: the materiality threshold. SDID regularises its weights, so on a large pool many donors sit just above zero."
    ) +
    theme_sdid()
}

#' @rdname plots
#' @export
plot_lambda_weights <- function(run) {
  pal <- sdid_palette()
  lambda <- run$fit$lambda
  if (!length(lambda) || all(lambda == 0)) {
    cli::cli_abort("This run has no time weights ({run$spec$estimator} does not use them).")
  }
  df <- tibble::tibble(period = run$periods[seq_along(lambda)], weight = as.numeric(lambda))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$period, y = .data$weight)) +
    ggplot2::geom_col(fill = pal$treated, width = 0.7) +
    ggplot2::labs(
      title = "Time weights",
      subtitle = sprintf("%s | heaviest pre-period: %s (%.2f)", spec_subtitle(run$spec),
                         format(df$period[which.max(df$weight)]), max(df$weight)),
      x = NULL, y = "Time weight (lambda)",
      caption = "Pre-treatment periods only: these are fitted to predict the post-period average."
    ) +
    theme_sdid()
}

#' @rdname plots
#' @export
plot_sdid_units <- function(run) {
  require_run_part(run, "est", "fitted object")
  p <- withCallingHandlers(
    synthdid::synthdid_units_plot(named_est(run)),
    warning = function(w) if (grepl("deprecated", conditionMessage(w))) invokeRestart("muffleWarning")
  )
  p + ggplot2::labs(title = "Donor contributions", subtitle = spec_subtitle(run$spec),
                    y = outcome_axis_label(run), x = NULL) +
    theme_sdid() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
}

#' @rdname plots
#' @export
plot_placebo_distribution <- function(run, bins = 50) {
  require_run_part(run, "placebo_draws", "placebo draws")
  pal <- sdid_palette()
  f <- run$fit
  draws <- f$placebo_draws
  if (identical(run$spec$estimator, "scm")) draws <- draws[!draws$treated, , drop = FALSE]
  kept <- draws[draws$keep, , drop = FALSE]
  z <- stats::qnorm(1 - (1 - run$spec$inference$conf_level) / 2)
  ggplot2::ggplot(draws, ggplot2::aes(x = .data$att)) +
    ggplot2::geom_histogram(bins = bins, fill = pal$neutral, alpha = 0.45) +
    ggplot2::geom_histogram(data = kept, bins = bins, fill = pal$treated) +
    ggplot2::geom_vline(xintercept = c(-z, z) * f$se, linetype = "dashed", colour = pal$zero, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = f$att, colour = pal$synthetic, linewidth = 1) +
    ggplot2::annotate("label", x = f$att, y = Inf, vjust = 1.2, label = sprintf("ATT %s", signif(f$att, 4)),
                      colour = pal$text_primary, fill = "white", size = 3.2) +
    ggplot2::labs(
      title = "Placebo distribution",
      subtitle = sprintf("%s | %d of %d draws kept | SE %s", spec_subtitle(run$spec),
                         f$placebos_kept, f$placebos_total, signif(f$se, 3)),
      x = paste0("Placebo estimate: ", outcome_axis_label(run)), y = "Draws",
      caption = paste0("Blue: draws kept; grey: trimmed for poor pre-period fit. Dashed: +/-", round(z, 2),
                       " SE. Orange: the estimate.\nPlacebo SEs are simulated and move a few percent between runs.")
    ) +
    theme_sdid()
}

#' @rdname plots
#' @export
plot_placebo_trim <- function(run, bins = 50) {
  require_run_part(run, "placebo_draws", "placebo draws")
  if (identical(run$spec$estimator, "scm")) cli::cli_abort("SCM permutation inference does not trim.")
  pal <- sdid_palette()
  f <- run$fit
  draws <- f$placebo_draws
  cutoff <- f$trim_k * f$treated_rmse
  p <- ggplot2::ggplot(draws, ggplot2::aes(x = .data$rmse, fill = .data$keep)) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = pal$treated, `FALSE` = pal$neutral),
                               labels = c(`TRUE` = "kept", `FALSE` = "trimmed"), breaks = c("TRUE", "FALSE")) +
    ggplot2::geom_vline(xintercept = f$treated_rmse, colour = pal$synthetic, linewidth = 1) +
    ggplot2::labs(
      title = "Placebo fit, and where the trim cut",
      subtitle = sprintf("%s | %d draw(s) trimmed at trim_k %s", spec_subtitle(run$spec),
                         sum(!draws$keep), format(f$trim_k)),
      x = "Pre-period trend RMSE", y = "Draws",
      caption = paste0("Orange: the treated unit's own pre-period fit. Dashed: the cutoff, trim_k times that.\n",
                       "Trimming restricts the inference distribution only; the counterfactual still uses every donor.")
    ) +
    theme_sdid()
  if (is.finite(cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = cutoff, colour = pal$zero, linetype = "dashed", linewidth = 0.6)
  }
  p
}

#' @rdname plots
#' @export
plot_donor_influence <- function(run, label_n = 6L) {
  pal <- sdid_palette()
  d <- donor_influence(run)
  net <- d$se_share - d$weight_share
  notable <- unique(c(utils::head(order(-net), label_n), utils::head(order(-d$weight_share), 2L)))
  d$lab <- ifelse(seq_len(nrow(d)) %in% notable, d$label, NA_character_)
  lim <- max(d$se_share, d$weight_share, na.rm = TRUE) * 1.08
  ggplot2::ggplot(d, ggplot2::aes(x = .data$weight_share, y = .data$se_share)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = pal$zero, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$n_trimmed), shape = 21, size = 3.2,
                        colour = "white", stroke = 0.6) +
    ggplot2::scale_fill_gradient(low = pal$treated, high = pal$neutral, name = "draws trimmed") +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), na.rm = TRUE, hjust = -0.12, size = 2.9,
                       colour = pal$text_secondary) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, lim)) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, lim)) +
    ggplot2::labs(
      title = "What each donor supplies, and what it costs",
      subtitle = spec_subtitle(run$spec),
      x = "Share of the counterfactual (unit weight)", y = "Share of the placebo variance",
      caption = paste0("Above the dashed line a donor supplies more uncertainty than estimate.\n",
                       "Dropping donors because the standard error falls is selection on the result; judge pools on pre-treatment fit, decided in advance.")
    ) +
    theme_sdid()
}

#' @rdname plots
#' @export
plot_scm_gaps <- function(run, corrected = FALSE) {
  if (!identical(run$spec$estimator, "scm")) cli::cli_abort("{.fn plot_scm_gaps} is for {.code estimator = \"scm\"} runs.")
  require_run_part(run, "placebo_gaps", "placebo gap paths")
  pal <- sdid_palette()
  # Uncorrected by default, as in Abadie, Diamond and Hainmueller: the bias
  # correction pins every predictor period to zero, which would flatten the
  # pre-period and hide exactly the fit this figure exists to show.
  g <- run$fit$placebo_gaps[[if (corrected) "corrected" else "raw"]]
  treated_row <- nrow(g)
  long <- tibble::tibble(
    unit = rep(rownames(g), times = ncol(g)),
    period = rep(run$periods, each = nrow(g)),
    gap = as.vector(g),
    treated = rep(seq_len(nrow(g)) == treated_row, times = ncol(g))
  )
  p_txt <- if (is.na(run$fit$rmspe_p)) "p undefined" else sprintf("rank %d of %d, p = %.3f",
                                                                 run$fit$rmspe_rank, nrow(g), run$fit$rmspe_p)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$period, y = .data$gap, group = .data$unit)) +
    ggplot2::geom_vline(xintercept = first_treated(run) - 0.5, colour = pal$zero, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0, colour = pal$zero, linewidth = 0.4) +
    ggplot2::geom_line(data = long[!long$treated, ], colour = pal$neutral, alpha = 0.5, linewidth = 0.4) +
    ggplot2::geom_line(data = long[long$treated, ], colour = pal$treated, linewidth = 1.1) +
    ggplot2::labs(
      title = "Treated gap against every placebo",
      subtitle = sprintf("%s | %s (post/pre RMSPE)", spec_subtitle(run$spec), p_txt),
      x = NULL, y = paste0("Gap: ", outcome_axis_label(run)),
      caption = paste0(
        "Blue: the treated unit. Grey: each donor refitted as if it were treated, with the treated unit excluded from its pool.\n",
        if (corrected) "Bias-corrected gaps: periods used as predictors are pinned near zero by construction."
        else "Uncorrected gaps, so the pre-period shows each synthetic control's actual fit. The p-value ranks corrected post-period gaps over this pre-period fit."
      )
    ) +
    theme_sdid()
}

#' Forest plot of several estimates
#'
#' @param df A data frame with `label`, `att`, `ci_lo` and `ci_hi` columns,
#'   such as [sdid_results_table()] output renamed, or built by hand.
#' @param highlight Labels to draw in the highlight colour, e.g. the
#'   pre-specified primary estimate.
#' @param sort_by `"none"` (input order), `"att"` or `"label"`.
#' @param xlab Axis label.
#'
#' @return A `ggplot`.
#'
#' @examples
#' plot_att_forest(data.frame(label = c("SDID", "SC", "DiD"), att = c(-15.6, -19.6, -27.3),
#'                            ci_lo = c(-35, -40, -45), ci_hi = c(4, 1, -9)),
#'                 highlight = "SDID")
#'
#' @export
plot_att_forest <- function(df, highlight = NULL, sort_by = c("none", "att", "label"),
                            xlab = "ATT") {
  sort_by <- match.arg(sort_by)
  missing <- setdiff(c("label", "att", "ci_lo", "ci_hi"), names(df))
  if (length(missing)) cli::cli_abort("{.arg df} is missing column{?s} {.field {missing}}.")
  pal <- sdid_palette()
  df <- tibble::as_tibble(df)
  df$is_highlight <- df$label %in% highlight
  levels <- switch(sort_by,
    none = rev(unique(as.character(df$label))),
    att = unique(as.character(df$label)[order(df$att)]),
    label = sort(unique(as.character(df$label)), decreasing = TRUE))
  df$label <- factor(df$label, levels = levels)
  n <- nrow(df)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$att, y = .data$label, colour = .data$is_highlight)) +
    ggplot2::geom_vline(xintercept = 0, colour = pal$zero, linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$ci_lo, xmax = .data$ci_hi), linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = pal$treated, `TRUE` = pal$synthetic),
                                 labels = c(`FALSE` = "other", `TRUE` = "highlighted"),
                                 guide = if (any(df$is_highlight)) "legend" else "none") +
    ggplot2::labs(
      title = "Estimates and intervals",
      subtitle = sprintf("%d specification%s; at alpha 0.05, P(at least one spurious rejection) ~ %d%%",
                         n, if (n == 1L) "" else "s", round(100 * (1 - 0.95^n))),
      x = xlab, y = NULL,
      caption = "When nothing is significant, the interval bound is the result."
    ) +
    theme_sdid()
}

#' The sdidlab ggplot2 theme
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @examples
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_sdid()
#' @export
theme_sdid <- function(base_size = 12) {
  pal <- sdid_palette()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.1),
      plot.subtitle = ggplot2::element_text(colour = pal$text_secondary, size = base_size * 0.85),
      plot.caption = ggplot2::element_text(colour = pal$text_muted, size = base_size * 0.72, hjust = 0),
      axis.title = ggplot2::element_text(colour = pal$text_secondary, size = base_size * 0.85),
      axis.text = ggplot2::element_text(colour = pal$text_secondary),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = pal$grid, linewidth = 0.3),
      legend.position = "top", legend.title = ggplot2::element_blank(),
      plot.title.position = "plot", plot.caption.position = "plot"
    )
}

sdid_palette <- function() {
  list(treated = "#2a78d6", synthetic = "#eb6834", neutral = "#8a8a85", grid = "#e3e3e0",
       zero = "#52514e", text_primary = "#0b0b0b", text_secondary = "#52514e", text_muted = "#8a8a85")
}

require_run_part <- function(run, part, what) {
  if (!is.null(run$fit[[part]])) return(invisible(TRUE))
  cli::cli_abort(c(
    "This run has no {what}.",
    i = "It was fitted with {.code keep = \"{run$fit$keep %||% 'none'}\"}; refit with {.code keep = \"all\"} to draw it."
  ), call = sys.call(-1))
}

outcome_axis_label <- function(run) {
  s <- run$spec
  if (!is.na(s$data$outcome_label)) return(s$data$outcome_label)
  name <- run$outcome %||% s$data$outcome
  if (is.null(name) || is.na(name)) return("Outcome")
  reg <- .registry$outcomes[[name]]
  if (is.null(reg)) return(name)
  if (is.na(reg$units)) reg$label else paste0(reg$label, " (", reg$units, ")")
}

first_treated <- function(run) {
  run$periods[run$fit$T0 + 1L]
}

# Relabel the fitted object's units with display names for synthdid's plots.
named_est <- function(run) {
  est <- run$fit$est
  setup <- attr(est, "setup")
  labels <- run_labels(run, rownames(setup$Y))
  rownames(setup$Y) <- make.unique(labels)
  attr(est, "setup") <- setup
  est
}
