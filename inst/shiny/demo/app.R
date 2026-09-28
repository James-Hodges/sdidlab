# sdidlab demo: build a synthetic difference-in-differences specification,
# fit it, inspect every diagnostic, compare specifications, and download the
# specification as YAML plus an R script that reproduces the fit exactly.
#
# Launch locally with sdidlab::run_demo().

library(shiny)
library(bslib)
library(sdidlab)

MAX_UPLOAD_MB <- 5
MAX_UNITS     <- 200
MAX_REPS      <- 100
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)

prop99 <- prop99_panel()

plot_tabs <- list(
  "Trajectory"           = function(run) plot_sdid_trajectory(run),
  "Gap path"             = function(run) plot_gap_path(run),
  "Event study"          = function(run) plot_event_study(run),
  "Donor weights"        = function(run) plot_donor_weights(run),
  "Time weights"         = function(run) plot_lambda_weights(run),
  "Placebo distribution" = function(run) plot_placebo_distribution(run),
  "Trim cutoff"          = function(run) plot_placebo_trim(run),
  "Donor influence"      = function(run) plot_donor_influence(run),
  "SCM placebo gaps"     = function(run) plot_scm_gaps(run)
)
# Which plots make sense for which estimator.
plots_for <- function(estimator) {
  switch(estimator,
    sdid = setdiff(names(plot_tabs), "SCM placebo gaps"),
    sc   = setdiff(names(plot_tabs), c("SCM placebo gaps", "Time weights")),
    did  = c("Trajectory", "Gap path", "Event study", "Placebo distribution", "Trim cutoff"),
    scm  = c("SCM placebo gaps", "Gap path", "Event study", "Donor weights", "Placebo distribution")
  )
}

ui <- page_navbar(
  title = "sdidlab demo",
  theme = bs_theme(bootswatch = "flatly"),
  # Pages scroll rather than squeezing every card into the viewport: a plot
  # compressed to a sliver is unreadable, and synthdid's arrows cannot be
  # drawn at all below a minimum height.
  fillable = FALSE,
  sidebar = sidebar(
    width = 340,
    accordion(
      open = c("Data", "Design"),
      accordion_panel(
        "Data",
        radioButtons("source", NULL, c("California Proposition 99" = "prop99", "Upload a CSV" = "csv")),
        conditionalPanel(
          "input.source == 'csv'",
          fileInput("file", NULL, accept = ".csv", placeholder = sprintf("Long format, up to %d MB", MAX_UPLOAD_MB)),
          selectInput("col_unit", "Unit column", NULL),
          selectInput("col_time", "Time column (numeric)", NULL),
          selectInput("col_outcome", "Outcome column", NULL),
          selectInput("col_label", "Display-name column (optional)", NULL)
        )
      ),
      accordion_panel(
        "Design",
        selectInput("treated", "Treated unit", NULL),
        numericInput("treat_time", "First treated period", NA),
        sliderInput("window", "Window", min = 0, max = 1, value = c(0, 1), step = 1, sep = ""),
        selectizeInput("drop_units", "Exclude donors", NULL, multiple = TRUE)
      ),
      accordion_panel(
        "Estimator",
        radioButtons("estimator", NULL, inline = TRUE,
                     c("SDID" = "sdid", "SC" = "sc", "DiD" = "did", "Classic SCM" = "scm")),
        conditionalPanel(
          "input.estimator == 'scm'",
          radioButtons("v_method", "Predictor weights", c("Equal" = "equal", "MSPE search" = "mspe"), inline = TRUE),
          checkboxInput("bias_correct", "Bias correction (Abadie-L'Hour)", TRUE),
          numericInput("penalty", "Penalty (Abadie-L'Hour; 0 = classic)", 0, min = 0, step = 0.001)
        ),
        conditionalPanel(
          "input.estimator != 'scm'",
          selectizeInput("covars", "Covariates (numeric columns)", NULL, multiple = TRUE),
          radioButtons("cov_mode", "Covariate route", inline = TRUE,
                       c("synthdid (optimized)" = "optimized", "Residualize first (projected)" = "projected"))
        )
      ),
      accordion_panel(
        "Inference",
        conditionalPanel(
          "input.estimator != 'scm'",
          sliderInput("reps", "Placebo replications", 20, MAX_REPS, 50, step = 10),
          numericInput("trim_k", "Trim placebos fitting worse than k x treated (blank = no trim)", 5, min = 1)
        ),
        conditionalPanel("input.estimator == 'scm'",
                         helpText("SCM uses exact in-space permutation: every donor refitted as if treated.")),
        textInput("label", "Label for this fit", "")
      )
    ),
    actionButton("fit", "Fit", class = "btn-primary w-100"),
    helpText("Fits are cached by specification hash: refitting an identical specification is instant.")
  ),

  nav_panel(
    "Result",
    uiOutput("summary"),
    navset_card_underline(id = "plot_tab", !!!lapply(names(plot_tabs), function(nm) {
      nav_panel(nm, plotOutput(paste0("p_", make.names(nm)), height = "520px"))
    })),
    card(card_header("Specification (YAML)"), verbatimTextOutput("spec_yaml")),
    card(card_header("Excluded candidates"), tableOutput("excluded"))
  ),

  nav_panel(
    "Compare",
    card(
      card_header("Every fit this session"),
      plotOutput("forest", height = "420px"),
      tableOutput("results"),
      div(class = "d-flex gap-2",
          downloadButton("dl_yaml", "Download specifications (YAML)"),
          downloadButton("dl_script", "Download R script"),
          actionButton("clear", "Clear", class = "btn-outline-danger"))
    )
  ),

  nav_panel(
    "About",
    card(
      markdown(paste(
        "This demo runs [sdidlab](https://github.com/James-Hodges/sdidlab), an R package built around",
        "the `synthdid` estimator of Arkhangelsky, Athey, Hirshberg, Imbens and Wager (2021).",
        "",
        "* **SDID, SC and DiD** estimates match `synthdid` exactly; the untrimmed placebo standard error is",
        "  `synthdid`'s own, draw for draw.",
        "* **Trimmed placebo inference** drops placebos whose pre-period fit is more than *k* times worse than",
        "  the treated unit's. It trims the inference distribution only, never the donor pool.",
        "* **Classic SCM** adds predictor matching, the Abadie-L'Hour bias correction and penalty, and",
        "  ranked-RMSPE permutation inference.",
        "",
        "Every fit can be downloaded as a YAML specification and an R script that reproduces it.",
        sprintf("Uploads are limited to %d MB and %d units, and placebo replications to %d, to keep the",
                MAX_UPLOAD_MB, MAX_UNITS, MAX_REPS),
        "hosted demo responsive; the package itself has no such limits."
      ))
    )
  )
)

server <- function(input, output, session) {
  uploaded <- reactive({
    req(input$source == "csv", input$file)
    d <- tryCatch(utils::read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE),
                  error = function(e) NULL)
    validate(need(is.data.frame(d) && nrow(d) > 0, "Could not read that file as a CSV."))
    d
  })

  observeEvent(uploaded(), {
    cols <- names(uploaded())
    nums <- cols[vapply(uploaded(), is.numeric, logical(1))]
    updateSelectInput(session, "col_unit", choices = cols, selected = cols[1])
    updateSelectInput(session, "col_time", choices = nums, selected = nums[1])
    updateSelectInput(session, "col_outcome", choices = nums, selected = nums[min(2, length(nums))])
    updateSelectInput(session, "col_label", choices = c("(none)" = "", cols))
  })

  panel <- reactive({
    if (input$source == "prop99") return(prop99)
    d <- uploaded()
    req(input$col_unit, input$col_time, input$col_outcome)
    validate(need(length(unique(d[[input$col_unit]])) <= MAX_UNITS,
                  sprintf("The demo is limited to %d units.", MAX_UNITS)))
    tryCatch(
      sdid_panel(d, unit = input$col_unit, time = input$col_time, outcome = input$col_outcome,
                 label = if (nzchar(input$col_label %||% "")) input$col_label else NULL),
      error = function(e) validate(conditionMessage(e))
    )
  })

  # Reset the design controls when the data change.
  observeEvent(panel(), {
    p <- panel()
    d <- p$data
    units <- sort(unique(d[[p$unit]]))
    times <- sort(unique(d[[p$time]]))
    treated <- if (!is.null(p$treatment)) unique(d[[p$unit]][d[[p$treatment]] %in% 1])[1] else units[1]
    tt <- if (!is.null(p$treatment)) min(d[[p$time]][d[[p$treatment]] %in% 1]) else times[ceiling(length(times) * 2 / 3)]
    updateSelectInput(session, "treated", choices = units, selected = treated)
    updateNumericInput(session, "treat_time", value = tt)
    updateSliderInput(session, "window", min = min(times), max = max(times), value = range(times))
    updateSelectizeInput(session, "drop_units", choices = units, selected = character(0))
    covs <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], c(p$time, p$outcome, p$treatment))
    updateSelectizeInput(session, "covars", choices = covs, selected = character(0))
  })

  current_spec <- reactive({
    trim <- if (is.na(input$trim_k %||% NA)) Inf else input$trim_k
    covars <- if (input$estimator == "scm") character(0) else input$covars %||% character(0)
    sdid_spec(
      label = if (nzchar(input$label)) input$label else NA,
      design = list(treated_unit = input$treated, treat_time = input$treat_time,
                    window = list(start = input$window[1], end = input$window[2])),
      pool = list(drop_units = setdiff(input$drop_units %||% character(0), input$treated)),
      covariates = list(vars = covars, mode = if (length(covars)) input$cov_mode else "none"),
      estimator = input$estimator,
      scm = list(v_method = input$v_method, bias_correct = isTRUE(input$bias_correct),
                 penalty = max(0, input$penalty %||% 0)),
      inference = list(replications = min(input$reps, MAX_REPS), trim_k = trim)
    )
  })

  cache <- reactiveVal(list())
  fits <- reactiveVal(list())
  current_run <- reactiveVal(NULL)

  observeEvent(input$fit, {
    spec <- tryCatch(current_spec(), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 8)
      NULL
    })
    req(spec)
    key <- paste(spec_hash(spec), input$source, digest::digest(panel()$data))
    run <- cache()[[key]]
    if (is.null(run)) {
      run <- withProgress(message = "Fitting", detail = "placebo inference refits the model many times", {
        tryCatch(fit_spec(spec, panel(), seed = seed_for(spec)), error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 10)
          NULL
        })
      })
      req(run)
      cache(utils::modifyList(cache(), stats::setNames(list(run), key)))
    }
    if (is.na(spec$label)) spec$label <- sprintf("Fit %d (%s)", length(fits()) + 1, spec$estimator)
    run$spec$label <- spec$label
    current_run(run)
    fits(c(fits(), list(run)))
  })

  observeEvent(input$clear, fits(list()))

  output$summary <- renderUI({
    run <- current_run()
    if (is.null(run)) {
      return(card(card_body("Choose a design in the sidebar and press Fit. The built-in data reproduce the ",
                            "published Proposition 99 estimate of -15.6 packs per capita.")))
    }
    f <- run$fit
    ci <- if (is.na(f$se)) "no interval" else sprintf("[%s, %s]", signif(f$ci_lo, 3), signif(f$ci_hi, 3))
    warn <- if (isTRUE(f$exact_match) && run$spec$scm$penalty == 0) {
      div(class = "alert alert-warning mt-2",
          "The synthetic control matches its predictors exactly, so its weights (and the estimate) are not unique. ",
          "Set a small penalty to select the unique Abadie-L'Hour solution.")
    }
    tagList(
      layout_column_wrap(
        width = 1 / 4,
        value_box("ATT", signif(f$att, 4), paste0(round(100 * run$spec$inference$conf_level), "% CI ", ci)),
        value_box("Standard error", if (is.na(f$se)) "-" else signif(f$se, 3),
                  if (identical(run$spec$estimator, "scm")) sprintf("permutation p = %s", signif(f$rmspe_p, 3))
                  else sprintf("%d of %d placebos kept", f$placebos_kept, f$placebos_total)),
        value_box("Pre-period fit", signif(f$trend_rmse, 3), "trend RMSE"),
        value_box("Donors", f$N0, sprintf("%d pre-periods, %s treated from %s", f$T0, run$treated_label, run$treat_time))
      ),
      warn
    )
  })

  lapply(names(plot_tabs), function(nm) {
    id <- paste0("p_", make.names(nm))
    # Re-execute on resize rather than replaying a recording: synthdid's
    # curved arrows recorded at a tab's zero-size initial layout cannot be
    # replayed at full size.
    output[[id]] <- renderPlot({
      req((session$clientData[[paste0("output_", id, "_width")]] %||% 0) > 200,
          (session$clientData[[paste0("output_", id, "_height")]] %||% 0) > 200)
      run <- current_run()
      req(run)
      validate(need(nm %in% plots_for(run$spec$estimator),
                    sprintf("Not applicable to the %s estimator.", run$spec$estimator)))
      tryCatch(plot_tabs[[nm]](run), error = function(e) validate(conditionMessage(e)))
    }, res = 96, execOnResize = TRUE)
  })

  output$spec_yaml <- renderText({
    run <- current_run()
    req(run)
    write_spec_yaml(run$spec, NULL)
  })

  output$excluded <- renderTable({
    run <- current_run()
    req(run)
    pool <- run$pool[!run$pool$included, c("label", "reason")]
    bal <- run$balance$dropped
    rbind(pool, if (nrow(bal)) data.frame(label = bal$unit, reason = bal$reason))
  })

  output$results <- renderTable({
    req(length(fits()))
    tab <- sdid_results_table(fits())
    tab[, c("label", "estimator", "window", "N0", "T0", "ATT", "SE", "ci_lo", "ci_hi", "trend_rmse", "placebos")]
  })

  output$forest <- renderPlot({
    req(length(fits()))
    tab <- sdid_results_table(fits())
    plot_att_forest(data.frame(label = make.unique(tab$label), att = tab$ATT, ci_lo = tab$ci_lo, ci_hi = tab$ci_hi))
  }, res = 96, execOnResize = TRUE)

  session_specs <- reactive({
    specs <- lapply(fits(), `[[`, "spec")
    specs[!duplicated(vapply(specs, function(s) paste(s$label, spec_hash(s)), ""))]
  })

  data_code <- reactive({
    if (input$source == "prop99") return("panel <- prop99_panel()")
    lab <- if (nzchar(input$col_label %||% "")) sprintf(", label = %s", deparse(input$col_label)) else ""
    sprintf("panel <- sdid_panel(read.csv(%s), unit = %s, time = %s, outcome = %s%s)",
            deparse(input$file$name %||% "your_data.csv"), deparse(input$col_unit),
            deparse(input$col_time), deparse(input$col_outcome), lab)
  })

  output$dl_yaml <- downloadHandler(
    filename = function() "sdidlab_specs.yml",
    content = function(file) {
      validate(need(length(session_specs()), "Fit something first."))
      write_specs_yaml(session_specs(), file, title = "Specifications from the sdidlab demo")
    }
  )
  output$dl_script <- downloadHandler(
    filename = function() "sdidlab_analysis.R",
    content = function(file) {
      validate(need(length(session_specs()), "Fit something first."))
      write_spec_script(session_specs(), file, data_code = data_code(),
                        title = "Analysis exported from the sdidlab demo")
    }
  )
}

shinyApp(ui, server)
