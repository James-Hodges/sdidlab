# Package index

## Data and specifications

- [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md)
  : Declare a long panel for sdidlab
- [`prop99_panel()`](https://james-hodges.github.io/sdidlab/reference/prop99_panel.md)
  : California's Proposition 99 as an sdid_panel
- [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md)
  : Build a specification
- [`spec_hash()`](https://james-hodges.github.io/sdidlab/reference/spec_hash.md)
  : A specification's stable identifier
- [`write_spec_yaml()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md)
  [`read_spec_yaml()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md)
  [`read_specs()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md)
  [`write_specs_yaml()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md)
  : Read and write specifications as YAML
- [`build_setup()`](https://james-hodges.github.io/sdidlab/reference/build_setup.md)
  : Build the estimator inputs for a specification
- [`resolve_pool()`](https://james-hodges.github.io/sdidlab/reference/resolve_pool.md)
  : Resolve a specification's donor pool
- [`register_outcome()`](https://james-hodges.github.io/sdidlab/reference/registry.md)
  [`register_pool()`](https://james-hodges.github.io/sdidlab/reference/registry.md)
  [`registered_outcomes()`](https://james-hodges.github.io/sdidlab/reference/registry.md)
  [`registered_pools()`](https://james-hodges.github.io/sdidlab/reference/registry.md)
  : Register derived outcomes and named donor pools

## Fitting

- [`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md)
  : Fit one specification end to end
- [`run_specs()`](https://james-hodges.github.io/sdidlab/reference/run_specs.md)
  : Fit a list of specifications
- [`seed_for()`](https://james-hodges.github.io/sdidlab/reference/seed_for.md)
  : A specification's own seed
- [`save_run()`](https://james-hodges.github.io/sdidlab/reference/run_store.md)
  [`load_run()`](https://james-hodges.github.io/sdidlab/reference/run_store.md)
  [`refresh_runs_index()`](https://james-hodges.github.io/sdidlab/reference/run_store.md)
  [`prune_runs()`](https://james-hodges.github.io/sdidlab/reference/run_store.md)
  : Save, load and index fitted runs
- [`write_spec_script()`](https://james-hodges.github.io/sdidlab/reference/write_spec_script.md)
  : Write a runnable R script for a set of specifications

## Estimators and inference

- [`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
  : Fit one synthetic-DiD specification and return its full result set
- [`placebo_se_trimmed()`](https://james-hodges.github.io/sdidlab/reference/placebo_se_trimmed.md)
  : Placebo standard error with poorly-fitting placebos trimmed
- [`adjusted_outcome()`](https://james-hodges.github.io/sdidlab/reference/adjusted_outcome.md)
  : The outcome a synthdid estimate was actually fitted on
- [`project_out_covariates()`](https://james-hodges.github.io/sdidlab/reference/project_out_covariates.md)
  : Residualise the outcome on covariates using untreated cells only
- [`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
  : Fit a bias-corrected synthetic control
- [`sc_weights_simplex()`](https://james-hodges.github.io/sdidlab/reference/sc_weights_simplex.md)
  : Simplex weights minimising a V-weighted predictor discrepancy
- [`scm_placebo_permutation()`](https://james-hodges.github.io/sdidlab/reference/scm_placebo_permutation.md)
  : Permute treatment in space across every donor
- [`rmspe_pvalue()`](https://james-hodges.github.io/sdidlab/reference/rmspe_pvalue.md)
  : Ranked-RMSPE permutation p-value for a synthetic control

## Results and plots

- [`sdid_results_table()`](https://james-hodges.github.io/sdidlab/reference/sdid_results_table.md)
  : Summarise fitted runs as a table
- [`donor_weight_table()`](https://james-hodges.github.io/sdidlab/reference/donor_weight_table.md)
  : Donor weights, largest first
- [`donor_influence()`](https://james-hodges.github.io/sdidlab/reference/donor_influence.md)
  : What each donor contributes to the estimate and to the standard
  error
- [`low_weight_units()`](https://james-hodges.github.io/sdidlab/reference/low_weight_units.md)
  : Donors carrying almost no weight
- [`plot_sdid_trajectory()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_gap_path()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_event_study()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_donor_weights()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_lambda_weights()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_sdid_units()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_placebo_distribution()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_placebo_trim()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_donor_influence()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  [`plot_scm_gaps()`](https://james-hodges.github.io/sdidlab/reference/plots.md)
  : Plot a fitted run
- [`plot_att_forest()`](https://james-hodges.github.io/sdidlab/reference/plot_att_forest.md)
  : Forest plot of several estimates
- [`theme_sdid()`](https://james-hodges.github.io/sdidlab/reference/theme_sdid.md)
  : The sdidlab ggplot2 theme

## Demo

- [`run_demo()`](https://james-hodges.github.io/sdidlab/reference/run_demo.md)
  : Launch the interactive demo

## Internals

- [`spec_defaults()`](https://james-hodges.github.io/sdidlab/reference/spec_defaults.md)
  : Default specification
- [`van_der_corput()`](https://james-hodges.github.io/sdidlab/reference/van_der_corput.md)
  : A deterministic low-discrepancy point in \[0, 1)
- [`sdidlab`](https://james-hodges.github.io/sdidlab/reference/sdidlab-package.md)
  [`sdidlab-package`](https://james-hodges.github.io/sdidlab/reference/sdidlab-package.md)
  : sdidlab: Specify, Run, Store and Plot Synthetic
  Difference-in-Differences Models
