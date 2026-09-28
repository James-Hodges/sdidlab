# sdidlab (development version)

* First public release, extracted from a higher-education research pipeline
  and generalised to any long panel.
* Estimation core verified against synthdid: SDID, SC and DiD estimates,
  weights and untrimmed placebo standard errors match exactly, with and
  without covariates.
* Fixes found during verification: placebo inference for `sc` and `did`
  estimates; fit diagnostics and the trimming rule now use the
  covariate-adjusted outcome; `project_out_covariates()` works with arrays;
  synthetic-control weights are solved exactly (active-set) rather than
  stopping short of convergence; exact pre-period fits are flagged.
* New: Abadie-L'Hour penalty in `fit_scm()`; `sdid_panel()`,
  `register_outcome()` and `register_pool()`; self-contained scripts from
  `write_spec_script()`; the interactive `run_demo()` app.
