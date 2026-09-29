# Changelog

## sdidlab 0.1.0

- First public release, extracted from a higher-education research
  pipeline and generalised to any long panel.
- Estimation core verified against synthdid: SDID, SC and DiD estimates,
  weights and untrimmed placebo standard errors match exactly, with and
  without covariates.
- Fixes found during verification: placebo inference for `sc` and `did`
  estimates; fit diagnostics and the trimming rule now use the
  covariate-adjusted outcome;
  [`project_out_covariates()`](https://james-hodges.github.io/sdidlab/reference/project_out_covariates.md)
  works with arrays; synthetic-control weights are solved exactly
  (active-set) rather than stopping short of convergence; exact
  pre-period fits are flagged.
- New: Abadie-L’Hour penalty in
  [`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md);
  [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md),
  [`register_outcome()`](https://james-hodges.github.io/sdidlab/reference/registry.md)
  and
  [`register_pool()`](https://james-hodges.github.io/sdidlab/reference/registry.md);
  self-contained scripts from
  [`write_spec_script()`](https://james-hodges.github.io/sdidlab/reference/write_spec_script.md);
  the interactive
  [`run_demo()`](https://james-hodges.github.io/sdidlab/reference/run_demo.md)
  app.
