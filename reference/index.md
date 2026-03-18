# Package index

## Data Synthesis

Functions for synthesising delay distribution data

- [`create_scenario()`](https://cm401.github.io/ddsynth/reference/create_scenario.md)
  : Create a scenario with specific characteristics
- [`generate_hierarchical_data_mixed()`](https://cm401.github.io/ddsynth/reference/generate_hierarchical_data_mixed.md)
  : Generate data from hierarchical model with mixed summary types and
  sample sizes
- [`generate_scenario_library()`](https://cm401.github.io/ddsynth/reference/generate_scenario_library.md)
  : Generate a comprehensive set of scenarios
- [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md)
  : Prepare Stan data from a list of dataset summaries

## Model Fitting

Functions for fitting hierarchical Bayesian models

- [`fit_model()`](https://cm401.github.io/ddsynth/reference/fit_model.md)
  : Fit Stan model to simulated data
- [`run_simulation_study_generalized()`](https://cm401.github.io/ddsynth/reference/run_simulation_study_generalized.md)
  : Run simulation study with generalized scenarios
- [`bsl_bridge_estimate()`](https://cm401.github.io/ddsynth/reference/bsl_bridge_estimate.md)
  : Compute a log-mean-exp bridge estimate of the marginal likelihood
- [`bsl_create_model()`](https://cm401.github.io/ddsynth/reference/bsl_create_model.md)
  : Create a BSL model object
- [`bsl_get_summary_stats()`](https://cm401.github.io/ddsynth/reference/bsl_get_summary_stats.md)
  : Compute predictive summary statistics from BSL posterior draws
- [`bsl_make_density_summary()`](https://cm401.github.io/ddsynth/reference/bsl_make_density_summary.md)
  : Compute posterior predictive density summary for plotting
- [`bsl_make_posterior_summary()`](https://cm401.github.io/ddsynth/reference/bsl_make_posterior_summary.md)
  : Compute posterior predictive density or CDF summary for plotting
- [`bsl_make_trace_df()`](https://cm401.github.io/ddsynth/reference/bsl_make_trace_df.md)
  : Build a tidy trace data frame from a multi-chain BSL fit list
- [`bsl_run_diagnostics()`](https://cm401.github.io/ddsynth/reference/bsl_run_diagnostics.md)
  : Run automated MCMC diagnostics on BSL fits and save plots/CSV
- [`bsl_summarise_posteriors()`](https://cm401.github.io/ddsynth/reference/bsl_summarise_posteriors.md)
  : Summarise BSL posteriors with point estimates and 95% credible
  intervals
- [`bsl_to_mcmc_list_multi()`](https://cm401.github.io/ddsynth/reference/bsl_to_mcmc_list_multi.md)
  : Convert a multi-chain BSL fit list to coda mcmc.list objects

## Results & Plotting

Functions for summarising and visualising results

## Simulation Metrics

Functions for evaluating simulation study performance

- [`check_coverage()`](https://cm401.github.io/ddsynth/reference/check_coverage.md)
  : Compute coverage for a parameter
- [`compute_iqd()`](https://cm401.github.io/ddsynth/reference/compute_iqd.md)
  : Compute Integrated Quadratic Distance (IQD)
- [`compute_median_bias()`](https://cm401.github.io/ddsynth/reference/compute_median_bias.md)
  : Compute median bias
- [`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md)
  : Compute posterior predictive CDF from a fitted Stan model
- [`extract_quantiles()`](https://cm401.github.io/ddsynth/reference/extract_quantiles.md)
  : Extract quantiles from a predictive CDF summary
- [`summarise_parameters()`](https://cm401.github.io/ddsynth/reference/summarise_parameters.md)
  : Summarise posterior parameter estimates across BSL models
