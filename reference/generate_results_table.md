# Generate a LaTeX summary table of incubation period estimates

Produces a booktabs-style LaTeX table with one row per pathogen and
(optionally) indented italic sub-rows for each subgroup analysis.

## Usage

``` r
generate_results_table(
  all_results,
  pathogens = NULL,
  model_weights = NULL,
  show_subgroups = TRUE,
  rhat_threshold = 1.05,
  ci_prob = 0.95,
  pathogen_labels = NULL,
  caption =
    paste0("Incubation period estimates by pathogen, grouped by epidemiological category. ",
    "Values are posterior medians with 95\\% credible intervals (CI). ",
    "$\\tau$ is the between-study heterogeneity SD. ",
    "Subgroup analyses are shown in italics beneath each pathogen."),
  label = "tab:incubation_periods"
)
```

## Arguments

- all_results:

  Nested list produced by `analysis/main_analysis.R`.

- model_weights:

  Pre-computed output from
  [`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md).
  Computed internally if `NULL` (slow — pre-compute and cache).

- show_subgroups:

  Logical. Include subgroup sub-rows (default `TRUE`).

- rhat_threshold:

  Numeric. Maximum acceptable Rhat for convergence; fits with max Rhat
  above this are skipped (default 1.05).

- ci_prob:

  Numeric. Width of the credible interval, e.g. 0.95 gives a 2.5%–97.5%
  interval (default 0.95).

- pathogen_labels:

  Optional named character vector overriding the display names used in
  the table rows.

- caption:

  LaTeX `\caption{}` string.

- label:

  LaTeX `\label{}` string.

## Value

A single character string containing the complete LaTeX table
environment, ready to paste into an Overleaf document.

## Details

For every row the table reports:

- the best-fitting distribution (selected with the same
  convergence-filtered LOO-weight logic as
  [`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md));

- the posterior predictive **median** with 95% credible interval;

- the posterior predictive **95th percentile** (\\p\_{95}\\) with 95%
  CI;

- the population-level location parameter \\\mu_0\\ with 95% CI;

- the shape parameter \\\phi\\ with 95% CI;

- the second shape parameter \\\kappa\\ with 95% CI (Burr XII and
  Generalised Gamma only; shown as \\-\\ for other families).

All intervals are formatted as `median (lower, upper)`.

The table requires the LaTeX packages booktabs and multirow in the
preamble.

## See also

[`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md),
[`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md)
