# Generate a LaTeX summary table split into three panels by data availability

Produces a single booktabs-style LaTeX table with the same columns as
[`generate_results_table()`](https://cm401.github.io/ddsynth/reference/generate_results_table.md),
but with pathogens divided into three labelled panels (A, B, C)
separated by horizontal rules and bold panel-label header rows.

## Usage

``` r
generate_results_table_split(
  all_results,
  group_a = c("SARS", "COVID_19", "Measles", "Mpox"),
  group_b = c("Nipah", "EVD", "MERS", "Cholera", "CCHF", "Dengue", "YFV"),
  group_c = c("MVD", "Lassa", "Zika", "RVF"),
  panel_labels = c("A", "B", "C"),
  model_weights = NULL,
  show_subgroups = TRUE,
  rhat_threshold = 1.05,
  ci_prob = 0.95,
  pathogen_labels = NULL,
  caption = paste0("Incubation period estimates by pathogen. ",
    "Values are posterior medians with 95\\% credible intervals (CI). ",
    "Pathogens are grouped by data availability: panel A (substantial data), ",
    "panel B (moderate data), and panel C (limited data). ",
    "Subgroup analyses are shown in italics beneath each pathogen."),
  label = "tab:incubation_periods"
)
```

## Arguments

- all_results:

  Nested list produced by `analysis/main_analysis.R`.

- group_a:

  Character vector of pathogen keys for panel A (high data).

- group_b:

  Character vector of pathogen keys for panel B (moderate data).

- group_c:

  Character vector of pathogen keys for panel C (low data).

- panel_labels:

  Length-3 character vector of panel labels (default
  `c("A", "B", "C")`).

- model_weights:

  Pre-computed output from
  [`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md).
  Computed once internally if `NULL`.

- show_subgroups:

  Logical. Include subgroup sub-rows (default `TRUE`).

- rhat_threshold:

  Numeric. Maximum acceptable Rhat (default 1.05).

- ci_prob:

  Numeric. Width of the credible interval (default 0.95).

- pathogen_labels:

  Optional named character vector overriding display names.

- caption:

  LaTeX `\caption{}` string.

- label:

  LaTeX `\label{}` string.

## Value

A single character string containing the complete LaTeX table
environment.

## See also

[`generate_results_table()`](https://cm401.github.io/ddsynth/reference/generate_results_table.md),
[`plot_main_figure_split()`](https://cm401.github.io/ddsynth/reference/plot_main_figure_split.md)
