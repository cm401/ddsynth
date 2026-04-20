# Split main figure into three panels by data availability

Calls
[`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md)
three times — once per pathogen group — and assembles the results into a
single vertically-stacked figure labelled A, B, and C.

## Usage

``` r
plot_main_figure_split(
  all_results,
  group_a = c("SARS", "COVID_19", "Flu", "Measles"),
  group_b = c("Nipah", "EVD", "MERS", "CCHF", "Cholera", "Typhoid", "Dengue", "YFV",
    "Mpox", "Smallpox"),
  group_c = c("MVD", "Lassa", "Zika", "RVF"),
  panel_labels = c("A", "B", "C"),
  model_weights = NULL,
  ncol = 4L,
  ...
)
```

## Arguments

- all_results:

  Nested results list (same as
  [`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md)).

- group_a:

  Character vector of pathogen keys for panel A (high data).

- group_b:

  Character vector of pathogen keys for panel B (moderate data).

- group_c:

  Character vector of pathogen keys for panel C (low data).

- panel_labels:

  Length-3 character vector of panel labels.

- model_weights:

  Pre-computed model weights from
  [`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md).
  If `NULL`, computed once and shared across all three panels.

- ncol:

  Number of columns within each panel group.

- ...:

  Additional arguments forwarded to
  [`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md).

## Value

A patchwork figure.
