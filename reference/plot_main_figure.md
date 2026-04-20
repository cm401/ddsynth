# Build the main incubation-period analysis figure

Produces a multi-panel figure with one tile per pathogen. Each tile
shows the posterior predictive CDF (95 % credible ribbon + median line)
for the best-fitting distribution selected by LOO-based model weights
(pseudo-Bayes factors; see
[`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md)).
Dashed segments at the P50 and P95 marks span the 95 % posterior
uncertainty of those quantiles.

When subgroup analyses are present (e.g. country or variant strata),
their median CDFs are overlaid as thin coloured lines without ribbons;
the overall "filtered" analysis remains the primary visual element. The
chosen distribution is labelled in the bottom-right corner together with
the posterior median of the mean delay and shape parameter(s). A
model-weight annotation (\\w\\) is shown in the top-left when one model
clearly dominates (\\w \ge 0.75\\).

## Usage

``` r
plot_main_figure(
  all_results,
  pathogens = NULL,
  pathogen_labels = NULL,
  x_max = NULL,
  n_draws = 500L,
  model_weights = NULL,
  show_subgroups = TRUE,
  max_subgroups = 8L,
  ncol = 4L,
  base_size = 9,
  show_title = TRUE,
  pathogen_highlights = NULL
)
```

## Arguments

- all_results:

  Nested list produced by `analysis/main_analysis.R`.

- pathogen_labels:

  Optional named character vector overriding display labels (e.g.
  `c(COVID_19 = "COVID-19")`).

- x_max:

  Upper limit of the x-axis (days). `NULL` (default) derives each
  panel's limit automatically as the nearest 5-day step at or above the
  99th percentile of the posterior predictive CDF (minimum 10 days).
  Alternatively, supply a single numeric applied to all panels, or a
  named numeric vector with per-pathogen overrides plus an optional
  `"default"` entry.

- n_draws:

  Integer. Posterior draws for
  [`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md)
  for the overall analysis. Subgroup CDFs use `ceiling(n_draws / 2)`.

- model_weights:

  Pre-computed output from
  [`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md).
  If `NULL`, weights are computed internally — this is slow (1–2 hours
  for a full set of pathogens). Pre-compute and cache this object.

- show_subgroups:

  Logical. Overlay subgroup median CDFs (default TRUE).

- max_subgroups:

  Integer. Maximum subgroup lines per panel (default 8).

- ncol:

  Integer. Columns in the assembled figure (default 4).

- base_size:

  Numeric. Base font size for
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  (default 9).

## Value

A patchwork ggplot object.

## See also

[`compute_pathogen_model_bayes_factors()`](https://cm401.github.io/ddsynth/reference/compute_pathogen_model_bayes_factors.md)
