# Compute the mean Weighted Interval Score for the posterior predictive distribution

Evaluates the posterior predictive distribution against test
observations drawn from the true marginal distribution, using the
Weighted Interval Score (Bracher et al. 2021) — a proper scoring rule
that rewards both calibration and sharpness.

## Usage

``` r
compute_wis(
  fit,
  dist_type,
  n_datasets,
  true_params,
  alpha_levels = c(0.5, 0.8, 0.9, 0.95),
  n_post_draws = 500L,
  n_mc_per_draw = 1000L,
  n_test = 200L,
  pred_q = NULL,
  true_median = NULL
)
```

## Arguments

- fit:

  A `stanfit` object.

- dist_type:

  Distribution type string.

- n_datasets:

  Number of studies in the fitted data.

- true_params:

  Named list with elements `mu0`, `tau`, `phi`, and optionally `kappa`.

- alpha_levels:

  Numeric vector of PI coverage levels (e.g. 0.95 for a 95% PI). Default
  `c(0.50, 0.80, 0.90, 0.95)`.

- n_post_draws:

  Number of posterior draws passed to
  `compute_posterior_predictive_quantile_ci`. Default 500.

- n_mc_per_draw:

  MC samples per posterior draw. Default 1000.

- n_test:

  Number of test observations drawn from the true marginal distribution.
  Default 200.

- pred_q:

  Optional pre-computed vector of posterior predictive median quantiles
  at the probability grid derived from `alpha_levels`, as returned by
  the `$median` element of
  [`compute_posterior_predictive_quantile_ci()`](https://cm401.github.io/ddsynth/reference/compute_posterior_predictive_quantile_ci.md).
  When supplied the internal call to that function is skipped, avoiding
  redundant computation.

- true_median:

  Optional pre-computed true marginal median (scalar), e.g. the `'50%'`
  element from
  [`compute_true_marginal_quantile()`](https://cm401.github.io/ddsynth/reference/compute_true_marginal_quantile.md).
  When supplied the internal call to that function is skipped.

## Value

A list with `wis` (mean WIS, same scale as the delay outcome) and
`rel_wis` (WIS divided by the true marginal median).
