# Compute credible interval for posterior predictive quantiles

For each requested probability `p`, draws from the posterior predictive
distribution and returns the 2.5th, 50th, and 97.5th percentiles of the
resulting quantile distribution across posterior draws.

## Usage

``` r
compute_posterior_predictive_quantile_ci(
  fit,
  dist_type,
  n_datasets,
  probs = c(0.5, 0.95),
  n_post_draws = 500L,
  n_mc_per_draw = 1000L
)
```

## Arguments

- fit:

  A `stanfit` object.

- dist_type:

  Distribution type string.

- n_datasets:

  Number of studies in the fitted data (used to decide whether to
  integrate over mu0/tau or use mean(loc_d)).

- probs:

  Probabilities for which predictive quantiles are computed. Default
  `c(0.5, 0.95)`.

- n_post_draws:

  Number of posterior draws to use. Default 500.

- n_mc_per_draw:

  Number of predictive samples per posterior draw. Default 1000.

## Value

A list with elements `lower`, `median`, and `upper`: the 2.5th, 50th,
and 97.5th percentiles of the posterior distribution of each quantile,
each a numeric vector of length `length(probs)`.
