# Compute posterior predictive CDF from a fitted Stan model

Integrates over posterior draws and between-study random effects to
produce a predictive CDF with pointwise credible bands.

## Usage

``` r
compute_predictive_cdf(
  fit,
  dist_name,
  x_seq = seq(0, 30, length.out = 500),
  n_draws = 500,
  L = 50
)
```

## Arguments

- fit:

  A `stanfit` object returned by
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

- dist_name:

  Character string: `"lognormal"`, `"gamma"`, or `"weibull"`.

- x_seq:

  Numeric vector of evaluation points (default: 500 points on
  `[0, 30]`).

- n_draws:

  Number of posterior draws to use (default: 500).

- L:

  Number of study-level locations to integrate over per draw (default:
  50).

## Value

A data frame with columns `x`, `median`, `mean`, `low`, `high`, and
`model`.
