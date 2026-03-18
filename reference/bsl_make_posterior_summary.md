# Compute posterior predictive density or CDF summary for plotting

Compute posterior predictive density or CDF summary for plotting

## Usage

``` r
bsl_make_posterior_summary(
  dist_name,
  post,
  x_seq = seq(0, 20, length.out = 400),
  n_draws = 200,
  L = 20,
  posterior_cdf = FALSE
)
```

## Arguments

- dist_name:

  Distribution name: `"lognormal"`, `"gamma"`, or `"weibull"`.

- post:

  Matrix or data frame of posterior draws.

- x_seq:

  Numeric vector of evaluation points.

- n_draws:

  Number of posterior draws to use.

- L:

  Number of study-level locations to integrate over per draw.

- posterior_cdf:

  Logical; if `TRUE`, returns cumulative density instead.

## Value

Data frame with columns `x`, `mean`, `low`, `high`, and `model`.
