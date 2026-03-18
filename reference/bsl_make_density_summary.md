# Compute posterior predictive density summary for plotting

Compute posterior predictive density summary for plotting

## Usage

``` r
bsl_make_density_summary(
  dist_name,
  post,
  x_seq = seq(0, 20, length.out = 300),
  n_draws = 200,
  L = 20
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

## Value

Data frame with columns `x`, `mean`, `low`, `high`, and `model`.
