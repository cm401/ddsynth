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
  L = 2000
)
```

## Arguments

- fit:

  A `stanfit` object returned by
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

- dist_name:

  Character string: `"lognormal"`, `"gamma"`, `"weibull"`, `"burr12"`
  (Burr Type XII), or `"gengamma"` (Generalised Gamma, Prentice
  parameterisation).

- x_seq:

  Numeric vector of evaluation points (default: 500 points on
  `[0, 30]`).

- n_draws:

  Number of posterior draws to use (default: 500).

- L:

  Number of study-level locations to integrate over per draw (default:
  2000). When `n_datasets < 5`, `mean(loc_d)` is used directly for all
  `L` locations (i.e. no between-study sampling) for consistency with
  the Stan generated quantities block — see
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md)
  for details.

## Value

A data frame with columns `x`, `median`, `mean`, `low`, `high`, and
`model`.
