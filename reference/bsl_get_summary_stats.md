# Compute predictive summary statistics from BSL posterior draws

Compute predictive summary statistics from BSL posterior draws

## Usage

``` r
bsl_get_summary_stats(dist_name, post, L = 50)
```

## Arguments

- dist_name:

  Distribution name: `"lognormal"`, `"gamma"`, or `"weibull"`.

- post:

  Matrix or data frame of posterior draws (columns: mu0, log_tau,
  log_phi).

- L:

  Number of study-level locations to integrate over per draw (default
  50).

## Value

Data frame with columns `mean`, `q90`, and `median` per posterior draw.
