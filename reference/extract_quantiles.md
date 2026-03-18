# Extract quantiles from a predictive CDF summary

For each requested probability, finds the x value where the CDF (and its
credible bounds) crosses that probability.

## Usage

``` r
extract_quantiles(cdf_summary, probs = c(0.5, 0.95))
```

## Arguments

- cdf_summary:

  A data frame produced by
  [`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md),
  with columns `x`, `median`, `low`, and `high`.

- probs:

  Numeric vector of probabilities to extract (default: `c(0.5, 0.95)`).

## Value

A data frame with columns `quantile`, `quantile_label`, `x_median`,
`x_low`, and `x_high`.
