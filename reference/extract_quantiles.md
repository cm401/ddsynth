# Extract quantiles from a predictive CDF summary

For each requested probability, finds the x value where the CDF (and its
credible bounds) crosses that probability.

## Usage

``` r
extract_quantiles(cdf_summary, probs = c(0.5, 0.95), cdf_mat = NULL)
```

## Arguments

- cdf_summary:

  A data frame produced by
  [`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md),
  with columns `x`, `median`, `low`, and `high`.

- probs:

  Numeric vector of probabilities to extract (default: `c(0.5, 0.95)`).

- cdf_mat:

  Numeric matrix of posterior CDF draws as returned by
  [`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md)
  (rows = posterior draws, columns = `x_seq` grid points). When
  supplied, the 95% prediction interval bounds (`x_low`, `x_high`) are
  computed by interpolating each draw's CDF to find the x at which it
  crosses `p`, then taking the 2.5% and 97.5% quantiles across draws.
  This is fully consistent with the ribbon in the CDF plot (both derive
  from the same `cdf_mat`). If `NULL`, falls back to inverting the
  summary credible bands, which can fail near the tails.

## Value

A data frame with columns `quantile`, `quantile_label`, `x_low`, and
`x_high`. `x_low` and `x_high` are the 2.5% and 97.5% bounds of the 95%
prediction interval for that quantile, on the same scale as the `x`
column of `cdf_summary`.
