# Prepare Stan data from a list of dataset summaries

Converts a list of dataset descriptors (each providing summary
statistics and a sample size) into the named list expected by the
`hierarchical_data_synthesis_summary_stats` Stan model.

## Usage

``` r
prepare_stan_data_from_datasets(
  datasets,
  dist_type = 1,
  use_custom_priors = 0,
  custom_priors = NULL
)
```

## Arguments

- datasets:

  A named list of lists. Each element must contain `n` (sample size) and
  one of the following combinations of summary statistics:

  `median`, `min`, `max`

  :   Median and range (summary type 1).

  `median`, `Q1`, `Q3`

  :   Median and inter-quartile range (summary type 2).

  `mean`, `sd`

  :   Mean and standard deviation (summary type 3).

  `freq_value`, `freq_count`

  :   Frequency table of (value, count) pairs (summary type 4). `n` is
      optional and defaults to `sum(freq_count)`.

- dist_type:

  Integer distribution code: `1` = log-normal, `2` = gamma, `3` =
  Weibull. Defaults to `1`.

- use_custom_priors:

  Integer flag (0 or 1) for custom prior use. Currently unused; reserved
  for future extension. Defaults to `0`.

- custom_priors:

  Optional list of custom prior values. Currently unused.

## Value

A named list suitable for passing to
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
as the `data` argument.
