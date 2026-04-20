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
  custom_priors = list()
)
```

## Arguments

- datasets:

  A named list of lists. Each element must contain one of the following
  combinations of summary statistics:

  `median`, `min`, `max`

  :   Median and range (summary type 1). `n` (sample size) is required.

  `median`, `Q1`, `Q3`

  :   Median and inter-quartile range (summary type 2). `n` is required.

  `mean`, `sd`

  :   Mean and standard deviation (summary type 3). `n` is required.

  `freq_value`, `freq_count`

  :   Frequency table of (value, count) pairs (summary type 4). `n` is
      optional and defaults to `sum(freq_count)`.

  `freq_lower`, `freq_upper`, `freq_count`

  :   Interval-censored frequency table (summary type 5). Each entry
      gives the lower and upper bound of the censoring interval and the
      count of individuals in that interval. When
      `freq_lower[i] == freq_upper[i]` the observation is treated as
      exact. `n` is optional and defaults to `sum(freq_count)`.

  Each element may also contain an optional `source` field — a free-text
  character string recording the bibliographic reference for that
  dataset (e.g. `"Surname (year), doi: doi.org/xyz"`). This field is
  ignored during Stan data preparation and is never passed to the model.

- dist_type:

  Integer distribution code: `1` = log-normal, `2` = gamma, `3` =
  Weibull. Defaults to `1`.

- use_custom_priors:

  Integer flag (0 or 1) for custom prior use. Currently unused; reserved
  for future extension. Defaults to `0`.

- custom_priors:

  Named list of prior overrides. Any values not supplied fall back to
  distribution-appropriate defaults (see Details). Recognised names:
  `mu0_sd`, `log_tau_mean`, `log_tau_sd`, `log_phi_mean`, `log_phi_sd`.

## Value

A named list suitable for passing to
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
as the `data` argument. The list always includes `freq_lower` and
`freq_upper` fields (populated with zeros for non-type-5 datasets), as
these are required by the Stan model regardless of which summary types
are present.

## Details

**Distribution-specific defaults for `log_phi`:**

Because `phi` has a different meaning in each distribution, the default
prior for `log_phi_mean` is chosen per `dist_type`:

|             |              |            |                        |                  |
|-------------|--------------|------------|------------------------|------------------|
| `dist_type` | Distribution | `phi`      | Default `log_phi_mean` | Prior median phi |
| 1           | Lognormal    | log-SD (σ) | -0.7                   | 0.50             |
| 2           | Gamma        | shape      | 2.5                    | 12.2             |
| 3           | Weibull      | shape      | 1.0                    | 2.7              |

Users can override any individual prior by passing only the relevant
element(s) in `custom_priors`, e.g.
`custom_priors = list(log_phi_mean = 3.0)` — all other priors will use
the distribution-appropriate defaults above.

## Note

**Backward compatibility:** The Stan model requires `freq_lower` and
`freq_upper` to be present in the data list for all runs, including
those that contain only type 1–4 datasets. This is handled automatically
when using this function. If you construct the Stan data list manually
(rather than via this function), you must include these fields
explicitly, e.g.:

    stan_data$freq_lower <- rep(0, stan_data$n_freq_total)
    stan_data$freq_upper <- rep(0, stan_data$n_freq_total)
