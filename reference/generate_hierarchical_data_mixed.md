# Generate data from hierarchical model with mixed summary types and sample sizes

Generate data from hierarchical model with mixed summary types and
sample sizes

## Usage

``` r
generate_hierarchical_data_mixed(
  n_datasets,
  n_obs,
  dist_type = c("lognormal", "gamma", "weibull", "burr12", "gengamma"),
  mu0,
  tau,
  phi,
  kappa = 1,
  summary_type = NULL
)
```

## Arguments

- n_datasets:

  Number of datasets to generate

- n_obs:

  Vector of sample sizes for each dataset (can vary)

- dist_type:

  Distribution type: "lognormal", "gamma", or "weibull"

- mu0:

  Population mean (location parameter)

- tau:

  Between-study standard deviation

- phi:

  Distribution-specific shape/scale parameter

- summary_type:

  Summary type specification. Can be: `NULL` (random mix of types 1-3),
  a single integer 1-4 (all datasets use that type), a vector of length
  `n_datasets` (one type per dataset), a length-3 probability vector
  (sample from types 1-3 with those probabilities), or a length-4
  probability vector (sample from types 1-4 with those probabilities).
  Type 4 produces a frequency table of observations rounded to the
  nearest whole day.

## Value

List containing true parameters and observed summary statistics
