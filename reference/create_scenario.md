# Create a scenario with specific characteristics

Create a scenario with specific characteristics

## Usage

``` r
create_scenario(
  scenario_name,
  dist_type,
  n_datasets,
  n_obs_config = c("fixed", "small_var", "large_var", "custom"),
  n_obs_values = NULL,
  summary_config = c("fixed", "mixed_balanced", "mixed_random",
    "mixed_balanced_with_freq", "mixed_random_with_freq", "custom"),
  summary_values = NULL,
  fixed_summary_type = 1,
  fixed_n_obs = 14,
  mu0 = 2,
  tau = 0.4,
  phi = NULL
)
```

## Arguments

- scenario_name:

  Descriptive name for the scenario

- dist_type:

  Distribution type

- n_datasets:

  Number of datasets

- n_obs_config:

  Configuration for sample sizes: "fixed", "small_var", "large_var",
  "custom"

- n_obs_values:

  Custom vector of sample sizes (if n_obs_config = "custom")

- summary_config:

  Configuration for summary types: "fixed", "mixed_balanced",
  "mixed_random", "custom"

- summary_values:

  Custom vector of summary types (if summary_config = "custom")

- fixed_summary_type:

  Summary type to use when summary_config = "fixed" (default 1).

- fixed_n_obs:

  Sample size to use when n_obs_config = "fixed" (default 14).

- mu0:

  Population mean

- tau:

  Between-study SD

- phi:

  Distribution-specific parameter

## Value

Scenario specification list
