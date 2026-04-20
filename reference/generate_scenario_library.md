# Generate a comprehensive set of scenarios

Builds the master scenario grid used by both `run_simulation_study.R`
and the `simulation_study` vignette. All four distribution families that
are identifiable from summary statistics (lognormal, gamma, Weibull,
Burr XII) are supported, plus the generalised gamma for identifiability
research.

## Usage

``` r
generate_scenario_library(
  include_homogeneous = TRUE,
  include_mixed = TRUE,
  include_varied_n = TRUE,
  include_freq_table = FALSE,
  include_burr12 = FALSE,
  include_gengamma = FALSE,
  include_gg_limitation = FALSE
)
```

## Arguments

- include_homogeneous:

  Include scenarios with fixed summary types

- include_mixed:

  Include scenarios with mixed summary types

- include_varied_n:

  Include scenarios with varied sample sizes

- include_freq_table:

  Include frequency-table summary scenarios

- include_burr12:

  Include Burr XII scenarios (`dist_type = "burr12"`)

- include_gengamma:

  Include generalised gamma standard scenarios
  (`dist_type = "gengamma"`)

- include_gg_limitation:

  Include generalised gamma limitation scenarios (only used when
  `include_gengamma = TRUE`)

## Value

Data frame with one row per scenario and columns: `scenario_name`,
`scenario_group`, `scenario_idx`, `dist_type` (one of `"lognormal"`,
`"gamma"`, `"weibull"`, `"burr12"`, `"gengamma"`), `n_datasets`, `mu0`,
`tau`, `phi`, `kappa`, `n_obs_mean`, `n_obs_sd`, `n_obs_min`,
`n_obs_max`, `summary_type_1_prop`–`summary_type_4_prop`,
`summary_type`, `vary_n`.
