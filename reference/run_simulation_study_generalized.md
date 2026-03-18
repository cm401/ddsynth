# Run simulation study with generalized scenarios

Run simulation study with generalized scenarios

## Usage

``` r
run_simulation_study_generalized(
  n_sim,
  scenarios_df,
  stan_model,
  seed = 123,
  save_name = "simulation_results_tmp_general.rds",
  n_cores = 1,
  chains = 4
)
```

## Arguments

- n_sim:

  Number of simulation replicates per scenario

- scenarios_df:

  Data frame from generate_scenario_library()

- stan_model:

  Compiled Stan model

- seed:

  Random seed

- save_name:

  File path to save final results as RDS.

- n_cores:

  Number of parallel workers requested. Automatically capped at
  `floor(detectCores() / chains)` to avoid CPU oversubscription.

- chains:

  Number of Stan chains per fit (default 4). Used to compute the worker
  cap.

## Value

Data frame with one row per simulation replicate and columns for
scenario metadata, true parameter values, coverage, bias, and IQD
metrics.
