# Run automated MCMC diagnostics on BSL fits and save plots/CSV

Run automated MCMC diagnostics on BSL fits and save plots/CSV

## Usage

``` r
bsl_run_diagnostics(fits_multi, diagnostic_dir = "bsl_diagnostics")
```

## Arguments

- fits_multi:

  Named list of chain lists (one per model).

- diagnostic_dir:

  Directory to save diagnostic outputs (created if absent).

## Value

Invisibly, a data frame of diagnostics across all models.
