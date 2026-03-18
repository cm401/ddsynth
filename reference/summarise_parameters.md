# Summarise posterior parameter estimates across BSL models

Summarise posterior parameter estimates across BSL models

## Usage

``` r
summarise_parameters(posterior_samples_list)
```

## Arguments

- posterior_samples_list:

  Named list of posterior sample matrices (columns: mu0, log_tau,
  log_phi).

## Value

Data frame with mean, median, and 95% CrI per parameter per model.
