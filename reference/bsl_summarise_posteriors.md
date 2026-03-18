# Summarise BSL posteriors with point estimates and 95% credible intervals

Summarise BSL posteriors with point estimates and 95% credible intervals

## Usage

``` r
bsl_summarise_posteriors(posterior_samples_list, L = 50)
```

## Arguments

- posterior_samples_list:

  Named list of posterior sample matrices.

- L:

  Number of study-level locations to integrate over (default 50).

## Value

Data frame with posterior summaries per model.
