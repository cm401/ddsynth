# Convert a multi-chain BSL fit list to coda mcmc.list objects

Convert a multi-chain BSL fit list to coda mcmc.list objects

## Usage

``` r
bsl_to_mcmc_list_multi(fits_multi)
```

## Arguments

- fits_multi:

  Named list of chain lists (one per model).

## Value

Named list of
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
objects.
