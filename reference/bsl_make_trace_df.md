# Build a tidy trace data frame from a multi-chain BSL fit list

Build a tidy trace data frame from a multi-chain BSL fit list

## Usage

``` r
bsl_make_trace_df(fits_multi)
```

## Arguments

- fits_multi:

  Named list of chain lists (one per model).

## Value

Data frame with columns `theta*`, `iter`, `chain`, and `model`.
