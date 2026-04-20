# Extract a tidy summary of all built-in datasets

Loops through all per-pathogen dataset lists and returns one row per
dataset entry. Sample size `n` is taken from the `n` field if present,
otherwise computed as `sum(freq_count)`.

## Usage

``` r
extract_dataset_summary()
```

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns `pathogen`, `pathogen_group`, `dataset_id`, `country`, `n`,
`n_log` (log1p of n), `subgroup`, and `source`.
