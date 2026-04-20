# Create a global map of dataset geographic distribution

Plots each dataset entry as a dot on a world map. Dot size is
proportional to `log(n+1)` and dot colour indicates pathogen. Datasets
with `country = "Mixed"` are excluded from the map and listed in an
inset box.

## Usage

``` r
create_data_map(
  data = NULL,
  exclude_no_n = FALSE,
  default_n_log = 1.5,
  alpha = 0.6,
  jitter_amount = 3,
  title = "Geographic distribution of delay distribution datasets"
)
```

## Arguments

- data:

  A data frame from
  [`extract_dataset_summary()`](https://cm401.github.io/ddsynth/reference/extract_dataset_summary.md).
  If `NULL`,
  [`extract_dataset_summary()`](https://cm401.github.io/ddsynth/reference/extract_dataset_summary.md)
  is called internally.

- exclude_no_n:

  Logical. If `TRUE`, entries with missing sample size are dropped. If
  `FALSE` (default), they appear at size `default_n_log`.

- default_n_log:

  Numeric. Log-scale dot size used when `n` is missing.

- alpha:

  Numeric (0–1). Point transparency.

- jitter_amount:

  Numeric. Maximum spatial jitter in degrees lat/lon, applied to
  separate overlapping points from the same country.

- title:

  Character. Plot title.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

For reproducible jitter, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before this function.
