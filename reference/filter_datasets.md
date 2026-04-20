# Filter a dataset list by subgroup and/or location

Retains only those entries whose `subgroup` and/or `location` fields
match the requested values. Entries that do not carry the field at all
are **kept** when the corresponding filter argument is `NULL` and
**dropped** when a filter is active (because their group membership is
unknown). Setting both arguments to `NULL` returns the list unchanged.

## Usage

``` r
filter_datasets(datasets, subgroup = NULL, location = NULL)
```

## Arguments

- datasets:

  A named list of datasets in the format accepted by
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).
  Each entry may optionally contain `subgroup` and/or `location`
  character fields.

- subgroup:

  Character vector of subgroup values to retain, or `NULL` (default) to
  skip subgroup filtering.

- location:

  Character vector of location values to retain, or `NULL` (default) to
  skip location filtering.

## Value

A named list containing only the datasets that satisfy all active filter
criteria.

## Examples

``` r
datasets <- list(
  d1 = list(mean = 4.0, sd = 2.4, n = 49,
            subgroup = "tick-bite", location = "Turkey"),
  d2 = list(mean = 6.0, sd = 3.1, n = 12,
            subgroup = "nosocomial", location = "Iran"),
  d3 = list(mean = 5.0, sd = 2.0, n = 30)   # no subgroup/location
)
filter_datasets(datasets, subgroup = "tick-bite")
#> $d1
#> $d1$mean
#> [1] 4
#> 
#> $d1$sd
#> [1] 2.4
#> 
#> $d1$n
#> [1] 49
#> 
#> $d1$subgroup
#> [1] "tick-bite"
#> 
#> $d1$location
#> [1] "Turkey"
#> 
#> 
filter_datasets(datasets, location = c("Turkey", "Iran"))
#> $d1
#> $d1$mean
#> [1] 4
#> 
#> $d1$sd
#> [1] 2.4
#> 
#> $d1$n
#> [1] 49
#> 
#> $d1$subgroup
#> [1] "tick-bite"
#> 
#> $d1$location
#> [1] "Turkey"
#> 
#> 
#> $d2
#> $d2$mean
#> [1] 6
#> 
#> $d2$sd
#> [1] 3.1
#> 
#> $d2$n
#> [1] 12
#> 
#> $d2$subgroup
#> [1] "nosocomial"
#> 
#> $d2$location
#> [1] "Iran"
#> 
#> 
filter_datasets(datasets, subgroup = "nosocomial", location = "Iran")
#> $d2
#> $d2$mean
#> [1] 6
#> 
#> $d2$sd
#> [1] 3.1
#> 
#> $d2$n
#> [1] 12
#> 
#> $d2$subgroup
#> [1] "nosocomial"
#> 
#> $d2$location
#> [1] "Iran"
#> 
#> 
```
