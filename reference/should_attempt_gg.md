# Check whether the Generalised Gamma is likely identifiable from a dataset

The Generalised Gamma (GG, dist_type = 5) has three parameters: location
(\\\mu\\), scale (\\\sigma\\/phi), and shape (\\Q\\/kappa). All datasets
share a single (\\\sigma\\, \\Q\\) pair, so the GG is only identifiable
when datasets consistently imply the same distributional shape. If
different studies show widely different coefficients of variation (CV =
SD/mean), the sampler cannot find a coherent (\\\sigma\\, \\Q\\) and
will exhibit poor mixing or divergences.

Two fast, pre-fit checks are applied:

- CV spread:

  Computes the CV for every dataset using moment approximations (same
  logic as
  [`update_phi_prior()`](https://cm401.github.io/ddsynth/reference/update_phi_prior.md)).
  If `max(CV) / min(CV) > cv_spread_threshold` the CVs are too
  inconsistent to identify the extra GG parameter.

- Information richness:

  The \\Q\\ parameter encodes tail behaviour beyond mean and variance.
  Datasets that supply only summary statistics (mean + SD, median + IQR,
  median + range) provide at most two moments and give weak leverage on
  \\Q\\. If the fraction of datasets with frequency-table or
  interval-censored data (summary types 4 and 5) is below
  `min_rich_fraction`, the shape is too poorly constrained.

## Usage

``` r
should_attempt_gg(
  datasets,
  cv_spread_threshold = 2.5,
  min_rich_fraction = 0.3,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A named list of datasets in the format accepted by
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).

- cv_spread_threshold:

  Numeric scalar (default 2.5). Maximum tolerated ratio of the largest
  to the smallest per-dataset CV. Increase to be more permissive,
  decrease to be stricter.

- min_rich_fraction:

  Numeric scalar in (0, 1\] (default 0.30). Minimum fraction of datasets
  that must be frequency-table or interval-censored (summary types 4 /
  5). Set to 0 to disable this check.

- verbose:

  Logical (default TRUE). Print a one-line verdict with the reason a
  check failed.

## Value

`TRUE` if both checks pass (GG fitting is worth attempting), `FALSE`
otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
should_attempt_gg(datasets_SARS)    # expected: FALSE
should_attempt_gg(datasets_Mpox)
} # }
```
