# Check whether Gamma can be reliably fitted from median + IQR summary statistics

The Gamma distribution (dist_type = 2) has a single shape parameter
\\\phi\\ that controls the coefficient of variation CV =
1/\\\sqrt{\phi}\\. When only median + IQR (summary type 2) data are
available, the likelihood gradient with respect to \\\phi\\ comes
exclusively from central order statistics at p = 0.25, 0.50, 0.75. These
central quantiles are most sensitive to location and spread, but carry
weak information about shape compared with tail quantiles (type-1:
min/max) or full frequency data (types 4/5).

As \\\phi\\ grows large the gamma distribution approaches a normal
distribution: the three central quantiles become nearly symmetric around
the mean, and the likelihood surface flattens in the \\\phi\\ direction.
MCMC consequently exhibits very small step sizes, high autocorrelation,
and inference that is dominated by the prior on `log_phi`.

Note that this problem does not apply to type-1 (median + range) or
type-3 (mean + SD) data. Type-1 uses extreme order statistics (min, max)
which lie in the tails where gamma shape sensitivity is largest. Type-3
directly identifies \\\phi\\ via \\\phi = (\text{mean}/\text{SD})^2\\,
giving a sharp, well-defined gradient regardless of how concentrated the
distribution is.

The heuristic uses the moment estimator \$\$\hat{\phi} =
\left(\frac{1.35 \times \text{median}}{\text{IQR}}\right)^2\$\$
(equivalent to CV\\^{-2}\\, the method-of-moments gamma shape estimate)
computed from each type-2 dataset. If the median implied shape across
type-2 datasets exceeds `max_implied_shape`, the function returns
`FALSE`.

## Usage

``` r
gamma_type2_reliable(
  datasets,
  max_implied_shape = 20,
  min_n = 50,
  verbose = TRUE
)
```

## Arguments

- datasets:

  A named list of datasets as accepted by
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).

- max_implied_shape:

  Numeric scalar (default 20). Maximum tolerated median implied shape
  across type-2 datasets. Corresponds to CV \\\approx\\ 0.22 and
  IQR/median \\\approx\\ 0.30. Reduce to be stricter; increase to be
  more permissive.

- min_n:

  Integer scalar (default 50). Minimum acceptable sample size for a
  type-2 dataset. If more than half of the type-2 datasets fall below
  this threshold, an advisory message is printed but `FALSE` is not
  returned — use this as a soft warning only.

- verbose:

  Logical (default `TRUE`). Print a one-line verdict with the reason a
  check failed.

## Value

`TRUE` if type-2 data appear adequate for gamma fitting, `FALSE` if the
implied shape is too large for reliable inference. Returns `TRUE`
silently when no type-2 datasets are present (the heuristic is not
relevant in that case).

## Examples

``` r
if (FALSE) { # \dontrun{
# Concentrated distribution — high implied shape, likely slow
ds_concentrated <- list(
  d1 = list(median = 10, Q1 = 9.2, Q3 = 10.8, n = 50),
  d2 = list(median = 12, Q1 = 11.1, Q3 = 12.9, n = 60)
)
gamma_type2_reliable(ds_concentrated)   # expected: FALSE

# Dispersed distribution — low implied shape, reliable
ds_dispersed <- list(
  d1 = list(median = 10, Q1 = 7, Q3 = 14, n = 80),
  d2 = list(median = 8,  Q1 = 5, Q3 = 12, n = 100)
)
gamma_type2_reliable(ds_dispersed)      # expected: TRUE
} # }
```
