# Matched-moments distribution comparison plot

Solves for the parameters of each distribution family such that all
curves share the same mean and coefficient of variation (CV) as the
reference Gamma scenario (mu0\\=\\2, phi\\=\\3.5). The resulting PDFs
are plotted on a faceted panel so that differences in distributional
*shape* — tail heaviness, skewness, modality — are visible even when the
first two moments are identical.

## Usage

``` r
generate_matched_moments_plot(burr_kappas = 2, gg_kappas = c(0.5, 1.5))
```

## Arguments

- burr_kappas:

  Numeric vector of kappa (k) values for Burr XII curves. Default:
  `c(2, 3, 5)`.

- gg_kappas:

  Numeric vector of Q values for the generalised gamma curves. Default:
  `c(0.1, 0.5, 1.0, 2.0, 3.0)`.

## Value

A `ggplot` object.
