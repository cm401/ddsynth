# Compute true marginal quantiles of the predictive distribution

Uses Monte Carlo integration over the study-level location hierarchy to
compute the true population-level quantiles of the predictive
distribution.

## Usage

``` r
compute_true_marginal_quantile(
  dist_type,
  mu0,
  tau,
  phi,
  kappa = 1,
  probs = c(0.5, 0.95),
  n_mc = 5000L
)
```

## Arguments

- dist_type:

  Distribution type: one of `"lognormal"`, `"gamma"`, `"weibull"`,
  `"burr12"`, `"gengamma"`.

- mu0:

  True population mean (location hyperparameter).

- tau:

  True between-study standard deviation.

- phi:

  True distribution-specific shape/scale parameter.

- kappa:

  True third distribution parameter (Burr XII k or GG Q). Ignored for
  2-parameter distributions.

- probs:

  Numeric vector of probabilities for which quantiles are computed.
  Default `c(0.5, 0.95)`.

- n_mc:

  Number of Monte Carlo draws. Default 5000.

## Value

A named numeric vector of quantiles (names are e.g. `"50%"`, `"95%"`).
