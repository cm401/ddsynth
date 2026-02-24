# Stan models

Place your Stan model files (`.stan`) in this directory.

## Convention

Each Stan model should be a self-contained `.stan` file, for example:

```
inst/stan/
├── lognormal_fit.stan   # fit a log-normal delay distribution
├── gamma_fit.stan       # fit a gamma delay distribution
└── weibull_fit.stan     # fit a Weibull delay distribution
```

## Accessing models at runtime

When the package is installed, everything under `inst/` is copied to the
package root, so a model can be located from R with:

```r
system.file("stan", "lognormal_fit.stan", package = "delaydistribution")
```

This path can then be passed to `rstan::stan()`, `cmdstanr::cmdstan_model()`,
or similar interfaces.
