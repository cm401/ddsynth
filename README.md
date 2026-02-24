# delaydistributions

**Methods to compute delay distributions from summary statistics**

The **delaydistributions** R package implements several methodologies for
synthesising a parametric delay distribution from the summary statistics
(mean, variance/SD, and optional sample size) of previous observations.

## Methods

| Function | Approach |
|---|---|
| `fit_gamma_mom()` | Method of moments – Gamma family |
| `fit_lognormal_mom()` | Method of moments – Log-Normal family |
| `fit_weibull_mom()` | Method of moments – Weibull family |
| `mle_from_summary()` | Maximum likelihood from summary statistics |
| `bayesian_synthesis()` | Bayesian hierarchical synthesis |

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("cm401/delay_distribution_data_synthesis")
```

## Quick start

```r
library(delaydistributions)

# Method of moments
fit_gamma_mom(mean = 5, variance = 4)
fit_lognormal_mom(mean = 5, sd = 2)

# MLE from multiple study summaries
mle_from_summary(
  mean = c(5.0, 5.2, 4.8),
  sd   = c(2.0, 1.8, 2.1),
  n    = c(100L, 80L, 120L),
  distribution = "gamma"
)

# Bayesian synthesis
bayesian_synthesis(
  mean = c(5.0, 5.2, 4.8),
  sd   = c(2.0, 1.8, 2.1),
  n    = c(100L, 80L, 120L)
)
```

## Package structure

```
R/
├── delaydistributions-package.R   # Package-level documentation
├── method_of_moments.R            # fit_gamma_mom, fit_lognormal_mom, fit_weibull_mom
├── mle_summary_stats.R            # mle_from_summary
├── bayesian_synthesis.R           # bayesian_synthesis
└── utils.R                        # Internal helpers
tests/
└── testthat/                      # Unit tests (testthat 3rd edition)
vignettes/
└── introduction.Rmd               # Getting-started vignette
```

## License

MIT
