# Compute LOO-based model weights for all pathogens

For each pathogen, computes Leave-One-Out cross-validation (LOO-CV) for
every valid distribution fit in the chosen `analysis` slot, then
converts the ELPD differences into pseudo-Bayes-factor model weights.

**Why not bridge sampling?**
[`bridgesampling::bridge_sampler()`](https://rdrr.io/pkg/bridgesampling/man/bridge_sampler.html)
requires the compiled C++ Stan model to be present in memory. When a
`stanfit` is saved to RDS and reloaded the internal model object is not
serialised and the call fails. LOO-CV operates entirely on the stored
MCMC draws (`log_lik` parameter) so it works correctly on reloaded fits.

**Method**: [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)
is called on the non-placeholder `log_lik` columns of each fit (the Stan
model stores zero-filled placeholder columns for unused
summary-statistic slots; these are stripped before passing to LOO). The
resulting ELPD estimates are converted to weights via softmax, giving
relative model probabilities under equal model priors — analogous to
Bayes factors computed from marginal likelihoods.

## Usage

``` r
compute_pathogen_model_bayes_factors(all_results, analysis = "filtered")
```

## Arguments

- all_results:

  Nested list produced by `analysis/main_analysis.R`.

- analysis:

  Character scalar. Which analysis slot to use (default `"filtered"`).

## Value

A named list, one entry per pathogen. Each entry is a named numeric
vector of model weights (summing to 1) sorted in decreasing order, keyed
by distribution name. Pathogens with fewer than one valid fit return
`NULL`.

## See also

[`plot_main_figure()`](https://cm401.github.io/ddsynth/reference/plot_main_figure.md)
