# Run pre-inference checks on a list of datasets

Performs a suite of fast, pre-MCMC checks to detect data issues that are
likely to cause convergence problems. Checks are run in order of
increasing computational cost and a summary is printed to the console.

## Usage

``` r
pre_inference_checks(
  datasets,
  stan_model,
  dist_type = 1,
  custom_priors = list(),
  phi_outlier_threshold = 5,
  phi_grid = seq(0.5, 50, by = 0.5),
  n_sim = 2000,
  loo_iter = 4000,
  loo_chains = 2,
  verbose = TRUE,
  filter = FALSE
)
```

## Arguments

- datasets:

  A named list of datasets in the format accepted by
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).

- stan_model:

  A compiled Stan model object from
  [`rstan::stan_model()`](https://mc-stan.org/rstan/reference/stan_model.html).

- dist_type:

  Integer distribution code: `1` = log-normal, `2` = gamma, `3` =
  Weibull. Defaults to `1`.

- custom_priors:

  Optional named list of prior overrides passed to
  [`prepare_stan_data_from_datasets()`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).

- phi_outlier_threshold:

  Multiplier used in the method-of-moments check. A dataset is flagged
  if its implied `phi` exceeds
  `phi_outlier_threshold * median(implied_phi)`. Defaults to `5`.

- phi_grid:

  Numeric vector of `phi` values for the log-likelihood surface scan.
  Defaults to `seq(0.5, 50, by = 0.5)`.

- n_sim:

  Number of draws for the prior predictive check. Defaults to `2000`.

- loo_iter:

  Number of MCMC iterations per chain for the leave-one-out
  single-dataset fits. Defaults to `4000`.

- loo_chains:

  Number of chains for the leave-one-out fits. Defaults to `2`.

- verbose:

  Logical. If `TRUE` (default), prints a formatted summary of all check
  results to the console.

- filter:

  Logical. If `TRUE`, any dataset flagged by at least one per-dataset
  check (method-of-moments outlier, outside prior predictive interval,
  or non-overlapping LOO phi posterior) is removed from the returned
  dataset list. Defaults to `FALSE`.

## Value

A named list with elements:

- `mom_consistency`:

  Data frame of implied `phi` per dataset with an `is_outlier` flag.

- `prior_predictive`:

  Data frame of 95% prior predictive intervals for the implied SD of
  each dataset, with an `outside_prior_pi` flag.

- `map_probe`:

  List with `phi_map` (MAP estimate of `phi`) and `map_converged`
  logical.

- `ll_surface`:

  Data frame of `phi` vs `log_prob` from the surface scan.

- `loo_fits`:

  Data frame of per-dataset `phi` posterior summaries from the
  leave-one-out fits.

- `datasets`:

  The input `datasets` list, filtered to remove flagged datasets when
  `filter = TRUE`, otherwise identical to the input.

## Details

The five checks performed are:

- 1\. Method-of-moments consistency:

  Estimates `phi` from each dataset individually using moment-based
  approximations and flags any dataset whose implied `phi` is more than
  `phi_outlier_threshold` times the median of all implied values.

- 2\. Prior predictive compatibility:

  Simulates summary statistics from the prior and checks whether each
  observed value falls within the 95% prior predictive interval.
  Datasets outside this range suggest a prior–data mismatch.

- 3\. MAP optimisation probe:

  Runs
  [`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)
  as a fast proxy for MCMC convergence. Failure or extreme `phi` at the
  MAP estimate is a reliable early warning that HMC will struggle.

- 4\. Log-likelihood surface scan:

  Evaluates the joint log-posterior over a grid of `phi` values (other
  parameters held at the MAP). A multimodal or sharply peaked surface
  explains treedepth exhaustion.

- 5\. Leave-one-out single-dataset fits:

  Fits the model to each dataset individually and compares the resulting
  `phi` posteriors. Non-overlapping credible intervals identify the
  specific datasets driving tension.
