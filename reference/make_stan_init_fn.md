# Create a Stan initialisation function from prior means

Returns a zero-argument function that initialises each chain at the
prior means stored in `stan_data`. Passing this to
[`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
via `init = make_stan_init_fn(stan_data)` avoids the default
`Uniform(-2, 2)` draws on the unconstrained scale, which send the GG and
Burr XII log-posteriors to \\-\infty\\ during initialisation.

## Usage

``` r
make_stan_init_fn(stan_data)
```

## Arguments

- stan_data:

  Stan data list as returned by
  [`prepare_stan_data_from_datasets`](https://cm401.github.io/ddsynth/reference/prepare_stan_data_from_datasets.md).

## Value

A zero-argument function suitable for `rstan::sampling(init = ...)`.
