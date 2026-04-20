# Fit Stan model to simulated data

Fit Stan model to simulated data

## Usage

``` r
fit_model(
  sim_data,
  stan_model,
  chains = 4L,
  iter = 10000L,
  warmup = 1000L,
  refresh = 0L,
  control = list(adapt_delta = 0.95, max_treedepth = 12L),
  ...
)
```

## Arguments

- sim_data:

  Simulated data from generate_hierarchical_data

- stan_model:

  Compiled Stan model

- chains:

  Number of Markov chains (default 4)

- iter:

  Total number of iterations per chain (default 10000)

- warmup:

  Number of warmup iterations per chain (default 1000)

- refresh:

  How often to print progress (default 0)

- control:

  List of control parameters passed to Stan (e.g. adapt_delta)

- ...:

  Additional arguments to pass to rstan::sampling()

## Value

Stan fit object
