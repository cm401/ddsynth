# Plot true CDFs for simulation study input scenarios

Generates a faceted CDF plot for all distribution families used in the
simulation study, coloured by family (Lancet palette) and differentiated
by linetype for Burr XII and generalised gamma sub-scenarios. CDFs are
evaluated at the population mean (`loc = mu0`); between-study
heterogeneity (`tau = 0.4`) is not shown.

## Usage

``` r
generate_cdf_plot_simulation_study_input()
```

## Value

A `ggplot` object.
