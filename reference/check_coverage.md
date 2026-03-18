# Compute coverage for a parameter

Compute coverage for a parameter

## Usage

``` r
check_coverage(fit, param_name, true_value, level = 0.95)
```

## Arguments

- fit:

  Stan fit object

- param_name:

  Parameter name

- true_value:

  True parameter value

- level:

  Credible interval level (default 0.95)

## Value

Logical indicating whether true value is in credible interval
