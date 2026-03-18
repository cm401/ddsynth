# Delay distribution from summary statistics
Methods to compute delay distributions from summary statistics

## Prerequisites

### R
R (≥ 4.0.0) is required. Download it from [CRAN](https://cran.r-project.org/).

### Stan
Several functions in this package use [Stan](https://mc-stan.org/) via the **rstan** package. Stan requires a working C++ toolchain:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and follow the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
- **macOS**: Install the Xcode Command Line Tools (`xcode-select --install`) and follow the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).
- **Linux**: Ensure `g++` and `make` are available (e.g. `sudo apt install build-essential` on Debian/Ubuntu), then install rstan normally.

Install rstan from CRAN once the toolchain is ready:

```r
install.packages("rstan")
```

Verify the installation works by running the built-in example:

```r
example(stan_model, package = "rstan", run.dontrun = TRUE)
```

## Installation

Install the development version of **ddsynth** directly from GitHub using either `remotes` or `pak`:

```r
# Using remotes (install remotes first if needed)
install.packages("remotes")
remotes::install_github("cm401/delay_distribution_data_synthesis")

# Alternatively, using pak
install.packages("pak")
pak::pkg_install("cm401/delay_distribution_data_synthesis")
```

## Quick start

```r
library(ddsynth)
```

See the package vignettes for detailed usage examples:

```r
browseVignettes("ddsynth")
```
