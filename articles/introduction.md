# Introduction to ddsynth

``` r
library(ddsynth)
```

## Overview

`ddsynth` provides methods to compute delay distributions from summary
statistics (e.g. a reported mean and standard deviation).

This vignette is a **placeholder** that will be expanded once the main
analysis functions have been added to the package.

## Planned workflow

1.  Supply summary statistics describing the delay (mean, SD, or
    quantiles).
2.  Choose a parametric family (e.g. log-normal, gamma, Weibull).
3.  Fit the distribution — either analytically or via a Stan model
    stored in `inst/stan/`.
4.  Inspect and use the fitted distribution.

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ddsynth_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          codetools_0.2-20 
#>  [5] fastmap_1.2.0     xfun_0.57         iterators_1.0.14  cachem_1.1.0     
#>  [9] knitr_1.51        htmltools_0.5.9   rmarkdown_2.31    lifecycle_1.0.5  
#> [13] cli_3.6.5         foreach_1.5.2     sass_0.4.10       pkgdown_2.2.0    
#> [17] textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2 compiler_4.5.3   
#> [21] tools_4.5.3       ragg_1.5.2        evaluate_1.0.5    bslib_0.10.0     
#> [25] yaml_2.3.12       jsonlite_2.0.0    rlang_1.1.7       fs_2.0.1
```
