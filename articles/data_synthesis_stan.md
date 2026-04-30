# Data Synthesis Using Summary Statistics with Stan

## Overview

This vignette demonstrates how to synthesise delay distribution
estimates from multiple studies, where each study reports only **summary
statistics** (median, range, IQR, or mean ± SD) rather than individual
observations.

The approach uses a **hierarchical Bayesian model** implemented in Stan.
Study-level location parameters are drawn from a shared population
distribution, allowing information to be pooled across studies while
still accounting for between-study heterogeneity.

Two parametric families are compared — **log-normal** and **Weibull** —
using both bridge sampling (log marginal likelihood) and leave-one-out
cross-validation (LOO-CV).

## Setup

``` r

library(ddsynth)
library(rstan)
library(bridgesampling)
library(loo)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(patchwork)

# Use all available cores and cache compiled Stan models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
```

## Compile the Stan Model

The Stan model file is shipped with the package under `inst/stan/`.

``` r

stan_model_path <- system.file(
  "stan", "hierarchical_data_synthesis_summary_stats.stan",
  package = "ddsynth"
)
stan_model_code <- stan_model(stan_model_path)
```

## Prepare Data

Each element of `datasets` describes one study via its reported summary
statistics and sample size. Three summary-statistic types are supported:

| `summary_type` | Fields required                                   |
|----------------|---------------------------------------------------|
| 1              | `median`, `min`, `max`                            |
| 2              | `median`, `Q1`, `Q3`                              |
| 3              | `mean`, `sd`                                      |
| 4              | `freq_value`, `freq_count` (for frequency tables) |

``` r

datasets <- list(
  d1  = list(median = 10.0, min =  9.0, max = 12.0, n =  4),
  d2  = list(median =  4.0, min =  2.0, max =  7.0, n =  6),
  d3  = list(median =  9.0, min =  6.0, max = 14.0, n = 11),
  d4  = list(median =  9.5, min =  4.0, max = 14.0, n = 22),
  d5  = list(median =  8.0, min =  3.0, max = 20.0, n = 15),
  d6  = list(median =  9.0, min =  6.0, max = 11.0, n = 14),
  d7  = list(median = 10.0, min =  8.0, max = 15.0, n = 11),
  d8  = list(median =  9.0, min =  6.0, max = 11.0, n = 11),
  d9  = list(median =  9.0, Q1  =  8.0, Q3  = 11.0, n = 82),
  d10 = list(mean   =  9.3, sd  =  1.9,             n = 18),
  d11 = list(
    freq_value = c(6, 7, 8, 9, 11, 12, 14),
    freq_count = c(1, 1, 3, 2,  1,  2,  1)
  )
)

custom_priors = list(mu0_sd = 1.0,
                     log_tau_mean = 0.2,
                     log_tau_sd = 0.5,
                     log_phi_mean = -0.7,
                     log_phi_sd   = 0.5)

checks   <- pre_inference_checks(datasets_Cholera, stan_model_code, dist_type = 2, filter = TRUE)
datasets_clean <- checks$datasets 
datasets_clean <- datasets_clean[!sapply(datasets_clean, is.null)]
stan_data <- prepare_stan_data_from_datasets(datasets_clean,custom_priors = custom_priors)
```

## Fit Models

We fit the hierarchical model under each candidate distribution family.

``` r

distributions <- c("lognormal", "gamma", "weibull", "burr", "gg" )
#distributions <- c("lognormal")
dist_codes    <- c(lognormal = 1, gamma = 2, weibull = 3, burr = 4, gg = 5)

fits           <- list()
bridge_samples <- list()

for (dist in distributions) {
  cat("\n\n========== Fitting", dist, "model ==========\n")

  stan_data$dist_type     <- dist_codes[dist]
  #stan_data$log_phi_mean  <- if (dist == "gamma") log(10) else log(1.5)
  stan_data               <- update_phi_prior(stan_data, datasets)  
  
  # Initialise from prior means to avoid log(0) rejections for GG / Burr XII
  init_fn <- make_stan_init_fn(stan_data)

  # GG has a strongly correlated (phi, kappa) posterior: many combinations of
  # sigma and Q fit typical summary statistics equally well, creating an
  # elongated ridge the diagonal mass matrix (diag_e) cannot traverse
  # efficiently.  Switching to a dense mass matrix (dense_e) lets NUTS learn
  # the phi-kappa correlation during warmup and precondition its leapfrog steps
  # to align with the ridge, resolving the max_treedepth saturation seen with
  # diag_e even at max_treedepth = 12.  Warmup is increased to 4000 to give
  # the mass-matrix adaptation enough iterations to converge.
  ctrl <- if (dist == "gg") {
    list(adapt_delta = 0.95, max_treedepth = 12, metric = "dense_e")
  } else {
    list(adapt_delta = 0.9,  max_treedepth = 10)
  }

  fit <- sampling(
    stan_model_code,
    data    = stan_data,
    chains  = 4,
    iter    = if (dist == "gg") 14000L else 12000L,
    warmup  = if (dist == "gg")  4000L else  2000L,
    thin    = 1,
    init    = init_fn,
    control = ctrl,
    seed    = 123
  )

  fits[[dist]]           <- fit
  bridge_samples[[dist]] <- bridge_sampler(fit, silent = TRUE)

  print(fit, pars = c("mu0", "tau", "phi",
                      "pred_mean", "pred_median", "pred_q90", "pred_q95"))
}
```

## Model Comparison

### Bridge Sampling (Log Marginal Likelihood)

``` r

log_ml    <- sapply(bridge_samples, function(x) x$logml)
rel_prob  <- exp(log_ml - max(log_ml))
rel_prob  <- rel_prob / sum(rel_prob)

model_comparison <- data.frame(
  model                   = names(log_ml),
  log_marginal_likelihood = log_ml,
  relative_probability    = rel_prob
)

cat("\n========== Model Comparison ==========\n")
print(model_comparison)

# Bayes factors relative to the best model
best_model <- names(which.max(log_ml))
bf_table   <- data.frame(
  model      = names(log_ml),
  BF_vs_best = exp(log_ml - max(log_ml))
)
cat("\nBayes Factors (relative to best model:", best_model, ")\n")
print(bf_table)
```

### LOO-CV Comparison

``` r

loo_list <- lapply(fits, loo, save_psis = TRUE)
loo_comparison <- loo_compare(loo_list)
cat("\n========== LOO Model Comparison ==========\n")
print(loo_comparison)
```

## Convergence Diagnostics

``` r

diagnostic_summary <- list()

for (dist in distributions) {
  cat("\n========== Diagnostics for", dist, "==========\n")

  summ      <- summary(fits[[dist]])$summary
  rhat_vals <- summ[, "Rhat"]
  n_eff     <- summ[, "n_eff"]

  cat("Max Rhat:", max(rhat_vals, na.rm = TRUE), "\n")
  cat("Parameters with Rhat > 1.05:", sum(rhat_vals > 1.05, na.rm = TRUE), "\n")
  cat("Min n_eff:", min(n_eff, na.rm = TRUE), "\n")
  cat("Parameters with n_eff < 400:", sum(n_eff < 400, na.rm = TRUE), "\n")

  diagnostic_summary[[dist]] <- data.frame(
    model      = dist,
    max_rhat   = max(rhat_vals, na.rm = TRUE),
    min_neff   = min(n_eff, na.rm = TRUE),
    mean_neff  = mean(n_eff, na.rm = TRUE)
  )
}

diagnostic_df <- bind_rows(diagnostic_summary)
print(diagnostic_df)
```

### Traceplots

``` r

tp <- lapply(distributions, function(dist) {
  traceplot(fits[[dist]], pars = c("mu0", "log_tau", "log_phi"),
            main = paste("Traceplots:", dist))
})
names(tp) <- distributions

tp[["lognormal"]] / tp[["weibull"]]
```

## Posterior Predictive Summaries

``` r

get_posterior_summaries <- function(fit) {
  sims <- rstan::extract(fit)
  data.frame(
    mean_mean   = mean(sims$pred_mean),
    mean_low    = quantile(sims$pred_mean,   0.025),
    mean_high   = quantile(sims$pred_mean,   0.975),
    median_mean = mean(sims$pred_median),
    median_low  = quantile(sims$pred_median, 0.025),
    median_high = quantile(sims$pred_median, 0.975),
    q90_mean    = mean(sims$pred_q90),
    q90_low     = quantile(sims$pred_q90,    0.025),
    q90_high    = quantile(sims$pred_q90,    0.975),
    q95_mean    = mean(sims$pred_q95),
    q95_low     = quantile(sims$pred_q95,    0.025),
    q95_high    = quantile(sims$pred_q95,    0.975)
  )
}

posterior_summaries <- bind_rows(
  lapply(names(fits), function(nm) {
    cbind(model = nm, get_posterior_summaries(fits[[nm]]))
  })
)

cat("\n========== Posterior Predictive Summaries ==========\n")
print(posterior_summaries)
```

## Parameter Summaries

``` r

summarise_parameters_stan <- function(fit, model_name) {
  sims <- rstan::extract(fit)

  mu0_draws <- sims$mu0
  tau_draws <- exp(sims$log_tau)
  phi_draws <- exp(sims$log_phi)

  data.frame(
    model     = model_name,
    parameter = c("mu0", "tau", "phi"),
    mean      = c(mean(mu0_draws), mean(tau_draws), mean(phi_draws)),
    median    = c(median(mu0_draws), median(tau_draws), median(phi_draws)),
    low_2.5   = c(quantile(mu0_draws, 0.025),
                  quantile(tau_draws, 0.025),
                  quantile(phi_draws, 0.025)),
    high_97.5 = c(quantile(mu0_draws, 0.975),
                  quantile(tau_draws, 0.975),
                  quantile(phi_draws, 0.975))
  )
}

param_summaries <- bind_rows(
  lapply(names(fits), function(nm) {
    summarise_parameters_stan(fits[[nm]], nm)
  })
)

cat("\n========== Parameter Summaries ==========\n")
print(param_summaries)
```

## Visualisations

### Posterior Predictive Summaries Plot

``` r

plot_data <- posterior_summaries %>%
  select(model, median_mean, median_low, median_high,
         q90_mean, q90_low, q90_high,
         q95_mean, q95_low, q95_high) %>%
  pivot_longer(cols = -model,
               names_to  = c("statistic", "measure"),
               names_pattern = "(.*)_(.*)",
               values_to = "value") %>%
  pivot_wider(names_from = measure, values_from = value)

p1 <- ggplot(plot_data, aes(x = model, y = mean, color = statistic)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = low, ymax = high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  scale_color_npg() +
  labs(title    = "Posterior Predictive Summaries by Model",
       subtitle = "Error bars show 95% credible intervals",
       x = "Distribution", y = "Value") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
```

### Parameter Estimates Plot

``` r

param_plot_data <- param_summaries %>%
  mutate(parameter = factor(parameter, levels = c("mu0", "tau", "phi")))

p2 <- ggplot(param_plot_data, aes(x = model, y = mean, color = model)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = low_2.5, ymax = high_97.5), width = 0.2) +
  facet_wrap(~parameter, scales = "free_y", ncol = 1) +
  scale_color_npg() +
  labs(title    = "Parameter Estimates by Model",
       subtitle = "Error bars show 95% credible intervals",
       x = "Distribution", y = "Parameter Value") +
  theme_minimal() +
  theme(legend.position = "none")

print(p2)
```

### Posterior Predictive CDFs

The
[`compute_predictive_cdf()`](https://cm401.github.io/ddsynth/reference/compute_predictive_cdf.md)
helper integrates over the posterior and the between-study random effect
to produce a predictive CDF with credible bands.

``` r

cdf_summaries <- list()

cdf_mats <- list()

for (dist in names(fits)) {
  cat("Computing CDF for", dist, "...\n")
  cdf_result          <- compute_predictive_cdf(
    fits[[dist]],
    dist,
    x_seq   = seq(0, 30, length.out = 500),
    n_draws = 500,
    L       = 2000
  )
  cdf_summaries[[dist]] <- cdf_result$summary
  cdf_mats[[dist]]      <- cdf_result$cdf_mat
}

cdf_df <- bind_rows(cdf_summaries)

# Extract key quantiles from the predictive CDF.
# Passing cdf_mat ensures the PI bounds are derived from the same posterior
# draws as the ribbon, making the dashed/dotted segments fully consistent
# with the shaded credible band.
quantile_summaries <- lapply(names(cdf_summaries), function(nm) {
  cbind(model = nm, extract_quantiles(cdf_summaries[[nm]],
                                      cdf_mat = cdf_mats[[nm]]))
})
quantile_df <- bind_rows(quantile_summaries)

cat("\n========== Extracted Quantiles from Predictive CDF ==========\n")
print(quantile_df)

# ---- Observed data summary ----
observed_summary <- do.call(rbind, lapply(names(datasets), function(nm) {
  d <- datasets[[nm]]

  has_median <- !is.null(d$median) && !is.na(d$median)
  has_mean   <- !is.null(d$mean)   && !is.na(d$mean)

  if (has_median) {
    center  <- d$median
    has_q1  <- !is.null(d$Q1) && !is.na(d$Q1)
    has_q3  <- !is.null(d$Q3) && !is.na(d$Q3)
    has_min <- !is.null(d$min) && !is.na(d$min)
    has_max <- !is.null(d$max) && !is.na(d$max)

    if (has_q1 && has_q3) {
      lower <- d$Q1; upper <- d$Q3; stat_type <- "Median (IQR)"
    } else if (has_min && has_max) {
      lower <- d$min; upper <- d$max; stat_type <- "Median (Range)"
    } else {
      lower <- center; upper <- center; stat_type <- "Median"
    }
  } else if (has_mean) {
    center    <- d$mean
    sd_val    <- if (!is.null(d$sd) && !is.na(d$sd)) d$sd else 0
    lower     <- center - sd_val
    upper     <- center + sd_val
    stat_type <- "Mean (\u00b1SD)"
  } else {
    return(NULL)
  }

  data.frame(dataset = nm, center_obs = center,
             lower_obs = lower, upper_obs = upper,
             stat_type = stat_type, n = d$n,
             stringsAsFactors = FALSE)
}))

observed_summary_compact <- observed_summary %>%
  arrange(center_obs) %>%
  mutate(y_pos = seq(0.1, 0.9, length.out = n()))

# ---- Observed data panel ----
p_obs <- ggplot(observed_summary_compact, aes(x = center_obs, y = y_pos)) +
  geom_point(aes(shape = stat_type, fill = stat_type),
             color = "black", size = 3, stroke = 0.8) +
  geom_errorbarh(aes(xmin = lower_obs, xmax = upper_obs, color = stat_type),
                 height = 0.05, linewidth = 0.7) +
  geom_text(aes(x = upper_obs + 0.5, label = dataset),
            size = 2.5, hjust = 0) +
  scale_shape_manual(
    values = c("Median (IQR)" = 21, "Mean (\u00b1SD)" = 22,
               "Median (Range)" = 23)) +
  scale_fill_manual(
    values = c("Median (IQR)" = "white", "Mean (\u00b1SD)" = "gray70",
               "Median (Range)" = "gray40")) +
  scale_color_manual(
    values = c("Median (IQR)" = "black", "Mean (\u00b1SD)" = "gray30",
               "Median (Range)" = "black")) +
  scale_x_continuous(limits = range(cdf_df$x)) +
  labs(y = "Observed", x = "Incubation period (days)",
       shape = "Statistic", fill = "Statistic", color = "Statistic") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), legend.position = "right")

# ---- Per-model CDF panels ----
model_colors <- c(lognormal = "darkblue", gamma = "darkgreen", weibull = "darkred")

plot_model_cdf <- function(model_name, model_color) {
  ggplot() +
    geom_ribbon(
      data = cdf_df %>% filter(model == model_name),
      aes(x = x, y = median, ymin = low, ymax = high),
      fill = model_color, alpha = 0.2
    ) +
    geom_line(
      data = cdf_df %>% filter(model == model_name),
      aes(x = x, y = median),
      color = model_color, linewidth = 1.2
    ) +
    geom_segment(
      data = quantile_df %>% filter(model == model_name, quantile == 0.5),
      aes(x = x_low, xend = x_high, y = quantile, yend = quantile),
      color = model_color, linetype = "dashed", linewidth = 1
    ) +
    geom_segment(
      data = quantile_df %>% filter(model == model_name, quantile == 0.95),
      aes(x = x_low, xend = x_high, y = quantile, yend = quantile),
      color = model_color, linetype = "dotted", linewidth = 1
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(title = tools::toTitleCase(model_name),
         y = "Cumulative probability", x = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
}

p_lognormal <- plot_model_cdf("lognormal", model_colors["lognormal"])
p_weibull   <- plot_model_cdf("weibull",   model_colors["weibull"])
p_gamma     <- plot_model_cdf("gamma",   model_colors["gamma"])

p_final <- (p_lognormal / p_gamma / p_weibull / p_obs) +
  plot_layout(heights = c(2, 2, 2, 1.5)) +
  plot_annotation(
    title    = "Posterior Predictive CDFs with Observed Data by Model",
    subtitle = paste("Dashed: predicted median",
                     "| Dotted: predicted 95th percentile",
                     "| Shaded: 95% CI"),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

print(p_final)
```

## Session Info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.0    tools_4.6.0       ragg_1.5.2        evaluate_1.0.5   
#> [21] bslib_0.10.0      yaml_2.3.12       jsonlite_2.0.0    rlang_1.2.0      
#> [25] fs_2.1.0
```
