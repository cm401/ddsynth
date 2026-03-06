# Internal helper utilities shared across methods.
# Functions defined here are not exported (prefix name with a dot, or use
# @noRd in the roxygen block to keep them package-internal).

#' Check that a numeric vector is strictly positive
#'
#' @param x A numeric vector.
#' @param arg Name of the argument (used in the error message).
#' @return `x` invisibly, or throws an error.
#' @noRd
check_positive <- function(x, arg = deparse(substitute(x))) {
  if (!is.numeric(x) || any(x <= 0, na.rm = TRUE)) {
    stop(sprintf("'%s' must be a numeric vector of strictly positive values.", arg),
         call. = FALSE)
  }
  invisible(x)
}

#' Check that a value is a single finite number
#'
#' @param x A scalar.
#' @param arg Name of the argument (used in the error message).
#' @return `x` invisibly, or throws an error.
#' @noRd
check_scalar <- function(x, arg = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric value.", arg),
         call. = FALSE)
  }
  invisible(x)
}


# -------------------------- 
# Function to compute predictive CDF
# -------------------------- 

#' Compute posterior predictive CDF from a fitted Stan model
#'
#' Integrates over posterior draws and between-study random effects to produce
#' a predictive CDF with pointwise credible bands.
#'
#' @param fit A `stanfit` object returned by [rstan::sampling()].
#' @param dist_name Character string: `"lognormal"`, `"gamma"`, or `"weibull"`.
#' @param x_seq Numeric vector of evaluation points (default: 500 points on
#'   `[0, 30]`).
#' @param n_draws Number of posterior draws to use (default: 500).
#' @param L Number of study-level locations to integrate over per draw
#'   (default: 50).
#'
#' @return A data frame with columns `x`, `median`, `mean`, `low`, `high`, and
#'   `model`.
#' @export
compute_predictive_cdf <- function(fit, dist_name, x_seq = seq(0, 30, length.out = 500), 
                                   n_draws = 500, L = 50) {
  
  # Extract posterior samples
  sims <- rstan::extract(fit)
  
  # Sample from posterior
  n_post <- length(sims$mu0)
  draws_idx <- sample(1:n_post, min(n_draws, n_post))
  
  # Storage for CDF values
  cdf_mat <- matrix(NA, nrow = length(draws_idx), ncol = length(x_seq))
  
  for (i in seq_along(draws_idx)) {
    idx <- draws_idx[i]
    mu0 <- sims$mu0[idx]
    tau <- sims$tau[idx]  
    phi <- sims$phi[idx]  
    
    # Integrate over L study-level locations
    locs <- rnorm(L, mean = mu0, sd = tau)
    
    # Compute CDF for each location and average
    cdf_l <- matrix(NA, nrow = L, ncol = length(x_seq))
    
    for (l in 1:L) {
      loc_d <- locs[l]
      
      if (dist_name == "lognormal") {
        cdf_l[l, ] <- plnorm(x_seq, meanlog = loc_d, sdlog = phi)
        
      } else if (dist_name == "gamma") {
        mean_d <- exp(loc_d)
        shape <- phi
        rate <- shape / mean_d
        cdf_l[l, ] <- pgamma(x_seq, shape = shape, rate = rate)
        
      } else if (dist_name == "weibull") {
        scale <- exp(loc_d)
        shape <- phi
        cdf_l[l, ] <- pweibull(x_seq, shape = shape, scale = scale)
      }
    }
    
    # Average over study-level locations
    cdf_mat[i, ] <- colMeans(cdf_l)
  }
  
  # Compute summary statistics
  data.frame(
    x = x_seq,
    median = apply(cdf_mat, 2, median, na.rm = TRUE),
    mean = apply(cdf_mat, 2, mean, na.rm = TRUE),
    low = apply(cdf_mat, 2, quantile, 0.025, na.rm = TRUE),
    high = apply(cdf_mat, 2, quantile, 0.975, na.rm = TRUE),
    model = dist_name
  )
}

# -------------------------- 
# Extract quantiles from CDF
# -------------------------- 

#' Extract quantiles from a predictive CDF summary
#'
#' For each requested probability, finds the x value where the CDF (and its
#' credible bounds) crosses that probability.
#'
#' @param cdf_summary A data frame produced by [compute_predictive_cdf()],
#'   with columns `x`, `median`, `low`, and `high`.
#' @param probs Numeric vector of probabilities to extract (default:
#'   `c(0.5, 0.95)`).
#'
#' @return A data frame with columns `quantile`, `quantile_label`, `x_median`,
#'   `x_low`, and `x_high`.
#' @export
extract_quantiles <- function(cdf_summary, probs = c(0.5, 0.95)) {  # CHANGED: added 0.95
  results <- list()
  
  for (p in probs) {
    # Find x value where CDF crosses probability p
    idx_median <- which.min(abs(cdf_summary$median - p))
    idx_low <- which.min(abs(cdf_summary$low - p))
    idx_high <- which.min(abs(cdf_summary$high - p))
    
    results[[paste0("q", p*100)]] <- data.frame(
      quantile = p,
      quantile_label = paste0("Q", p*100),  # NEW: for labeling
      x_median = cdf_summary$x[idx_median],
      x_low = cdf_summary$x[idx_low],
      x_high = cdf_summary$x[idx_high]
    )
  }
  
  dplyr::bind_rows(results)
}


#' Prepare Stan data from a list of dataset summaries
#'
#' Converts a list of dataset descriptors (each providing summary statistics
#' and a sample size) into the named list expected by the
#' `hierarchical_data_synthesis_summary_stats` Stan model.
#'
#' @param datasets A named list of lists.  Each element must contain `n` (sample
#'   size) and one of the following combinations of summary statistics:
#'   \describe{
#'     \item{`median`, `min`, `max`}{Median and range (summary type 1).}
#'     \item{`median`, `Q1`, `Q3`}{Median and inter-quartile range (summary type 2).}
#'     \item{`mean`, `sd`}{Mean and standard deviation (summary type 3).}
#'   }
#' @param dist_type Integer distribution code: `1` = log-normal, `2` = gamma,
#'   `3` = Weibull. Defaults to `1`.
#' @param use_custom_priors Integer flag (0 or 1) for custom prior use.
#'   Currently unused; reserved for future extension. Defaults to `0`.
#' @param custom_priors Optional list of custom prior values. Currently unused.
#'
#' @return A named list suitable for passing to [rstan::sampling()] as the
#'   `data` argument.
#' @export
prepare_stan_data_from_datasets <- function(datasets, dist_type = 1,
                                            use_custom_priors = 0,
                                            custom_priors = NULL) {
  
  n_datasets <- length(datasets)
  
  # Initialize vectors
  n_obs_vec <- integer(n_datasets)
  summary_type <- integer(n_datasets)
  obs_stat1 <- numeric(n_datasets)
  obs_stat2 <- numeric(n_datasets)
  obs_stat3 <- numeric(n_datasets)
  
  # Process each dataset
  for (i in seq_along(datasets)) {
    d <- datasets[[i]]
    n_obs_vec[i] <- d$n
    
    # Determine summary type and extract statistics
    if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) {
      # Type 1: median + range (min, max)
      summary_type[i] <- 1
      obs_stat1[i] <- d$median
      obs_stat2[i] <- d$min
      obs_stat3[i] <- d$max
      
    } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) {
      # Type 2: median + IQR (Q1, Q3)
      summary_type[i] <- 2
      obs_stat1[i] <- d$median
      obs_stat2[i] <- d$Q1
      obs_stat3[i] <- d$Q3
      
    } else if (!is.null(d$mean) && !is.null(d$sd)) {
      # Type 3: mean + sd
      summary_type[i] <- 3
      obs_stat1[i] <- d$mean
      obs_stat2[i] <- d$sd
      obs_stat3[i] <- 0  # placeholder
      
    } else {
      stop(paste("Dataset", i, "does not have recognized summary statistics"))
    }
  }
  
  # Create base Stan data
  stan_data <- list(
    n_datasets   = n_datasets,
    n_obs        = as.array(n_obs_vec),
    summary_type = as.array(summary_type),
    dist_type    = dist_type,
    obs_stat1    = as.array(obs_stat1),
    obs_stat2    = as.array(obs_stat2),
    obs_stat3    = as.array(obs_stat3)
  )
  
  stan_data$mu0_mean <- log(mean(obs_stat1))  # log of overall central estimates as prior mean for mu0
  stan_data$mu0_sd <- 1.0
  stan_data$log_tau_mean <- 0.2
  stan_data$log_tau_sd <- 0.5
  stan_data$log_phi_mean <- 0.2
  stan_data$log_phi_sd <- 0.5
  
  return(stan_data)
}



# Function for simulation studies -----------------------------------------


#' Create a scenario with specific characteristics
#'
#' @param scenario_name Descriptive name for the scenario
#' @param dist_type Distribution type
#' @param n_datasets Number of datasets
#' @param n_obs_config Configuration for sample sizes: "fixed", "small_var", "large_var", "custom"
#' @param n_obs_values Custom vector of sample sizes (if n_obs_config = "custom")
#' @param summary_config Configuration for summary types: "fixed", "mixed_balanced", "mixed_random", "custom"
#' @param summary_values Custom vector of summary types (if summary_config = "custom")
#' @param fixed_summary_type Summary type to use when summary_config = "fixed" (default 1).
#' @param fixed_n_obs Sample size to use when n_obs_config = "fixed" (default 14).
#' @param mu0 Population mean
#' @param tau Between-study SD
#' @param phi Distribution-specific parameter
#' @return Scenario specification list
#' @export
create_scenario <- function(scenario_name,
                            dist_type,
                            n_datasets,
                            n_obs_config = c("fixed", "small_var", "large_var", "custom"),
                            n_obs_values = NULL,
                            summary_config = c("fixed", "mixed_balanced", "mixed_random", "custom"),
                            summary_values = NULL,
                            fixed_summary_type = 1,
                            fixed_n_obs = 14,
                            mu0 = 2.0,
                            tau = 0.4,
                            phi = NULL) {
  
  n_obs_config <- match.arg(n_obs_config)
  summary_config <- match.arg(summary_config)
  
  # Generate sample sizes based on configuration
  if (n_obs_config == "fixed") {
    n_obs <- rep(fixed_n_obs, n_datasets)
    
  } else if (n_obs_config == "small_var") {
    # Small variation: +/-20% around mean
    mean_n <- fixed_n_obs
    n_obs <- round(runif(n_datasets, mean_n * 0.6, mean_n * 1.4))
    n_obs <- pmax(n_obs, 5)  # Minimum of 10
    
  } else if (n_obs_config == "large_var") {
    # Large variation: realistic range from small to large studies
    n_obs <- sample(c(5,10, 15, 20, 30, 50, 75, 100), 
                    n_datasets, replace = TRUE)
    
  } else if (n_obs_config == "custom") {
    if (is.null(n_obs_values) || length(n_obs_values) != n_datasets) {
      stop("For custom n_obs_config, must provide n_obs_values with length = n_datasets")
    }
    n_obs <- n_obs_values
  }
  
  # Generate summary types based on configuration
  if (summary_config == "fixed") {
    summary_type <- rep(fixed_summary_type, n_datasets)
    
  } else if (summary_config == "mixed_balanced") {
    # Equal representation of all three types
    summary_type <- rep(1:3, length.out = n_datasets)
    summary_type <- sample(summary_type)  # Shuffle
    
  } else if (summary_config == "mixed_random") {
    # Random mix with realistic probabilities
    # Median+range more common, mean+sd less common
    summary_type <- sample(1:3, n_datasets, replace = TRUE, 
                           prob = c(0.5, 0.2, 0.3))
    
  } else if (summary_config == "custom") {
    if (is.null(summary_values) || length(summary_values) != n_datasets) {
      stop("For custom summary_config, must provide summary_values with length = n_datasets")
    }
    summary_type <- summary_values
  }
  
  # Set default phi if not provided
  if (is.null(phi)) {
    phi <- switch(dist_type,
                  "lognormal" = 0.35,
                  "gamma" = 3.5,
                  "weibull" = 2.0)
  }
  
  list(
    scenario_name = scenario_name,
    dist_type = dist_type,
    n_datasets = n_datasets,
    n_obs = n_obs,
    summary_type = summary_type,
    mu0 = mu0,
    tau = tau,
    phi = phi
  )
}



#' Generate a comprehensive set of scenarios
#'
#' @param include_homogeneous Include scenarios with fixed summary types
#' @param include_mixed Include scenarios with mixed summary types
#' @param include_varied_n Include scenarios with varied sample sizes
#' @return Data frame of scenarios
#' @export
generate_scenario_library <- function(include_homogeneous = TRUE,
                                      include_mixed = TRUE,
                                      include_varied_n = TRUE) {
  
  scenarios <- list()
  idx <- 1
  
  distributions <- c("lognormal", "gamma", "weibull")
  
  # ===== HOMOGENEOUS SCENARIOS =====
  if (include_homogeneous) {
    for (dist in distributions) {
      for (summary_type in 1:3) {
        for (n_datasets in c(5, 10, 20)) {
          for (n_obs in c(5, 10, 20, 50)) {
            scenarios[[idx]] <- create_scenario(
              scenario_name = sprintf("Homog_%s_ST%d_D%d_N%d", 
                                      dist, summary_type, n_datasets, n_obs),
              dist_type = dist,
              n_datasets = n_datasets,
              n_obs_config = "fixed",
              fixed_n_obs = n_obs,
              summary_config = "fixed",
              fixed_summary_type = summary_type
            )
            idx <- idx + 1
          }
        }
      }
    }
  }
  
  # ===== MIXED SUMMARY TYPE SCENARIOS =====
  if (include_mixed) {
    for (dist in distributions) {
      for (n_datasets in c(15, 30)) {
        # Balanced mix
        scenarios[[idx]] <- create_scenario(
          scenario_name = sprintf("Mixed_Balanced_%s_D%d", dist, n_datasets),
          dist_type = dist,
          n_datasets = n_datasets,
          n_obs_config = "fixed",
          fixed_n_obs = 30,
          summary_config = "mixed_balanced"
        )
        idx <- idx + 1
        
        # Random mix
        scenarios[[idx]] <- create_scenario(
          scenario_name = sprintf("Mixed_Random_%s_D%d", dist, n_datasets),
          dist_type = dist,
          n_datasets = n_datasets,
          n_obs_config = "fixed",
          fixed_n_obs = 30,
          summary_config = "mixed_random"
        )
        idx <- idx + 1
      }
    }
  }
  
  # ===== VARIED SAMPLE SIZE SCENARIOS =====
  if (include_varied_n) {
    for (dist in distributions) {
      for (summary_config in c("fixed", "mixed_balanced")) {
        # Small variation in sample size
        scenarios[[idx]] <- create_scenario(
          scenario_name = sprintf("VarN_Small_%s_%s", 
                                  dist, 
                                  ifelse(summary_config == "fixed", "ST1", "Mixed")),
          dist_type = dist,
          n_datasets = 20,
          n_obs_config = "small_var",
          fixed_n_obs = 40,
          summary_config = summary_config,
          fixed_summary_type = 1
        )
        idx <- idx + 1
        
        # Large variation in sample size
        scenarios[[idx]] <- create_scenario(
          scenario_name = sprintf("VarN_Large_%s_%s", 
                                  dist,
                                  ifelse(summary_config == "fixed", "ST1", "Mixed")),
          dist_type = dist,
          n_datasets = 20,
          n_obs_config = "large_var",
          summary_config = summary_config,
          fixed_summary_type = 1
        )
        idx <- idx + 1
      }
    }
  }
  
  # Convert to data frame for easier handling
  scenarios_df <- dplyr::bind_rows(lapply(scenarios, function(s) {
    data.frame(
      scenario_name = s$scenario_name,
      dist_type = s$dist_type,
      n_datasets = s$n_datasets,
      mu0 = s$mu0,
      tau = s$tau,
      phi = s$phi,
      n_obs_mean = mean(s$n_obs),
      n_obs_sd = sd(s$n_obs),
      n_obs_min = min(s$n_obs),
      n_obs_max = max(s$n_obs),
      summary_type_1_prop = mean(s$summary_type == 1),
      summary_type_2_prop = mean(s$summary_type == 2),
      summary_type_3_prop = mean(s$summary_type == 3),
      stringsAsFactors = FALSE
    )
  }))
  
  # Store full scenario details as attribute
  attr(scenarios_df, "full_scenarios") <- scenarios
  
  return(scenarios_df)
}

#' Fit Stan model to simulated data
#'
#' @param sim_data Simulated data from generate_hierarchical_data
#' @param stan_model Compiled Stan model
#' @param ... Additional arguments to pass to sampling()
#' @return Stan fit object
#' @export
fit_model <- function(sim_data, stan_model, ...) {
  rstan::sampling(
    stan_model,
    data = sim_data$obs_data,
    chains = 4,
    iter = 10000,
    warmup = 1000,
    ...
  )
}

#' Compute coverage for a parameter
#'
#' @param fit Stan fit object
#' @param param_name Parameter name
#' @param true_value True parameter value
#' @param level Credible interval level (default 0.95)
#' @return Logical indicating whether true value is in credible interval
#' @export
check_coverage <- function(fit, param_name, true_value, level = 0.95) {
  draws <- rstan::extract(fit)
  
  if (param_name %in% names(draws)) {
    param_draws <- draws[[param_name]]
    ci <- quantile(param_draws, probs = c((1-level)/2, 1-(1-level)/2))
    return(true_value >= ci[1] & true_value <= ci[2])
  } else {
    warning(paste("Parameter", param_name, "not found"))
    return(NA)
  }
}

#' Compute median bias
#'
#' @param fit Stan fit object
#' @param param_name Parameter name
#' @param true_value True parameter value
#' @return Median bias (median estimate - true value)
#' @export
compute_median_bias <- function(fit, param_name, true_value) {
  draws <- rstan::extract(fit)
  
  if (param_name %in% names(draws)) {
    param_draws <- draws[[param_name]]
    return(median(param_draws) - true_value)
  } else {
    warning(paste("Parameter", param_name, "not found"))
    return(NA)
  }
}

#' Compute Integrated Quadratic Distance (IQD)
#'
#' @param fit Stan fit object
#' @param true_params List of true parameters
#' @param dist_type Distribution type
#' @param x_grid Grid of x values for integration
#' @return IQD value
#' @export
compute_iqd <- function(fit, true_params, dist_type, x_grid = NULL) {

  draws <- rstan::extract(fit)
  
  # Extract posterior samples
  mu0_samples <- draws$mu0
  tau_samples <- exp(draws$log_tau)
  phi_samples <- exp(draws$log_phi)
  
  n_samples <- length(mu0_samples)
  
  # Create grid if not provided
  if (is.null(x_grid)) {
    if (dist_type == "lognormal") {
      x_grid <- seq(0.01, exp(true_params$mu0 + 3*true_params$tau), length.out = 200)
    } else if (dist_type == "gamma") {
      mean_max <- exp(true_params$mu0 + 3*true_params$tau)
      x_grid <- seq(0.01, mean_max * 3, length.out = 200)
    } else if (dist_type == "weibull") {
      scale_max <- exp(true_params$mu0 + 3*true_params$tau)
      x_grid <- seq(0.01, scale_max * 3, length.out = 200)
    }
  }
  
  # Compute true predictive density
  true_density <- numeric(length(x_grid))
  for (i in seq_along(x_grid)) {
    x <- x_grid[i]
    
    # Integrate over random effects distribution
    integrand <- function(loc) {
      if (dist_type == "lognormal") {
        dlnorm(x, meanlog = loc, sdlog = true_params$phi) * 
          dnorm(loc, mean = true_params$mu0, sd = true_params$tau)
      } else if (dist_type == "gamma") {
        mean_d <- exp(loc)
        shape <- true_params$phi
        rate <- shape / mean_d
        dgamma(x, shape = shape, rate = rate) * 
          dnorm(loc, mean = true_params$mu0, sd = true_params$tau)
      } else if (dist_type == "weibull") {
        scale <- exp(loc)
        shape <- true_params$phi
        dweibull(x, shape = shape, scale = scale) * 
          dnorm(loc, mean = true_params$mu0, sd = true_params$tau)
      }
    }
    
    true_density[i] <- integrate(integrand, 
                                 lower = true_params$mu0 - 5*true_params$tau,
                                 upper = true_params$mu0 + 5*true_params$tau)$value
  }
  
  # Compute estimated predictive density (average over posterior samples)
  est_density <- numeric(length(x_grid))
  
  # Subsample for computational efficiency
  sample_idx <- sample(1:n_samples, min(500, n_samples))
  
  for (i in seq_along(x_grid)) {
    x <- x_grid[i]
    
    density_samples <- numeric(length(sample_idx))
    for (s in seq_along(sample_idx)) {
      idx <- sample_idx[s]
      
      # Integrate over random effects for this posterior sample
      integrand <- function(loc) {
        if (dist_type == "lognormal") {
          dlnorm(x, meanlog = loc, sdlog = phi_samples[idx]) * 
            dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
        } else if (dist_type == "gamma") {
          mean_d <- exp(loc)
          shape <- phi_samples[idx]
          rate <- shape / mean_d
          dgamma(x, shape = shape, rate = rate) * 
            dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
        } else if (dist_type == "weibull") {
          scale <- exp(loc)
          shape <- phi_samples[idx]
          dweibull(x, shape = shape, scale = scale) * 
            dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
        }
      }
      
      density_samples[s] <- integrate(integrand,
                                      lower = mu0_samples[idx] - 5*tau_samples[idx],
                                      upper = mu0_samples[idx] + 5*tau_samples[idx])$value
    }
    
    est_density[i] <- mean(density_samples)
  }
  
  # Compute IQD using trapezoidal rule
  dx <- diff(x_grid)
  squared_diff <- (true_density - est_density)^2
  iqd <- sum((squared_diff[-1] + squared_diff[-length(squared_diff)]) / 2 * dx)
  
  return(iqd)
}

#' Run simulation study with generalized scenarios
#'
#' @param n_sim Number of simulation replicates per scenario
#' @param scenarios_df Data frame from generate_scenario_library()
#' @param stan_model Compiled Stan model
#' @param seed Random seed
#' @param save_name File path for intermediate RDS save after each scenario.
#' @param parallel Use parallel processing (requires future package)
#' @return Data frame with one row per simulation replicate and columns for
#'   scenario metadata, true parameter values, coverage, bias, and IQD metrics.
#' @export
run_simulation_study_generalized <- function(n_sim, 
                                             scenarios_df, 
                                             stan_model, 
                                             seed = 123,
                                             save_name = "simulation_results_tmp_general.rds",
                                             parallel = FALSE) {
  
  set.seed(seed)
  
  # Get full scenario details
  full_scenarios <- attr(scenarios_df, "full_scenarios")
  
  if (is.null(full_scenarios)) {
    stop("scenarios_df must be created by generate_scenario_library()")
  }
  
  results <- list()
  result_idx <- 1
  
  for (scenario_idx in seq_len(nrow(scenarios_df))) {
    
    scenario <- scenarios_df[scenario_idx,] #full_scenarios[[scenario_idx]]
    
    cat("\n========================================\n")
    cat("Scenario", rownames(scenarios_df[scenario_idx,]), "of", length(full_scenarios), "\n")
    cat("Name:", scenario$scenario_name, "\n")
    cat("Distribution:", scenario$dist_type, "\n")
    cat("n_datasets:", scenario$n_datasets, "\n")
    cat("Sample sizes: ", paste(range(scenario$n_obs_mean), collapse = "-"), 
        " (mean:", round(mean(scenario$n_obs_mean), 1), ")\n")
    cat("Summary types:", paste(unique(scenario$summary_type), collapse = ", "), "\n")
    cat("========================================\n\n")
    
    for (sim in 1:n_sim) {
      
      cat("  Simulation", sim, "of", n_sim, "\n")
      
      if (startsWith(scenario$scenario_name, "VarN"))
      {
        n_obs_in <- pmax(pmin(round(rnorm( scenario$n_datasets, scenario$n_obs_mean, scenario$n_obs_sd )),scenario$n_obs_max),scenario$n_obs_min)         
      } else {
        n_obs_in <- scenario$n_obs_mean        
      }
      
      # Generate data using the mixed function
      sim_data <- generate_hierarchical_data_mixed(
        n_datasets = scenario$n_datasets,
        n_obs = n_obs_in,
        dist_type = scenario$dist_type,
        mu0 = scenario$mu0,
        tau = scenario$tau,
        phi = scenario$phi,
        summary_type = ifelse(scenario$summary_type_1_prop == 1, 1,
                              ifelse(scenario$summary_type_2_prop == 1, 2,
                                     ifelse(scenario$summary_type_3_prop == 1, 3, c(scenario$summary_type_1_prop,scenario$summary_type_2_prop,scenario$summary_type_3_prop))))
      )
      
      # Fit model
      tryCatch({
        fit <- fit_model(sim_data, stan_model, 
                         refresh = 0, 
                         control = list(adapt_delta = 0.95, max_treedepth = 12))
        
        # Check convergence
        fit_summary <- summary(fit)$summary
        max_rhat <- max(fit_summary[, "Rhat"], na.rm = TRUE)
        min_neff <- min(fit_summary[, "n_eff"], na.rm = TRUE)
        
        if (max_rhat > 1.1 || min_neff < 100) {
          warning(paste("Convergence issues in scenario", scenario_idx, 
                        "sim", sim, ": Rhat =", round(max_rhat, 3),
                        ", min n_eff =", round(min_neff, 0)))
        }
        
        # Compute metrics
        coverage_mu0 <- check_coverage(fit, "mu0", sim_data$true_params$mu0)
        coverage_tau <- check_coverage(fit, "tau", sim_data$true_params$tau)
        coverage_phi <- check_coverage(fit, "phi", sim_data$true_params$phi)
        
        bias_mu0 <- compute_median_bias(fit, "mu0", sim_data$true_params$mu0)
        bias_tau <- compute_median_bias(fit, "tau", sim_data$true_params$tau)
        bias_phi <- compute_median_bias(fit, "phi", sim_data$true_params$phi)
        
        # Compute relative bias
        rel_bias_mu0 <- bias_mu0 / sim_data$true_params$mu0
        rel_bias_tau <- bias_tau / sim_data$true_params$tau
        rel_bias_phi <- bias_phi / sim_data$true_params$phi
        
        iqd <- compute_iqd(fit, sim_data$true_params, scenario$dist_type)
        
        # Store results
        results[[result_idx]] <- data.frame(
          scenario_idx = scenario_idx,
          scenario_name = scenario$scenario_name,
          sim = sim,
          dist_type = scenario$dist_type,
          n_datasets = scenario$n_datasets,
          n_obs_mean = mean(scenario$n_obs),
          n_obs_sd = sd(scenario$n_obs),
          n_obs_min = min(scenario$n_obs),
          n_obs_max = max(scenario$n_obs),
          summary_type_diversity = length(unique(scenario$summary_type)),
          prop_summary_type_1 = mean(scenario$summary_type == 1),
          prop_summary_type_2 = mean(scenario$summary_type == 2),
          prop_summary_type_3 = mean(scenario$summary_type == 3),
          true_mu0 = scenario$mu0,
          true_tau = scenario$tau,
          true_phi = scenario$phi,
          coverage_mu0 = coverage_mu0,
          coverage_tau = coverage_tau,
          coverage_phi = coverage_phi,
          bias_mu0 = bias_mu0,
          bias_tau = bias_tau,
          bias_phi = bias_phi,
          rel_bias_mu0 = rel_bias_mu0,
          rel_bias_tau = rel_bias_tau,
          rel_bias_phi = rel_bias_phi,
          iqd = iqd,
          max_rhat = max_rhat,
          min_neff = min_neff,
          converged = max_rhat <= 1.1 & min_neff >= 100,
          stringsAsFactors = FALSE
        )
        
        result_idx <- result_idx + 1
        
      }, error = function(e) {
        warning(paste("Error in scenario", scenario_idx, "sim", sim, ":", e$message))

        results[[result_idx]] <<- data.frame(
          scenario_idx = scenario_idx,
          scenario_name = scenario$scenario_name,
          sim = sim,
          dist_type = scenario$dist_type,
          n_datasets = scenario$n_datasets,
          n_obs_mean = mean(scenario$n_obs),
          n_obs_sd = sd(scenario$n_obs),
          n_obs_min = min(scenario$n_obs),
          n_obs_max = max(scenario$n_obs),
          summary_type_diversity = length(unique(scenario$summary_type)),
          prop_summary_type_1 = mean(scenario$summary_type == 1),
          prop_summary_type_2 = mean(scenario$summary_type == 2),
          prop_summary_type_3 = mean(scenario$summary_type == 3),
          true_mu0 = scenario$mu0,
          true_tau = scenario$tau,
          true_phi = scenario$phi,
          coverage_mu0 = NA,
          coverage_tau = NA,
          coverage_phi = NA,
          bias_mu0 = NA,
          bias_tau = NA,
          bias_phi = NA,
          rel_bias_mu0 = NA,
          rel_bias_tau = NA,
          rel_bias_phi = NA,
          iqd = NA,
          max_rhat = NA,
          min_neff = NA,
          converged = FALSE,
          stringsAsFactors = FALSE
        )

        result_idx <<- result_idx + 1
      })
    }
    
    saveRDS(results, save_name)
  }
  
  # Combine results
  dplyr::bind_rows(results)
}

