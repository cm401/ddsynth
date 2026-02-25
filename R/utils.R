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