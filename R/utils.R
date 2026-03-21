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
#'   (default: 2000). When `n_datasets < 5`, `mean(loc_d)` is used directly
#'   for all `L` locations (i.e. no between-study sampling) for consistency
#'   with the Stan generated quantities block — see
#'   [prepare_stan_data_from_datasets()] for details.
#'
#' @return A data frame with columns `x`, `median`, `mean`, `low`, `high`, and
#'   `model`.
#' @export
compute_predictive_cdf <- function(fit, dist_name, x_seq = seq(0, 30, length.out = 500),
                                   n_draws = 500, L = 2000) {

  # Extract posterior samples
  sims <- rstan::extract(fit)

  # Determine number of datasets from the study-level location parameter
  n_datasets <- dim(sims$loc_d)[2]

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

    # Integrate over L study-level locations.
    # When n_datasets < 5, mu0 is confounded with tau * loc_d_raw and is not
    # directly identified by the data — only their sum (loc_d) is. Using mu0
    # therefore produces overly wide predictive intervals. Instead we use
    # mean(loc_d_draw), the mean of the study-level location estimates, which
    # is the quantity the data actually constrains. This is consistent with
    # the Stan generated quantities block.
    #   n_datasets == 1 : mean(loc_d) == loc_d[1], tightly identified
    #   n_datasets 2-4  : sample mean of study-level estimates
    if (n_datasets < 5) {
      loc_d_draw <- sims$loc_d[idx, ]
      locs <- rep(mean(loc_d_draw), L)
    } else {
      locs <- rnorm(L, mean = mu0, sd = tau)
    }

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
#'     \item{`freq_value`, `freq_count`}{Frequency table of (value, count) pairs
#'       (summary type 4). `n` is optional and defaults to `sum(freq_count)`.}
#'     \item{`freq_lower`, `freq_upper`, `freq_count`}{Interval-censored frequency
#'       table (summary type 5). Each entry gives the lower and upper bound of the
#'       censoring interval and the count of individuals in that interval.  When
#'       `freq_lower[i] == freq_upper[i]` the observation is treated as exact.
#'       `n` is optional and defaults to `sum(freq_count)`.}
#'   }
#' @param dist_type Integer distribution code: `1` = log-normal, `2` = gamma,
#'   `3` = Weibull. Defaults to `1`.
#' @param use_custom_priors Integer flag (0 or 1) for custom prior use.
#'   Currently unused; reserved for future extension. Defaults to `0`.
#' @param custom_priors Named list of prior overrides. Any values not supplied
#'   fall back to distribution-appropriate defaults (see Details). Recognised
#'   names: `mu0_sd`, `log_tau_mean`, `log_tau_sd`, `log_phi_mean`,
#'   `log_phi_sd`.
#'
#' @details
#' **Distribution-specific defaults for `log_phi`:**
#'
#' Because `phi` has a different meaning in each distribution, the default
#' prior for `log_phi_mean` is chosen per `dist_type`:
#'
#' | `dist_type` | Distribution | `phi` | Default `log_phi_mean` | Prior median phi |
#' |---|---|---|---|---|
#' | 1 | Lognormal | log-SD (σ) | -0.7 | 0.50 |
#' | 2 | Gamma | shape | 2.5 | 12.2 |
#' | 3 | Weibull | shape | 1.0 | 2.7 |
#'
#' Users can override any individual prior by passing only the relevant
#' element(s) in `custom_priors`, e.g.
#' `custom_priors = list(log_phi_mean = 3.0)` — all other priors will use
#' the distribution-appropriate defaults above.
#'
#' @return A named list suitable for passing to [rstan::sampling()] as the
#'   `data` argument. The list always includes `freq_lower` and `freq_upper`
#'   fields (populated with zeros for non-type-5 datasets), as these are
#'   required by the Stan model regardless of which summary types are present.
#'
#' @note **Backward compatibility:** The Stan model requires `freq_lower` and
#'   `freq_upper` to be present in the data list for all runs, including those
#'   that contain only type 1--4 datasets. This is handled automatically when
#'   using this function. If you construct the Stan data list manually (rather
#'   than via this function), you must include these fields explicitly, e.g.:
#'   ```r
#'   stan_data$freq_lower <- rep(0, stan_data$n_freq_total)
#'   stan_data$freq_upper <- rep(0, stan_data$n_freq_total)
#'   ```
#' @export
prepare_stan_data_from_datasets <- function(datasets, dist_type = 1,
                                            use_custom_priors = 0,
                                            custom_priors = list()) {

  # Apply distribution-specific defaults for log_phi_mean/log_phi_sd.
  # phi has a different meaning in each distribution:
  #   lognormal : phi = log-SD (sigma),  typical range 0.2-0.7  -> log_phi_mean = -0.7
  #   gamma     : phi = shape,           typical range 5-30     -> log_phi_mean =  2.5
  #   weibull   : phi = shape,           typical range 2-6      -> log_phi_mean =  1.0
  # All other priors share the same sensible defaults regardless of dist_type.
  phi_defaults <- list(
    `1` = list(log_phi_mean = -0.7, log_phi_sd = 0.5),   # lognormal
    `2` = list(log_phi_mean =  2.5, log_phi_sd = 0.5),   # gamma
    `3` = list(log_phi_mean =  1.0, log_phi_sd = 0.5)    # weibull
  )[[as.character(dist_type)]]

  defaults <- list(
    mu0_sd       = 1.0,
    log_tau_mean = 0.2,
    log_tau_sd   = 0.5,
    log_phi_mean = phi_defaults$log_phi_mean,
    log_phi_sd   = phi_defaults$log_phi_sd
  )

  # User-supplied values in custom_priors override defaults; anything not
  # supplied falls back to the distribution-appropriate default above.
  custom_priors <- modifyList(defaults, custom_priors)

  n_datasets <- length(datasets)

  if (n_datasets < 5) {
    warning(
      "n_datasets = ", n_datasets, " (< 5): the between-study heterogeneity ",
      "parameter tau cannot be reliably identified from so few studies and will ",
      "be largely determined by its prior. Predicted quantities (pred_mean, ",
      "pred_median, etc.) are therefore computed at the population mean mu0 ",
      "rather than averaging over the predictive distribution for new studies. ",
      "See Higgins & Thompson (2002) doi:10.1002/sim.1186, ",
      "Gelman (2006) doi:10.1214/06-BA117A, ",
      "and Rover et al. (2021) doi:10.1002/jrsm.1475.",
      call. = FALSE
    )
  }

  # Initialize vectors
  n_obs_vec      <- integer(n_datasets)
  summary_type   <- integer(n_datasets)
  obs_stat1      <- numeric(n_datasets)
  obs_stat2      <- numeric(n_datasets)
  obs_stat3      <- numeric(n_datasets)
  # Frequency table flat arrays (for summary_type == 4 and 5)
  freq_value_all <- numeric(0)
  freq_lower_all <- numeric(0)
  freq_upper_all <- numeric(0)
  freq_count_all <- integer(0)
  freq_start_vec <- integer(n_datasets)
  freq_len_vec   <- integer(n_datasets)
  running_start  <- 1L

  # Process each dataset
  for (i in seq_along(datasets)) {
    d <- datasets[[i]]

    # Determine summary type and extract statistics
    if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) {
      # Type 1: median + range (min, max)
      n_obs_vec[i]    <- d$n
      summary_type[i] <- 1
      obs_stat1[i]    <- d$median
      obs_stat2[i]    <- d$min
      obs_stat3[i]    <- d$max

    } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) {
      # Type 2: median + IQR (Q1, Q3)
      n_obs_vec[i]    <- d$n
      summary_type[i] <- 2
      obs_stat1[i]    <- d$median
      obs_stat2[i]    <- d$Q1
      obs_stat3[i]    <- d$Q3

    } else if (!is.null(d$mean) && !is.null(d$sd)) {
      # Type 3: mean + sd
      n_obs_vec[i]    <- d$n
      summary_type[i] <- 3
      obs_stat1[i]    <- d$mean
      obs_stat2[i]    <- d$sd
      obs_stat3[i]    <- 0  # placeholder

    } else if (!is.null(d$freq_value) && !is.null(d$freq_count)) {
      # Type 4: frequency table (e.g. delays rounded to whole days)
      n_obs_vec[i]      <- if (!is.null(d$n)) d$n else sum(d$freq_count)
      summary_type[i]   <- 4
      obs_stat1[i]      <- 0  # placeholder
      obs_stat2[i]      <- 0  # placeholder
      obs_stat3[i]      <- 0  # placeholder
      freq_start_vec[i] <- running_start
      freq_len_vec[i]   <- length(d$freq_value)
      freq_value_all    <- c(freq_value_all, as.numeric(d$freq_value))
      freq_lower_all    <- c(freq_lower_all, rep(0, length(d$freq_value)))  # unused for type 4
      freq_upper_all    <- c(freq_upper_all, rep(0, length(d$freq_value)))  # unused for type 4
      freq_count_all    <- c(freq_count_all, as.integer(d$freq_count))
      running_start     <- running_start + freq_len_vec[i]

    } else if (!is.null(d$freq_lower) && !is.null(d$freq_upper) && !is.null(d$freq_count)) {
      # Type 5: interval-censored frequency table
      if (length(d$freq_lower) != length(d$freq_upper) ||
          length(d$freq_lower) != length(d$freq_count)) {
        stop(paste("Dataset", i, ": freq_lower, freq_upper and freq_count must all have the same length"))
      }
      if (any(d$freq_lower > d$freq_upper)) {
        stop(paste("Dataset", i, ": all freq_lower values must be <= their corresponding freq_upper values"))
      }
      n_obs_vec[i]      <- if (!is.null(d$n)) d$n else sum(d$freq_count)
      summary_type[i]   <- 5
      obs_stat1[i]      <- 0  # placeholder
      obs_stat2[i]      <- 0  # placeholder
      obs_stat3[i]      <- 0  # placeholder
      freq_start_vec[i] <- running_start
      freq_len_vec[i]   <- length(d$freq_lower)
      freq_value_all    <- c(freq_value_all, rep(0, length(d$freq_lower)))  # unused for type 5
      freq_lower_all    <- c(freq_lower_all, as.numeric(d$freq_lower))
      freq_upper_all    <- c(freq_upper_all, as.numeric(d$freq_upper))
      freq_count_all    <- c(freq_count_all, as.integer(d$freq_count))
      running_start     <- running_start + freq_len_vec[i]

    } else {
      stop(paste("Dataset", i, "does not have recognized summary statistics"))
    }
  }

  # Compute central estimates for mu0 prior (use weighted mean from freq table for types 4 and 5)
  central_estimates <- numeric(n_datasets)
  for (i in seq_len(n_datasets)) {
    if (summary_type[i] %in% c(1L, 2L, 3L)) {
      central_estimates[i] <- obs_stat1[i]
    } else if (summary_type[i] == 4L && freq_len_vec[i] > 0) {
      s  <- freq_start_vec[i]
      ln <- freq_len_vec[i]
      fv <- freq_value_all[s:(s + ln - 1)]
      fc <- freq_count_all[s:(s + ln - 1)]
      central_estimates[i] <- sum(fv * fc) / sum(fc)
    } else if (summary_type[i] == 5L && freq_len_vec[i] > 0) {
      s   <- freq_start_vec[i]
      ln  <- freq_len_vec[i]
      fl  <- freq_lower_all[s:(s + ln - 1)]
      fu  <- freq_upper_all[s:(s + ln - 1)]
      fc  <- freq_count_all[s:(s + ln - 1)]
      mid <- (fl + fu) / 2
      central_estimates[i] <- sum(mid * fc) / sum(fc)
    }
  }
  valid_centrals <- central_estimates[central_estimates > 0]

  # Create base Stan data
  stan_data <- list(
    n_datasets   = n_datasets,
    n_obs        = as.array(n_obs_vec),
    summary_type = as.array(summary_type),
    dist_type    = dist_type,
    obs_stat1    = as.array(obs_stat1),
    obs_stat2    = as.array(obs_stat2),
    obs_stat3    = as.array(obs_stat3),
    n_freq_total = length(freq_value_all),
    freq_value   = freq_value_all,
    freq_lower   = freq_lower_all,
    freq_upper   = freq_upper_all,
    freq_count   = freq_count_all,
    freq_start   = as.array(freq_start_vec),
    freq_len     = as.array(freq_len_vec)
  )

  stan_data$mu0_mean     <- if (length(valid_centrals) > 0) log(mean(valid_centrals)) else 0
  stan_data$mu0_sd       <- custom_priors$mu0_sd
  stan_data$log_tau_mean <- custom_priors$log_tau_mean
  stan_data$log_tau_sd   <- custom_priors$log_tau_sd
  stan_data$log_phi_mean <- custom_priors$log_phi_mean
  stan_data$log_phi_sd   <- custom_priors$log_phi_sd

  return(stan_data)
}



# Pre-inference checks -----------------------------------------------------

#' Run pre-inference checks on a list of datasets
#'
#' Performs a suite of fast, pre-MCMC checks to detect data issues that are
#' likely to cause convergence problems. Checks are run in order of increasing
#' computational cost and a summary is printed to the console.
#'
#' The five checks performed are:
#' \describe{
#'   \item{1. Method-of-moments consistency}{Estimates `phi` from each dataset
#'     individually using moment-based approximations and flags any dataset
#'     whose implied `phi` is more than `phi_outlier_threshold` times the
#'     median of all implied values.}
#'   \item{2. Prior predictive compatibility}{Simulates summary statistics from
#'     the prior and checks whether each observed value falls within the 95\%
#'     prior predictive interval. Datasets outside this range suggest a
#'     prior--data mismatch.}
#'   \item{3. MAP optimisation probe}{Runs [rstan::optimizing()] as a fast
#'     proxy for MCMC convergence. Failure or extreme `phi` at the MAP
#'     estimate is a reliable early warning that HMC will struggle.}
#'   \item{4. Log-likelihood surface scan}{Evaluates the joint log-posterior
#'     over a grid of `phi` values (other parameters held at the MAP). A
#'     multimodal or sharply peaked surface explains treedepth exhaustion.}
#'   \item{5. Leave-one-out single-dataset fits}{Fits the model to each dataset
#'     individually and compares the resulting `phi` posteriors. Non-overlapping
#'     credible intervals identify the specific datasets driving tension.}
#' }
#'
#' @param datasets A named list of datasets in the format accepted by
#'   [prepare_stan_data_from_datasets()].
#' @param stan_model A compiled Stan model object from [rstan::stan_model()].
#' @param dist_type Integer distribution code: `1` = log-normal, `2` = gamma,
#'   `3` = Weibull. Defaults to `1`.
#' @param custom_priors Optional named list of prior overrides passed to
#'   [prepare_stan_data_from_datasets()].
#' @param phi_outlier_threshold Multiplier used in the method-of-moments check.
#'   A dataset is flagged if its implied `phi` exceeds
#'   `phi_outlier_threshold * median(implied_phi)`. Defaults to `5`.
#' @param phi_grid Numeric vector of `phi` values for the log-likelihood
#'   surface scan. Defaults to `seq(0.5, 50, by = 0.5)`.
#' @param n_sim Number of draws for the prior predictive check. Defaults to
#'   `2000`.
#' @param loo_iter Number of MCMC iterations per chain for the leave-one-out
#'   single-dataset fits. Defaults to `4000`.
#' @param loo_chains Number of chains for the leave-one-out fits. Defaults to
#'   `2`.
#' @param verbose Logical. If `TRUE` (default), prints a formatted summary of
#'   all check results to the console.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`mom_consistency`}{Data frame of implied `phi` per dataset with
#'       an `is_outlier` flag.}
#'     \item{`prior_predictive`}{Data frame of 95\% prior predictive intervals
#'       for the implied SD of each dataset, with an `outside_prior_pi` flag.}
#'     \item{`map_probe`}{List with `phi_map` (MAP estimate of `phi`) and
#'       `map_converged` logical.}
#'     \item{`ll_surface`}{Data frame of `phi` vs `log_prob` from the surface
#'       scan.}
#'     \item{`loo_fits`}{Data frame of per-dataset `phi` posterior summaries
#'       from the leave-one-out fits.}
#'   }
#' @export
pre_inference_checks <- function(datasets,
                                 stan_model,
                                 dist_type             = 1,
                                 custom_priors         = list(),
                                 phi_outlier_threshold = 5,
                                 phi_grid              = seq(0.5, 50, by = 0.5),
                                 n_sim                 = 2000,
                                 loo_iter              = 4000,
                                 loo_chains            = 2,
                                 verbose               = TRUE) {

  dist_name <- c("1" = "lognormal", "2" = "gamma", "3" = "weibull")[[as.character(dist_type)]]
  stan_data <- prepare_stan_data_from_datasets(datasets, dist_type = dist_type,
                                               custom_priors = custom_priors)

  results <- list()

  # ── Check 1: Method-of-moments consistency ────────────────────────────────
  # Estimate implied phi from each dataset individually. For each distribution:
  #   lognormal : phi = log-SD -> approximate as sd(log(x)); use CV approximation
  #   gamma     : phi = shape  -> (mean/sd)^2
  #   weibull   : phi = shape  -> approximate from CV via Newton iteration
  mom_df <- purrr::imap_dfr(datasets, function(d, name) {
    mean_est <- sd_est <- NA_real_

    if (!is.null(d$mean) && !is.null(d$sd)) {
      mean_est <- d$mean
      sd_est   <- d$sd
    } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) {
      mean_est <- d$median
      sd_est   <- (d$Q3 - d$Q1) / 1.35
    } else if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) {
      mean_est <- d$median
      sd_est   <- (d$max - d$min) / 4
    } else if (!is.null(d$freq_value) && !is.null(d$freq_count)) {
      w        <- d$freq_count / sum(d$freq_count)
      mean_est <- sum(d$freq_value * w)
      sd_est   <- sqrt(sum(w * (d$freq_value - mean_est)^2))
    } else if (!is.null(d$freq_lower) && !is.null(d$freq_upper) && !is.null(d$freq_count)) {
      mid      <- (d$freq_lower + d$freq_upper) / 2
      w        <- d$freq_count / sum(d$freq_count)
      mean_est <- sum(mid * w)
      sd_est   <- sqrt(sum(w * (mid - mean_est)^2))
    }

    implied_phi <- NA_real_
    if (!is.na(mean_est) && !is.na(sd_est) && sd_est > 0) {
      implied_phi <- switch(dist_name,
        lognormal = log(1 + (sd_est / mean_est)^2),        # approx log-variance
        gamma     = (mean_est / sd_est)^2,                 # shape = (mean/sd)^2
        weibull   = {                                       # invert CV numerically
          cv <- sd_est / mean_est
          # CV^2 = Gamma(1+2/k)/Gamma(1+1/k)^2 - 1; solve for k
          obj <- function(k) sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1) - cv
          tryCatch(stats::uniroot(obj, c(0.1, 200))$root, error = function(e) NA_real_)
        }
      )
    }

    tibble::tibble(dataset = name, mean_est = mean_est, sd_est = sd_est,
                   implied_phi = implied_phi)
  })

  med_phi  <- stats::median(mom_df$implied_phi, na.rm = TRUE)
  mom_df   <- dplyr::mutate(
    mom_df,
    is_outlier = !is.na(implied_phi) &
      (implied_phi > phi_outlier_threshold * med_phi |
       implied_phi < med_phi / phi_outlier_threshold)
  )
  results$mom_consistency <- mom_df

  # ── Check 2: Prior predictive compatibility ────────────────────────────────
  mu0_draws  <- stats::rnorm(n_sim, stan_data$mu0_mean, stan_data$mu0_sd)
  tau_draws  <- exp(stats::rnorm(n_sim, stan_data$log_tau_mean, stan_data$log_tau_sd))
  phi_draws  <- exp(stats::rnorm(n_sim, stan_data$log_phi_mean, stan_data$log_phi_sd))
  loc_draws  <- stats::rnorm(n_sim, mu0_draws, tau_draws)

  sim_sd <- switch(dist_name,
    lognormal = exp(loc_draws) * sqrt(exp(phi_draws^2) - 1),
    gamma     = exp(loc_draws) / sqrt(phi_draws),
    weibull   = exp(loc_draws) * sqrt(gamma(1 + 2/phi_draws) - gamma(1 + 1/phi_draws)^2)
  )

  prior_pi <- stats::quantile(sim_sd, c(0.025, 0.975), na.rm = TRUE)

  prior_df <- dplyr::mutate(
    mom_df,
    prior_sd_lo        = prior_pi[[1]],
    prior_sd_hi        = prior_pi[[2]],
    outside_prior_pi   = !is.na(sd_est) &
      (sd_est < prior_pi[[1]] | sd_est > prior_pi[[2]])
  )
  results$prior_predictive <- prior_df

  # ── Check 3: MAP optimisation probe ───────────────────────────────────────
  map_result <- tryCatch({
    opt        <- rstan::optimizing(stan_model, data = stan_data, hessian = FALSE,
                                    refresh = 0)
    phi_map    <- exp(opt$par[["log_phi"]])
    list(phi_map = phi_map, map_converged = opt$return_code == 0,
         return_code = opt$return_code)
  }, error = function(e) {
    list(phi_map = NA_real_, map_converged = FALSE, return_code = NA_integer_,
         error_msg = conditionMessage(e))
  })
  results$map_probe <- map_result

  # ── Check 4: Log-likelihood surface scan ──────────────────────────────────
  mu0_init      <- stan_data$mu0_mean
  log_tau_init  <- stan_data$log_tau_mean
  loc_d_raw_init <- rep(0, stan_data$n_datasets)

  ll_surface <- purrr::map_dfr(phi_grid, function(phi_val) {
    pars <- list(mu0       = mu0_init,
                 log_tau   = log_tau_init,
                 log_phi   = log(phi_val),
                 loc_d_raw = loc_d_raw_init)
    lp <- tryCatch({
      rstan::log_prob(stan_model,
                      rstan::unconstrain_pars(stan_model, data = stan_data, pars = pars),
                      adjust_transform = TRUE)
    }, error = function(e) NA_real_)
    tibble::tibble(phi = phi_val, log_prob = lp)
  })
  results$ll_surface <- ll_surface

  # ── Check 5: Leave-one-out single-dataset fits ────────────────────────────
  loo_fits <- purrr::imap_dfr(datasets, function(d, name) {
    # suppressWarnings() is intentional here: n_datasets = 1 is expected for
    # each LOO fit (triggering the tau identifiability warning), and divergent
    # transitions / treedepth warnings from individual fits are uninformative
    # in this diagnostic context.
    single_data <- tryCatch(
      suppressWarnings(
        prepare_stan_data_from_datasets(stats::setNames(list(d), name),
                                        dist_type     = dist_type,
                                        custom_priors = custom_priors)
      ),
      error = function(e) NULL
    )
    if (is.null(single_data)) {
      return(tibble::tibble(dataset = name, phi_mean = NA_real_,
                            phi_lo = NA_real_, phi_hi = NA_real_,
                            rhat = NA_real_,   n_eff = NA_real_))
    }

    fit <- tryCatch(
      suppressWarnings(
        rstan::sampling(stan_model, data = single_data,
                        iter = loo_iter, chains = loo_chains,
                        refresh = 0, show_messages = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble::tibble(dataset = name, phi_mean = NA_real_,
                            phi_lo = NA_real_, phi_hi = NA_real_,
                            rhat = NA_real_,   n_eff = NA_real_))
    }

    s <- rstan::summary(fit, pars = "phi")$summary
    tibble::tibble(dataset  = name,
                   phi_mean = s[, "mean"],
                   phi_lo   = s[, "2.5%"],
                   phi_hi   = s[, "97.5%"],
                   rhat     = s[, "Rhat"],
                   n_eff    = s[, "n_eff"])
  })
  results$loo_fits <- loo_fits

  # ── Verbose summary ───────────────────────────────────────────────────────
  if (verbose) {
    cli::cli_h1("Pre-inference checks ({dist_name}, {length(datasets)} datasets)")

    cli::cli_h2("1. Method-of-moments consistency (phi outlier threshold: {phi_outlier_threshold}x median)")
    print(dplyr::select(mom_df, dataset, implied_phi, is_outlier))
    n_out <- sum(mom_df$is_outlier, na.rm = TRUE)
    if (n_out > 0) {
      cli::cli_alert_warning("{n_out} dataset(s) have an implied phi far from the others: {mom_df$dataset[mom_df$is_outlier]}")
    } else {
      cli::cli_alert_success("All implied phi values are broadly consistent")
    }

    cli::cli_h2("2. Prior predictive compatibility (95% PI for SD: [{round(prior_pi[[1]], 2)}, {round(prior_pi[[2]], 2)}])")
    outside <- dplyr::filter(prior_df, outside_prior_pi)
    if (nrow(outside) > 0) {
      cli::cli_alert_warning("{nrow(outside)} dataset(s) have SD outside the 95% prior predictive interval: {outside$dataset}")
    } else {
      cli::cli_alert_success("All observed SDs are within the 95% prior predictive interval")
    }

    cli::cli_h2("3. MAP optimisation probe")
    if (!map_result$map_converged) {
      cli::cli_alert_danger("MAP optimisation failed (return code {map_result$return_code}) — MCMC likely to struggle")
    } else {
      cli::cli_alert_success("MAP converged; phi_MAP = {round(map_result$phi_map, 2)}")
      if (!is.na(map_result$phi_map) && (map_result$phi_map > 50 || map_result$phi_map < 0.1)) {
        cli::cli_alert_warning("phi_MAP = {round(map_result$phi_map, 2)} is extreme — check prior and data consistency")
      }
    }

    cli::cli_h2("4. Log-likelihood surface scan")
    finite_ll <- dplyr::filter(ll_surface, is.finite(log_prob))
    if (nrow(finite_ll) > 0) {
      peak_phi <- finite_ll$phi[which.max(finite_ll$log_prob)]
      cli::cli_alert_info("Surface peak at phi ~ {peak_phi}. Plot with: plot(results$ll_surface$phi, results$ll_surface$log_prob, type = 'l')")
    } else {
      cli::cli_alert_danger("Log-likelihood surface is entirely non-finite — severe model/data mismatch")
    }

    cli::cli_h2("5. Leave-one-out single-dataset phi posteriors")
    print(loo_fits)
    phi_ranges_overlap <- function(df) {
      # Flag any dataset whose 95% CI does not overlap the majority
      med_lo <- stats::median(df$phi_lo, na.rm = TRUE)
      med_hi <- stats::median(df$phi_hi, na.rm = TRUE)
      dplyr::mutate(df,
        no_overlap = !is.na(phi_lo) & (phi_hi < med_lo | phi_lo > med_hi))
    }
    loo_flagged <- phi_ranges_overlap(loo_fits)
    n_no_overlap <- sum(loo_flagged$no_overlap, na.rm = TRUE)
    if (n_no_overlap > 0) {
      cli::cli_alert_warning(
        "{n_no_overlap} dataset(s) have phi posteriors that do not overlap the majority: {loo_flagged$dataset[loo_flagged$no_overlap]}"
      )
    } else {
      cli::cli_alert_success("All per-dataset phi posteriors broadly overlap")
    }
  }

  invisible(results)
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
                            summary_config = c("fixed", "mixed_balanced", "mixed_random",
                                               "mixed_balanced_with_freq", "mixed_random_with_freq",
                                               "custom"),
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

  } else if (summary_config == "mixed_balanced_with_freq") {
    # Equal representation of all four summary types (1-4)
    summary_type <- rep(1:4, length.out = n_datasets)
    summary_type <- sample(summary_type)  # Shuffle

  } else if (summary_config == "mixed_random_with_freq") {
    # Random mix including frequency table type
    summary_type <- sample(1:4, n_datasets, replace = TRUE,
                           prob = c(0.4, 0.15, 0.25, 0.2))

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
                                      include_varied_n = TRUE,
                                      include_freq_table = FALSE) {
  
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
  
  # ===== HOMOGENEOUS FREQ TABLE SCENARIOS =====
  if (include_homogeneous && include_freq_table) {
    for (dist in distributions) {
      for (n_datasets_val in c(1,3)) {
        for (n_obs_val in c(10, 40)) {
          scenarios[[idx]] <- create_scenario(
            scenario_name      = sprintf("Homog_%s_ST4_D%d_N%d", dist, n_datasets_val, n_obs_val),
            dist_type          = dist,
            n_datasets         = n_datasets_val,
            n_obs_config       = "fixed",
            fixed_n_obs        = n_obs_val,
            summary_config     = "fixed",
            fixed_summary_type = 4
          )
          idx <- idx + 1
        }
      }
    }
  }
  
  # ===== MIXED SCENARIOS WITH FREQUENCY TABLE =====
  if (include_mixed && include_freq_table) {
    for (dist in distributions) {
      for (n_datasets_val in c(15, 30)) {
        scenarios[[idx]] <- create_scenario(
          scenario_name  = sprintf("Mixed_Balanced_Freq_%s_D%d", dist, n_datasets_val),
          dist_type      = dist,
          n_datasets     = n_datasets_val,
          n_obs_config   = "fixed",
          fixed_n_obs    = 30,
          summary_config = "mixed_balanced_with_freq"
        )
        idx <- idx + 1

        scenarios[[idx]] <- create_scenario(
          scenario_name  = sprintf("Mixed_Random_Freq_%s_D%d", dist, n_datasets_val),
          dist_type      = dist,
          n_datasets     = n_datasets_val,
          n_obs_config   = "fixed",
          fixed_n_obs    = 30,
          summary_config = "mixed_random_with_freq"
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
      summary_type_4_prop = mean(s$summary_type == 4),
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
  mu0_samples  <- draws$mu0
  tau_samples  <- exp(draws$log_tau)
  phi_samples  <- exp(draws$log_phi)
  loc_d_samples <- draws$loc_d          # [n_samples x n_datasets]

  n_samples  <- length(mu0_samples)
  n_datasets <- dim(loc_d_samples)[2]
  
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

      if (n_datasets < 5) {
        # When n_datasets < 5, tau is not identified and mu0 is confounded
        # with tau * loc_d_raw. Use mean(loc_d) as a point estimate of the
        # study-level location, consistent with the Stan generated quantities
        # block and compute_predictive_cdf().
        loc_point <- mean(loc_d_samples[idx, ])
        phi_s     <- phi_samples[idx]

        density_samples[s] <- if (dist_type == "lognormal") {
          dlnorm(x, meanlog = loc_point, sdlog = phi_s)
        } else if (dist_type == "gamma") {
          mean_d <- exp(loc_point)
          dgamma(x, shape = phi_s, rate = phi_s / mean_d)
        } else if (dist_type == "weibull") {
          dweibull(x, shape = phi_s, scale = exp(loc_point))
        }

      } else {
        # n_datasets >= 5: integrate over Normal(mu0, tau) random effects
        integrand <- function(loc) {
          if (dist_type == "lognormal") {
            dlnorm(x, meanlog = loc, sdlog = phi_samples[idx]) *
              dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
          } else if (dist_type == "gamma") {
            mean_d <- exp(loc)
            shape  <- phi_samples[idx]
            dgamma(x, shape = shape, rate = shape / mean_d) *
              dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
          } else if (dist_type == "weibull") {
            scale <- exp(loc)
            shape <- phi_samples[idx]
            dweibull(x, shape = shape, scale = scale) *
              dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
          }
        }

        density_samples[s] <- integrate(integrand,
                                        lower = mu0_samples[idx] - 5 * tau_samples[idx],
                                        upper = mu0_samples[idx] + 5 * tau_samples[idx])$value
      }
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
#' @importFrom foreach %dopar%
#' @param n_sim Number of simulation replicates per scenario
#' @param scenarios_df Data frame from generate_scenario_library()
#' @param stan_model Compiled Stan model
#' @param seed Random seed
#' @param save_name File path to save final results as RDS.
#' @param n_cores Number of parallel workers requested. Automatically capped at
#'   `floor(detectCores() / chains)` to avoid CPU oversubscription.
#' @param chains Number of Stan chains per fit (default 4). Used to compute the
#'   worker cap.
#' @return Data frame with one row per simulation replicate and columns for
#'   scenario metadata, true parameter values, coverage, bias, and IQD metrics.
#' @export
run_simulation_study_generalized <- function(n_sim,
                                             scenarios_df,
                                             stan_model,
                                             seed = 123,
                                             save_name = "simulation_results_tmp_general.rds",
                                             n_cores = 1,
                                             chains = 4) {

  full_scenarios <- attr(scenarios_df, "full_scenarios")
  if (is.null(full_scenarios)) {
    stop("scenarios_df must be created by generate_scenario_library()")
  }

  # Cap workers to avoid oversubscription: each worker runs `chains` Stan chains
  max_workers <- max(1L, floor(parallel::detectCores() / chains))
  n_workers   <- min(n_cores, max_workers)

  if (n_cores > max_workers) {
    message("Capping n_cores from ", n_cores, " to ", max_workers,
            " (floor(", parallel::detectCores(), " cores / ", chains, " chains))")
  }

  # Flat task list: one row per (scenario, sim) pair
  tasks <- expand.grid(scenario_idx = seq_len(nrow(scenarios_df)),
                       sim          = seq_len(n_sim))

  # Register parallel or sequential backend
  if (n_workers > 1) {
    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    message("Running ", nrow(tasks), " tasks on ", n_workers, " workers ",
            "(", chains, " chains each).")
  } else {
    foreach::registerDoSEQ()
  }

  run_one <- function(scenario_idx, sim, scenario, stan_model, seed, chains, n_sim) {
    # Reproducible per-task seed
    set.seed(seed + scenario_idx * n_sim + sim)

    if (startsWith(scenario$scenario_name, "VarN")) {
      n_obs_in <- pmax(
        pmin(round(stats::rnorm(scenario$n_datasets, scenario$n_obs_mean, scenario$n_obs_sd)),
             scenario$n_obs_max),
        scenario$n_obs_min)
    } else {
      n_obs_in <- scenario$n_obs_mean
    }

    st4_prop <- if (!is.null(scenario$summary_type_4_prop)) scenario$summary_type_4_prop else 0
    summary_type_arg <- if (scenario$summary_type_1_prop == 1) {
      1L
    } else if (scenario$summary_type_2_prop == 1) {
      2L
    } else if (scenario$summary_type_3_prop == 1) {
      3L
    } else if (st4_prop == 1) {
      4L
    } else {
      c(scenario$summary_type_1_prop,
        scenario$summary_type_2_prop,
        scenario$summary_type_3_prop,
        st4_prop)
    }

    sim_data <- generate_hierarchical_data_mixed(
      n_datasets   = scenario$n_datasets,
      n_obs        = n_obs_in,
      dist_type    = scenario$dist_type,
      mu0          = scenario$mu0,
      tau          = scenario$tau,
      phi          = scenario$phi,
      summary_type = summary_type_arg
    )

    # Number of distinct summary types present in this scenario
    summary_type_diversity <- sum(
      c(scenario$summary_type_1_prop,
        scenario$summary_type_2_prop,
        scenario$summary_type_3_prop,
        st4_prop) > 0)

    tryCatch({
      fit <- fit_model(sim_data, stan_model,
                       chains  = chains,
                       refresh = 0,
                       control = list(adapt_delta = 0.95, max_treedepth = 12))

      fit_summary <- summary(fit)$summary
      max_rhat    <- max(fit_summary[, "Rhat"],  na.rm = TRUE)
      min_neff    <- min(fit_summary[, "n_eff"], na.rm = TRUE)

      if (max_rhat > 1.1 || min_neff < 100) {
        warning(paste("Convergence issues in scenario", scenario_idx,
                      "sim", sim, ": Rhat =", round(max_rhat, 3),
                      ", min n_eff =", round(min_neff, 0)))
      }

      coverage_mu0 <- check_coverage(fit, "mu0", sim_data$true_params$mu0)
      coverage_tau <- check_coverage(fit, "tau", sim_data$true_params$tau)
      coverage_phi <- check_coverage(fit, "phi", sim_data$true_params$phi)

      bias_mu0 <- compute_median_bias(fit, "mu0", sim_data$true_params$mu0)
      bias_tau <- compute_median_bias(fit, "tau", sim_data$true_params$tau)
      bias_phi <- compute_median_bias(fit, "phi", sim_data$true_params$phi)

      iqd <- compute_iqd(fit, sim_data$true_params, scenario$dist_type)

      data.frame(
        scenario_idx           = scenario_idx,
        scenario_name          = scenario$scenario_name,
        sim                    = sim,
        dist_type              = scenario$dist_type,
        n_datasets             = scenario$n_datasets,
        n_obs_mean             = scenario$n_obs_mean,
        n_obs_sd               = scenario$n_obs_sd,
        n_obs_min              = scenario$n_obs_min,
        n_obs_max              = scenario$n_obs_max,
        summary_type_diversity = summary_type_diversity,
        prop_summary_type_1    = scenario$summary_type_1_prop,
        prop_summary_type_2    = scenario$summary_type_2_prop,
        prop_summary_type_3    = scenario$summary_type_3_prop,
        prop_summary_type_4    = st4_prop,
        true_mu0               = scenario$mu0,
        true_tau               = scenario$tau,
        true_phi               = scenario$phi,
        coverage_mu0           = coverage_mu0,
        coverage_tau           = coverage_tau,
        coverage_phi           = coverage_phi,
        bias_mu0               = bias_mu0,
        bias_tau               = bias_tau,
        bias_phi               = bias_phi,
        rel_bias_mu0           = bias_mu0 / sim_data$true_params$mu0,
        rel_bias_tau           = bias_tau / sim_data$true_params$tau,
        rel_bias_phi           = bias_phi / sim_data$true_params$phi,
        iqd                    = iqd,
        max_rhat               = max_rhat,
        min_neff               = min_neff,
        converged              = max_rhat <= 1.1 & min_neff >= 100,
        stringsAsFactors       = FALSE
      )
    }, error = function(e) {
      warning(paste("Error in scenario", scenario_idx, "sim", sim, ":", e$message))
      data.frame(
        scenario_idx           = scenario_idx,
        scenario_name          = scenario$scenario_name,
        sim                    = sim,
        dist_type              = scenario$dist_type,
        n_datasets             = scenario$n_datasets,
        n_obs_mean             = scenario$n_obs_mean,
        n_obs_sd               = scenario$n_obs_sd,
        n_obs_min              = scenario$n_obs_min,
        n_obs_max              = scenario$n_obs_max,
        summary_type_diversity = summary_type_diversity,
        prop_summary_type_1    = scenario$summary_type_1_prop,
        prop_summary_type_2    = scenario$summary_type_2_prop,
        prop_summary_type_3    = scenario$summary_type_3_prop,
        prop_summary_type_4    = st4_prop,
        true_mu0               = scenario$mu0,
        true_tau               = scenario$tau,
        true_phi               = scenario$phi,
        coverage_mu0           = NA,
        coverage_tau           = NA,
        coverage_phi           = NA,
        bias_mu0               = NA,
        bias_tau               = NA,
        bias_phi               = NA,
        rel_bias_mu0           = NA,
        rel_bias_tau           = NA,
        rel_bias_phi           = NA,
        iqd                    = NA,
        max_rhat               = NA,
        min_neff               = NA,
        converged              = FALSE,
        stringsAsFactors       = FALSE
      )
    })
  }

  results <- foreach::foreach(
    i          = seq_len(nrow(tasks)),
    .combine   = dplyr::bind_rows,
    .packages  = c("rstan", "ddsynth"),
    .export    = c("stan_model", "run_one")
  ) %dopar% {
    run_one(
      scenario_idx = tasks$scenario_idx[i],
      sim          = tasks$sim[i],
      scenario     = scenarios_df[tasks$scenario_idx[i], ],
      stan_model   = stan_model,
      seed         = seed,
      chains       = chains,
      n_sim        = n_sim
    )
  }

  saveRDS(results, save_name)
  results
}

run_simulation_study_generalized_non_parallel <- function(n_sim, 
                                             scenarios_df, 
                                             stan_model, 
                                             seed = 123,
                                             save_name = "simulation_results_tmp_general.rds") {
  
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
      st4_prop_np <- if (!is.null(scenario$summary_type_4_prop)) scenario$summary_type_4_prop else 0
      summary_type_np <- if (scenario$summary_type_1_prop == 1) {
        1L
      } else if (scenario$summary_type_2_prop == 1) {
        2L
      } else if (scenario$summary_type_3_prop == 1) {
        3L
      } else if (st4_prop_np == 1) {
        4L
      } else {
        c(scenario$summary_type_1_prop,
          scenario$summary_type_2_prop,
          scenario$summary_type_3_prop,
          st4_prop_np)
      }

      sim_data <- generate_hierarchical_data_mixed(
        n_datasets = scenario$n_datasets,
        n_obs = n_obs_in,
        dist_type = scenario$dist_type,
        mu0 = scenario$mu0,
        tau = scenario$tau,
        phi = scenario$phi,
        summary_type = summary_type_np
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
          prop_summary_type_4 = mean(scenario$summary_type == 4),
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
          prop_summary_type_4 = mean(scenario$summary_type == 4),
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

#' Generate data from hierarchical model with mixed summary types and sample sizes
#'
#' @param n_datasets Number of datasets to generate
#' @param n_obs Vector of sample sizes for each dataset (can vary)
#' @param dist_type Distribution type: "lognormal", "gamma", or "weibull"
#' @param mu0 Population mean (location parameter)
#' @param tau Between-study standard deviation
#' @param phi Distribution-specific shape/scale parameter
#' @param summary_type Summary type specification. Can be: `NULL` (random mix of
#'   types 1-3), a single integer 1-4 (all datasets use that type), a vector of
#'   length `n_datasets` (one type per dataset), a length-3 probability vector
#'   (sample from types 1-3 with those probabilities), or a length-4 probability
#'   vector (sample from types 1-4 with those probabilities). Type 4 produces a
#'   frequency table of observations rounded to the nearest whole day.
#' @return List containing true parameters and observed summary statistics
generate_hierarchical_data_mixed <- function(n_datasets,
                                             n_obs,
                                             dist_type = c("lognormal", "gamma", "weibull"),
                                             mu0,
                                             tau,
                                             phi,
                                             summary_type = NULL) {
  
  dist_type <- match.arg(dist_type)
  
  # If n_obs is a single value, replicate it
  if (length(n_obs) == 1) {
    n_obs <- rep(n_obs, n_datasets)
  }
  
  # Validate n_obs length
  if (length(n_obs) != n_datasets) {
    stop("Length of n_obs must equal n_datasets or be a single value")
  }
  
  # If summary_type is NULL or single value, handle appropriately
  if (is.null(summary_type)) {
    # Default: random mix of types 1-3 (no freq table)
    summary_type <- sample(1:3, n_datasets, replace = TRUE)
  } else if (length(summary_type) == 4) {
    summary_type <- sample(1:4, n_datasets, replace = TRUE, prob = summary_type)
  } else if (length(summary_type) == 3) {
    summary_type <- sample(1:3, n_datasets, replace = TRUE, prob = summary_type)
  } else if (length(summary_type) == 1) {
    summary_type <- rep(summary_type, n_datasets)
  }
  
  # Validate summary_type length
  if (length(summary_type) != n_datasets) {
    stop("Length of summary_type must equal n_datasets or be a single value")
  }
  
  # Generate study-specific location parameters
  loc_d <- rnorm(n_datasets, mean = mu0, sd = tau)
  
  # Initialize storage
  obs_stat1   <- numeric(n_datasets)
  obs_stat2   <- numeric(n_datasets)
  obs_stat3   <- numeric(n_datasets)
  freq_tables <- vector("list", n_datasets)
  
  # Generate data for each dataset
  for (d in 1:n_datasets) {
    n <- n_obs[d]
    loc <- loc_d[d]
    st <- summary_type[d]
    
    # Generate raw data based on distribution type
    if (dist_type == "lognormal") {
      data_d <- rlnorm(n, meanlog = loc, sdlog = phi)
      
    } else if (dist_type == "gamma") {
      mean_d <- exp(loc)
      shape <- phi
      rate <- shape / mean_d
      data_d <- rgamma(n, shape = shape, rate = rate)
      
    } else if (dist_type == "weibull") {
      scale <- exp(loc)
      shape <- phi
      data_d <- rweibull(n, shape = shape, scale = scale)
    }
    
    # Compute summary statistics based on type for this specific dataset
    if (st == 1) {  # median + range
      obs_stat1[d] <- median(data_d)
      obs_stat2[d] <- min(data_d)
      obs_stat3[d] <- max(data_d)
      
    } else if (st == 2) {  # median + IQR
      obs_stat1[d] <- median(data_d)
      obs_stat2[d] <- quantile(data_d, 0.25)
      obs_stat3[d] <- quantile(data_d, 0.75)
      
    } else if (st == 3) {  # mean + sd
      obs_stat1[d] <- mean(data_d)
      obs_stat2[d] <- sd(data_d)
      obs_stat3[d] <- 0  # placeholder

    } else if (st == 4) {  # frequency table (rounded to full days)
      data_d_rounded <- pmax(round(data_d), 1L)
      freq_tbl       <- table(data_d_rounded)
      freq_tables[[d]] <- list(
        value = as.numeric(names(freq_tbl)),
        count = as.integer(freq_tbl)
      )
      obs_stat1[d] <- 0  # placeholder
      obs_stat2[d] <- 0  # placeholder
      obs_stat3[d] <- 0  # placeholder
    }
  }

  # Build flat frequency-table arrays required by the Stan model
  freq_value_flat <- numeric(0)
  freq_count_flat <- integer(0)
  freq_start_vec  <- integer(n_datasets)
  freq_len_vec    <- integer(n_datasets)
  running_start   <- 1L

  for (d in seq_len(n_datasets)) {
    if (summary_type[d] == 4 && !is.null(freq_tables[[d]])) {
      ft                <- freq_tables[[d]]
      freq_start_vec[d] <- running_start
      freq_len_vec[d]   <- length(ft$value)
      freq_value_flat   <- c(freq_value_flat, ft$value)
      freq_count_flat   <- c(freq_count_flat, ft$count)
      running_start     <- running_start + freq_len_vec[d]
    }
  }

  list(
    true_params = list(
      mu0 = mu0,
      tau = tau,
      phi = phi,
      loc_d = loc_d
    ),
    obs_data = list(
      n_datasets   = n_datasets,
      n_obs        = as.array(n_obs),
      summary_type = as.array(summary_type),
      dist_type    = switch(dist_type,
                            "lognormal" = 1,
                            "gamma"     = 2,
                            "weibull"   = 3),
      obs_stat1    = as.array(obs_stat1),
      obs_stat2    = as.array(obs_stat2),
      obs_stat3    = as.array(obs_stat3),
      # Frequency table fields (populated only when summary_type == 4)
      n_freq_total = length(freq_value_flat),
      freq_value   = freq_value_flat,
      freq_count   = freq_count_flat,
      freq_start   = as.array(freq_start_vec),
      freq_len     = as.array(freq_len_vec),
      # Default priors
      mu0_mean     = 1,
      mu0_sd       = 2,
      log_tau_mean = 0.2,
      log_tau_sd   = 0.5,
      log_phi_mean = ifelse(dist_type == "lognormal", 0.2,
                            ifelse(dist_type == "gamma", 1.0, 1.0)),
      log_phi_sd   = 1
    )
  )
}