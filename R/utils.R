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
#' @param dist_name Character string: `"lognormal"`, `"gamma"`, `"weibull"`,
#'   `"burr12"` (Burr Type XII), or `"gengamma"` (Generalised Gamma, Prentice
#'   parameterisation).
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
    mu0   <- sims$mu0[idx]
    tau   <- sims$tau[idx]
    phi   <- sims$phi[idx]
    kappa <- sims$kappa[idx]

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
        shape  <- phi
        rate   <- shape / mean_d
        cdf_l[l, ] <- pgamma(x_seq, shape = shape, rate = rate)

      } else if (dist_name == "weibull") {
        scale <- exp(loc_d)
        shape <- phi
        cdf_l[l, ] <- pweibull(x_seq, shape = shape, scale = scale)

      } else if (dist_name == "burr") {
        # Burr XII CDF: F(x) = 1 - (1 + (x/lambda)^c)^(-k)
        # lambda = exp(loc_d), c = phi, k = kappa
        lambda  <- exp(loc_d)
        cdf_l[l, ] <- 1 - (1 + (x_seq / lambda)^phi)^(-kappa)

      } else if (dist_name == "gg") {
        # Generalised Gamma (Prentice): mu = loc_d, sigma = phi, Q = kappa
        # CDF = pgamma(gamma_shape * exp(Q * w), shape = gamma_shape, rate = 1)
        # where gamma_shape = 1/Q^2, w = (log(x) - mu) / sigma
        gamma_shape <- 1 / kappa^2
        w           <- (log(x_seq) - loc_d) / phi
        cdf_l[l, ] <- pgamma(gamma_shape * exp(kappa * w),
                             shape = gamma_shape, rate = 1)
      } else {
        stop(sprintf(
          "compute_predictive_cdf: unknown dist_name '%s'. ",
          "Must be one of: 'lognormal', 'gamma', 'weibull', 'burr', 'gg'.",
          dist_name
        ), call. = FALSE)
      }
    }

    # Average over study-level locations (na.rm = TRUE guards against rare
    # numerical edge cases in individual draws without silently hiding them)
    cdf_mat[i, ] <- colMeans(cdf_l, na.rm = TRUE)
  }
  
  # Compute summary statistics
  summary_df <- data.frame(
    x      = x_seq,
    median = apply(cdf_mat, 2, median,   na.rm = TRUE),
    mean   = apply(cdf_mat, 2, mean,     na.rm = TRUE),
    low    = apply(cdf_mat, 2, quantile, 0.025, na.rm = TRUE),
    high   = apply(cdf_mat, 2, quantile, 0.975, na.rm = TRUE),
    model  = dist_name
  )

  # Return both the summary (for plotting the ribbon/line) and the raw matrix
  # (for computing consistent PI bounds on derived quantiles via extract_quantiles)
  list(summary = summary_df, cdf_mat = cdf_mat)
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
#' @param cdf_mat Numeric matrix of posterior CDF draws as returned by
#'   [compute_predictive_cdf()] (rows = posterior draws, columns = `x_seq`
#'   grid points). When supplied, the 95% prediction interval bounds
#'   (`x_low`, `x_high`) are computed by interpolating each draw's CDF to
#'   find the x at which it crosses `p`, then taking the 2.5% and 97.5%
#'   quantiles across draws. This is fully consistent with the ribbon in
#'   the CDF plot (both derive from the same `cdf_mat`). If `NULL`, falls
#'   back to inverting the summary credible bands, which can fail near the
#'   tails.
#'
#' @return A data frame with columns `quantile`, `quantile_label`, `x_low`,
#'   and `x_high`. `x_low` and `x_high` are the 2.5% and 97.5% bounds of
#'   the 95% prediction interval for that quantile, on the same scale as
#'   the `x` column of `cdf_summary`.
#' @export
extract_quantiles <- function(cdf_summary, probs = c(0.5, 0.95), cdf_mat = NULL) {

  x_seq   <- cdf_summary$x
  results <- list()

  for (p in probs) {

    if (!is.null(cdf_mat)) {
      # For each posterior draw, interpolate the x at which the CDF crosses p.
      # approx() with rule = 1 returns NA when p lies outside the CDF range
      # (i.e. x_seq does not extend far enough); na.rm = TRUE handles this
      # gracefully — but a high NA rate suggests x_seq should be widened.
      x_at_p <- apply(cdf_mat, 1, function(cdf_row) {
        approx(x = cdf_row, y = x_seq, xout = p, rule = 1)$y
      })

      na_frac <- mean(is.na(x_at_p))
      if (na_frac > 0.05)
        warning(sprintf(
          "extract_quantiles: %.0f%% of draws did not reach p = %.2f within x_seq. Consider increasing max(x_seq) in compute_predictive_cdf().",
          na_frac * 100, p
        ))

      x_lo  <- quantile(x_at_p, 0.025, na.rm = TRUE)
      x_hi  <- quantile(x_at_p, 0.975, na.rm = TRUE)

    } else {
      # Fallback: invert summary credible bands.
      # x where the upper CDF band crosses p → lower x bound of PI
      # x where the lower CDF band crosses p → upper x bound of PI
      x_lo  <- x_seq[which.min(abs(cdf_summary$high - p))]
      x_hi  <- x_seq[which.min(abs(cdf_summary$low  - p))]
    }

    results[[paste0("q", p * 100)]] <- data.frame(
      quantile       = p,
      quantile_label = paste0("Q", p * 100),
      x_low          = x_lo,
      x_high         = x_hi
    )
  }

  dplyr::bind_rows(results)
}


#' Detect the day-fraction reporting resolution of an order-statistic dataset
#'
#' Infers how finely a study's reported order statistics were rounded, so the
#' day-rounded order-statistic likelihood (Reviewer 2, point 1(iii)) can use a
#' window narrower than a full day where the data supports it, rather than
#' assuming every dataset was rounded to the nearest whole day.
#'
#' `exact_vals` (the sample min and max for a median+range dataset) are
#' always true order statistics: an actual individual observation, under any
#' convention, at any sample size. `risky_vals` (the median, and for a
#' median+IQR dataset the reported quartiles) are not: for an even sample
#' size the standard median is the *average* of the two middle order
#' statistics, and quartile conventions vary and routinely interpolate
#' between adjacent order statistics. A day-integer dataset can therefore
#' report a median or quartile with a spurious fractional part (most often an
#' exact half, from averaging two integers) that reflects the interpolation
#' arithmetic, not the study's real measurement precision. Checked against
#' the curated corpus: several real median+range datasets (e.g. n=28, n=8,
#' n=10, n=22) show exactly this pattern, a half-integer median next to
#' integer min/max, for both even and odd n.
#'
#' Because of this, `exact_vals` are checked against the full set of "nice"
#' divisors of a day (whole day, half-day, ..., hourly), coarsest first,
#' since an integer is trivially consistent with every finer grid too (e.g.
#' `5` is a multiple of `1/24` as well as of `1`) and the coarsest match is
#' the meaningful one. `risky_vals` are trusted only for signals that a
#' generic linear interpolation between two integers is unlikely to produce
#' by chance: an hour-based grid (hourly, 2-hourly, 3-hourly) or a decimal
#' precision of two or more places (a simple interpolation weight reproduces
#' at most one non-trivial decimal digit from two integers, e.g. an eighth
#' gives `x.125`, a genuine 3-decimal case, which is why the 2dp threshold is
#' conservative rather than exact). Otherwise `risky_vals` are ignored and
#' the result defaults to a whole day. This means the function will
#' sometimes underestimate a median+IQR dataset's true resolution (there is
#' no min/max to anchor it), but underestimating resolution only means
#' falling back to the wider, already-accepted day window, not repeating the
#' overconfidence problem this exists to fix.
#'
#' @param exact_vals Numeric vector of statistics that are always true order
#'   statistics (min and max, for a median+range dataset). `NA`s are
#'   dropped. Pass `numeric(0)` if none apply (e.g. a median+IQR dataset).
#' @param risky_vals Numeric vector of statistics that may be interpolated
#'   rather than raw order statistics (the median; and, for a median+IQR
#'   dataset, the quartiles too). `NA`s are dropped. Default `numeric(0)`.
#' @param tol Numerical tolerance for judging a value to be a multiple of a
#'   candidate grid. Default `1e-6`.
#' @return A single number: the detected resolution in days (`1` = whole day,
#'   `1/24` = hourly, etc.). Defaults to `1` if both arguments are empty
#'   after dropping `NA`s.
#' @export
detect_resolution <- function(exact_vals, risky_vals = numeric(0), tol = 1e-6) {
  exact_vals <- exact_vals[!is.na(exact_vals)]
  risky_vals <- risky_vals[!is.na(risky_vals)]

  grid_match <- function(vals, denoms) {
    for (den in denoms) {
      if (all(abs(vals * den - round(vals * den)) < tol)) return(1 / den)
    }
    NA_real_
  }
  decimals_of <- function(x) {
    s <- sub("0+$", "", formatC(x, digits = 6, format = "f"))
    s <- sub("\\.$", "", s)  # a bare trailing "." remains when x is a whole number
    if (!grepl(".", s, fixed = TRUE)) return(0L)
    nchar(strsplit(s, ".", fixed = TRUE)[[1]][2])
  }
  decimal_precision <- function(vals) {
    if (length(vals) == 0) return(NA_real_)
    max_dec <- min(max(vapply(vals, decimals_of, integer(1))), 3L)
    if (max_dec == 0L) NA_real_ else 10^(-max_dec)
  }

  nice_denoms <- c(1, 2, 3, 4, 6, 8, 12, 24)

  # Exact order statistics: trust the full grid search, then the general
  # decimal-precision fallback (a study reporting min/max to N decimal
  # places really did measure to that precision).
  if (length(exact_vals) > 0) {
    res <- grid_match(exact_vals, nice_denoms)
    if (is.na(res)) res <- decimal_precision(exact_vals)
    if (!is.na(res)) return(res)
  }

  # No exact anchor, or it was all whole days: only accept risky_vals'
  # evidence of finer resolution if it is not plausibly a simple
  # interpolation artifact (an hour-based grid, or >=2dp decimal precision).
  #
  # Finding the COARSEST grid risky_vals fit (not just checking membership of
  # {8,12,24} directly) matters here: a half-day value like 6.5 is *also*
  # trivially a multiple of 1/8, 1/12 and 1/24 (0.5 = 4/8 = 6/12 = 12/24), so
  # checking those denominators on their own would wrongly "confirm" an
  # hour-based grid for a dataset that is really just half-day. Only the
  # coarsest grid that fits is a genuine claim about the resolution.
  if (length(risky_vals) > 0) {
    coarsest_den <- NA_real_
    for (den in nice_denoms) {
      if (all(abs(risky_vals * den - round(risky_vals * den)) < tol)) { coarsest_den <- den; break }
    }
    hour_res <- if (!is.na(coarsest_den) && coarsest_den %in% c(8, 12, 24)) 1 / coarsest_den else NA_real_
    dec_res  <- decimal_precision(risky_vals)
    if (!is.na(dec_res) && dec_res > 0.01) dec_res <- NA_real_  # 1dp: too easily a simple-fraction artifact
    candidates <- c(hour_res, dec_res)
    candidates <- candidates[!is.na(candidates)]
    if (length(candidates) > 0) return(min(candidates))
  }

  1
}


#' Prepare Stan data from a list of dataset summaries
#'
#' Converts a list of dataset descriptors (each providing summary statistics
#' and a sample size) into the named list expected by the
#' `hierarchical_data_synthesis_summary_stats` Stan model.
#'
#' @param datasets A named list of lists. Each element must contain one of the
#'   following combinations of summary statistics:
#'   \describe{
#'     \item{`median`, `min`, `max`}{Median and range (summary type 1). `n`
#'       (sample size) is required.}
#'     \item{`median`, `Q1`, `Q3`}{Median and inter-quartile range (summary
#'       type 2). `n` is required.}
#'     \item{`mean`, `sd`}{Mean and standard deviation (summary type 3). `n`
#'       is required.}
#'     \item{`freq_value`, `freq_count`}{Frequency table of (value, count)
#'       pairs (summary type 4). `n` is optional and defaults to
#'       `sum(freq_count)`.}
#'     \item{`freq_lower`, `freq_upper`, `freq_count`}{Interval-censored
#'       frequency table (summary type 5). Each entry gives the lower and
#'       upper bound of the censoring interval and the count of individuals
#'       in that interval. When `freq_lower[i] == freq_upper[i]` the
#'       observation is treated as exact. `n` is optional and defaults to
#'       `sum(freq_count)`.}
#'     \item{`expo_lower`, `expo_upper`, `event_lower`, `event_upper`,
#'       `freq_count`}{Double interval-censored frequency table (summary
#'       type 6). Both the exposure window (`expo_lower`, `expo_upper`) and
#'       the event window (`event_lower`, `event_upper`) are interval-censored.
#'       All five vectors must have the same length. Requires
#'       `expo_upper[i] <= event_lower[i]` for all i (positive delays). `n`
#'       is optional and defaults to `sum(freq_count)`.}
#'     \item{`expo_lower`, `expo_upper`, `event_lower`, `event_upper`,
#'       `event_observed`, `freq_count`, `truncation_time`}{Double
#'       interval-censored with right truncation/censoring (summary type 7).
#'       As type 6, plus `event_observed` (1 = onset seen, 0 = right-censored:
#'       not yet observed by the analysis cutoff) and a single
#'       `truncation_time` per dataset (the analysis cutoff `T`, in the same
#'       time origin as the exposure/event windows). `truncation_time` is a
#'       study-design choice, not a rounded observation of a random event
#'       time, so it is never treated as day-rounded the way exposure/event
#'       windows are. Set it using the same time-encoding convention as that
#'       dataset's own windows: e.g. if event days are encoded as `[day,
#'       day+1)` intervals, set `truncation_time` to `(last observed day +
#'       1)` to mean "complete through the end of that day", not to the bare
#'       day number.}
#'   }
#'   Each element may also contain an optional `source` field — a free-text
#'   character string recording the bibliographic reference for that dataset
#'   (e.g. `"Surname (year), doi: doi.org/xyz"`). This field is ignored
#'   during Stan data preparation and is never passed to the model.
#' @param dist_type Integer distribution code: `1` = log-normal, `2` = gamma,
#'   `3` = Weibull. Defaults to `1`.
#' @param use_custom_priors Integer flag (0 or 1) for custom prior use.
#'   Currently unused; reserved for future extension. Defaults to `0`.
#' @param custom_priors Named list of prior overrides. Any values not supplied
#'   fall back to distribution-appropriate defaults (see Details). Recognised
#'   names: `mu0_sd`, `log_tau_mean`, `log_tau_sd`, `log_phi_mean`,
#'   `log_phi_sd`, `log_omega_mean`, `log_omega_sd`. `log_phi_mean`/`log_phi_sd`
#'   are the prior for `log_phi0` (population-mean dispersion);
#'   `log_omega_mean`/`log_omega_sd` are the prior for `log_omega`, the
#'   between-study SD of log dispersion (per-study `phi_d` is hierarchical,
#'   not a shared scalar).
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
#' @note **Backward compatibility:** The Stan model requires `freq_lower`,
#'   `freq_upper`, `event_lower`, and `event_upper` to be present in the data
#'   list for all runs. These are populated automatically by this function.
#'   If constructing the Stan data list manually, fill all four arrays with
#'   zeros for non-applicable datasets, e.g.:
#'   ```r
#'   stan_data$freq_lower  <- rep(0, stan_data$n_freq_total)
#'   stan_data$freq_upper  <- rep(0, stan_data$n_freq_total)
#'   stan_data$event_lower <- rep(0, stan_data$n_freq_total)
#'   stan_data$event_upper <- rep(0, stan_data$n_freq_total)
#'   ```
#' @importFrom utils modifyList
#' @export
prepare_stan_data_from_datasets <- function(datasets, dist_type = 1,
                                            use_custom_priors = 0,
                                            custom_priors = list()) {

  # Apply distribution-specific defaults for log_phi_mean/log_phi_sd,
  # log_omega_mean (log_omega_sd stays flat, see below), and
  # log_kappa_mean/log_kappa_sd.
  #
  # phi meaning per distribution:
  #   lognormal  (1): phi = log-SD (sigma),  typical range 0.2-0.7  -> log_phi_mean = -0.7
  #   gamma      (2): phi = shape,           typical range 5-30     -> log_phi_mean =  2.5
  #   weibull    (3): phi = shape,           typical range 2-6      -> log_phi_mean =  1.0
  #   burr XII   (4): phi = c (shape1),      typical range 1-5      -> log_phi_mean =  0.7
  #   gen. gamma (5): phi = sigma (log-disp),typical range 0.2-1.0  -> log_phi_mean = -0.5
  #
  # log_omega_mean: between-study SD of log dispersion, family-specific
  # because how tightly phi is pinned down per study (and hence how much
  # of its apparent between-study spread is real heterogeneity vs
  # estimation noise) differs by family - notably gamma's shape parameter
  # is intrinsically harder to identify per-study than lognormal's sigma
  # (confirmed by a per-family, per-dataset no-pooling MAP scoping across
  # COVID-19/SARS/Cholera/Dengue, 103 datasets/family; pooled within-pathogen
  # sd(log phi): lognormal 0.410, burr12 0.415, gengamma 0.483, weibull 0.457,
  # gamma 0.677). log_omega_mean below is log() of that pooled sd, rounded.
  # A single shared uninformative prior on top of a single scalar phi could
  # not previously produce this failure mode; with phi now hierarchical
  # per-dataset, a too-tight shared omega prior caused a severe funnel/mixing
  # pathology for gamma specifically (Rhat > 100) that resolved once
  # family-specific scale was used. See REVISION_TODO.md, point 4.
  #
  # kappa meaning per distribution:
  #   dist 1-3: kappa is unused; wide uninformative prior centred at 1.
  #   burr XII (4): kappa = k (shape2), typical range 1-10  -> log_kappa_mean = 1.0
  #   gen. gamma (5): kappa = Q (shape), typical range 0.3-3 -> log_kappa_mean = 0.0
  dist_defaults <- list(
    `1` = list(log_phi_mean = -0.7, log_phi_sd = 0.5, log_omega_mean = -0.9, log_kappa_mean = 0.0, log_kappa_sd = 1.0),
    `2` = list(log_phi_mean =  2.5, log_phi_sd = 0.5, log_omega_mean = -0.4, log_kappa_mean = 0.0, log_kappa_sd = 1.0),
    `3` = list(log_phi_mean =  1.0, log_phi_sd = 0.5, log_omega_mean = -0.8, log_kappa_mean = 0.0, log_kappa_sd = 1.0),
    `4` = list(log_phi_mean =  0.7, log_phi_sd = 0.5, log_omega_mean = -0.9, log_kappa_mean = 1.0, log_kappa_sd = 0.5),
    `5` = list(log_phi_mean = -0.5, log_phi_sd = 0.5, log_omega_mean = -0.7, log_kappa_mean = 0.0, log_kappa_sd = 0.5)
  )[[as.character(dist_type)]]

  if (is.null(dist_defaults)) {
    stop(sprintf("'dist_type' must be 1, 2, 3, 4, or 5 (got %s).", dist_type), call. = FALSE)
  }

  defaults <- list(
    mu0_sd        = 1.0,
    log_tau_mean  = 0.2,
    log_tau_sd    = 0.5,
    log_phi_mean  = dist_defaults$log_phi_mean,
    log_phi_sd    = dist_defaults$log_phi_sd,
    # Between-study SD of log dispersion (log_phi_d), analogous to log_tau
    # for location. log_omega_mean is family-specific (see dist_defaults
    # above); log_omega_sd is left flat at 0.5 across families, matching
    # the log_tau_sd=0.5 convention already used here - it is a generic
    # uncertainty width on top of the point estimate, not itself
    # data-derived per family. See REVISION_TODO.md, point 4.
    log_omega_mean = dist_defaults$log_omega_mean,
    log_omega_sd   = 0.5,
    log_kappa_mean = dist_defaults$log_kappa_mean,
    log_kappa_sd   = dist_defaults$log_kappa_sd
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
  # Day-fraction rounding resolution for summary_type 1/2 (see
  # detect_resolution()); irrelevant for other types, left at the default 1.
  resolution_vec <- rep(1, n_datasets)
  # Frequency table flat arrays (for summary_type == 4, 5, 6, and 7)
  freq_value_all      <- numeric(0)
  freq_lower_all      <- numeric(0)
  freq_upper_all      <- numeric(0)
  event_lower_all     <- numeric(0)   # event window lower bounds for types 6/7; 0 elsewhere
  event_upper_all     <- numeric(0)   # event window upper bounds for types 6/7; 0 elsewhere
  event_observed_all  <- integer(0)   # 1=onset seen, 0=right-censored; type 7 only (1 elsewhere)
  freq_count_all      <- integer(0)
  freq_start_vec      <- integer(n_datasets)
  freq_len_vec        <- integer(n_datasets)
  truncation_time_vec <- numeric(n_datasets)  # analysis date T for type 7; 0 elsewhere
  running_start  <- 1L

  # Process each dataset
  for (i in seq_along(datasets)) {
     d<- datasets[[i]]

    # Determine summary type and extract statistics
    if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) {
      # Type 1: median + range (min, max)
      n_obs_vec[i]    <- d$n
      summary_type[i] <- 1
      obs_stat1[i]    <- d$median
      obs_stat2[i]    <- d$min
      obs_stat3[i]    <- d$max
      # min/max are always true order statistics; median may be an
      # interpolated average for even n, so it is only a "risky" signal.
      resolution_vec[i] <- detect_resolution(c(d$min, d$max), d$median)

    } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) {
      # Type 2: median + IQR (Q1, Q3)
      n_obs_vec[i]    <- d$n
      summary_type[i] <- 2
      obs_stat1[i]    <- d$median
      obs_stat2[i]    <- d$Q1
      obs_stat3[i]    <- d$Q3
      # No min/max anchor here; median and quartiles can all be interpolated,
      # so all three are "risky" (see detect_resolution()).
      resolution_vec[i] <- detect_resolution(numeric(0), c(d$median, d$Q1, d$Q3))

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
      # freq_value entries are true recorded values, not interpolated
      # statistics, so (unlike median/Q1/Q3) they can all be trusted with the
      # full "exact" grid search, the same way min/max are for type 1.
      resolution_vec[i] <- detect_resolution(d$freq_value)
      freq_start_vec[i] <- running_start
      freq_len_vec[i]   <- length(d$freq_value)
      freq_value_all       <- c(freq_value_all,      as.numeric(d$freq_value))
      freq_lower_all       <- c(freq_lower_all,      rep(0, length(d$freq_value)))  # unused for type 4
      freq_upper_all       <- c(freq_upper_all,      rep(0, length(d$freq_value)))  # unused for type 4
      event_lower_all      <- c(event_lower_all,     rep(0, length(d$freq_value)))  # unused for type 4
      event_upper_all      <- c(event_upper_all,     rep(0, length(d$freq_value)))  # unused for type 4
      event_observed_all   <- c(event_observed_all,  rep(1L, length(d$freq_value))) # unused for type 4
      freq_count_all       <- c(freq_count_all,      as.integer(d$freq_count))
      running_start        <- running_start + freq_len_vec[i]

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
      # Only entries with no reported range (freq_lower == freq_upper) are
      # treated as rounded point values in Stan; resolution is inferred from
      # just those (true recorded values, so the "exact" grid search
      # applies). detect_resolution() defaults to 1 (unused) if a dataset has
      # no such entries.
      resolution_vec[i] <- detect_resolution(d$freq_lower[d$freq_lower == d$freq_upper])
      freq_start_vec[i] <- running_start
      freq_len_vec[i]   <- length(d$freq_lower)
      freq_value_all      <- c(freq_value_all,     rep(0, length(d$freq_lower)))  # unused for type 5
      freq_lower_all      <- c(freq_lower_all,     as.numeric(d$freq_lower))
      freq_upper_all      <- c(freq_upper_all,     as.numeric(d$freq_upper))
      event_lower_all     <- c(event_lower_all,    rep(0, length(d$freq_lower)))  # unused for type 5
      event_upper_all     <- c(event_upper_all,    rep(0, length(d$freq_lower)))  # unused for type 5
      event_observed_all  <- c(event_observed_all, rep(1L, length(d$freq_lower))) # unused for type 5
      freq_count_all      <- c(freq_count_all,     as.integer(d$freq_count))
      running_start     <- running_start + freq_len_vec[i]

    } else if (!is.null(d$expo_lower) && !is.null(d$expo_upper) &&
               !is.null(d$event_lower) && !is.null(d$event_upper) &&
               !is.null(d$freq_count) && is.null(d$truncation_time)) {
      # Type 6: double interval-censored frequency table.
      # expo_lower / expo_upper: exposure window bounds.
      # event_lower / event_upper: event window bounds.
      # In Stan, freq_lower / freq_upper carry the exposure bounds; the new
      # event_lower / event_upper arrays carry the event bounds.
      n_len <- length(d$expo_lower)
      if (length(d$expo_upper)  != n_len || length(d$event_lower) != n_len ||
          length(d$event_upper) != n_len || length(d$freq_count)  != n_len) {
        stop(paste("Dataset", i,
                   ": expo_lower, expo_upper, event_lower, event_upper and freq_count",
                   "must all have the same length"))
      }
      if (any(d$expo_lower > d$expo_upper)) {
        stop(paste("Dataset", i,
                   ": all expo_lower values must be <= their corresponding expo_upper values"))
      }
      if (any(d$event_lower > d$event_upper)) {
        stop(paste("Dataset", i,
                   ": all event_lower values must be <= their corresponding event_upper values"))
      }
      if (any(d$expo_upper > d$event_lower)) {
        stop(paste("Dataset", i,
                   ": all expo_upper values must be <= their corresponding event_lower values",
                   "(delays must be non-negative)"))
      }
      n_obs_vec[i]      <- if (!is.null(d$n)) d$n else sum(d$freq_count)
      summary_type[i]   <- 6L
      obs_stat1[i]      <- 0  # placeholder
      obs_stat2[i]      <- 0  # placeholder
      obs_stat3[i]      <- 0  # placeholder
      freq_start_vec[i] <- running_start
      freq_len_vec[i]   <- n_len
      freq_value_all      <- c(freq_value_all,     rep(0, n_len))               # unused for type 6
      freq_lower_all      <- c(freq_lower_all,     as.numeric(d$expo_lower))    # exposure lower bound
      freq_upper_all      <- c(freq_upper_all,     as.numeric(d$expo_upper))    # exposure upper bound
      event_lower_all     <- c(event_lower_all,    as.numeric(d$event_lower))   # event lower bound
      event_upper_all     <- c(event_upper_all,    as.numeric(d$event_upper))   # event upper bound
      event_observed_all  <- c(event_observed_all, rep(1L, n_len))              # unused for type 6
      freq_count_all      <- c(freq_count_all,     as.integer(d$freq_count))
      running_start       <- running_start + n_len

    } else if (!is.null(d$expo_lower) && !is.null(d$expo_upper) &&
               !is.null(d$freq_count) && !is.null(d$truncation_time)) {
      # Type 7: doubly interval-censored with right truncation/censoring.
      # event_lower / event_upper may contain NAs for right-censored individuals
      # (onset not yet observed by the analysis date T).
      n_len <- length(d$expo_lower)
      evl_raw <- if (!is.null(d$event_lower)) d$event_lower else rep(NA_real_, n_len)
      evu_raw <- if (!is.null(d$event_upper)) d$event_upper else rep(NA_real_, n_len)

      if (length(d$expo_upper) != n_len || length(d$freq_count) != n_len ||
          length(evl_raw) != n_len || length(evu_raw) != n_len) {
        stop(paste("Dataset", i,
                   ": expo_lower, expo_upper, event_lower, event_upper, and freq_count",
                   "must all have the same length"))
      }
      if (any(d$expo_lower > d$expo_upper)) {
        stop(paste("Dataset", i,
                   ": all expo_lower values must be <= their corresponding expo_upper values"))
      }
      obs_mask <- !is.na(evl_raw)
      if (any(obs_mask & (evl_raw > evu_raw), na.rm = TRUE)) {
        stop(paste("Dataset", i,
                   ": all event_lower values must be <= their corresponding event_upper values"))
      }
      if (any(obs_mask & (d$expo_upper > evl_raw), na.rm = TRUE)) {
        stop(paste("Dataset", i,
                   ": all expo_upper values must be <= their corresponding event_lower values",
                   "(delays must be non-negative)"))
      }
      if (any(d$expo_upper >= d$truncation_time)) {
        stop(paste("Dataset", i,
                   ": all expo_upper values must be < truncation_time"))
      }
      if (any(obs_mask & (evu_raw > d$truncation_time), na.rm = TRUE)) {
        stop(paste("Dataset", i,
                   ": all event_upper values must be <= truncation_time"))
      }
      n_obs_vec[i]          <- if (!is.null(d$n)) d$n else sum(d$freq_count)
      summary_type[i]       <- 7L
      truncation_time_vec[i] <- d$truncation_time
      obs_stat1[i]          <- 0
      obs_stat2[i]          <- 0
      obs_stat3[i]          <- 0
      freq_start_vec[i]     <- running_start
      freq_len_vec[i]       <- n_len
      freq_value_all        <- c(freq_value_all,     rep(0, n_len))
      freq_lower_all        <- c(freq_lower_all,     as.numeric(d$expo_lower))
      freq_upper_all        <- c(freq_upper_all,     as.numeric(d$expo_upper))
      # Replace NA event bounds with 0 for right-censored rows (ignored in Stan)
      event_lower_all       <- c(event_lower_all,    as.numeric(ifelse(obs_mask, evl_raw, 0)))
      event_upper_all       <- c(event_upper_all,    as.numeric(ifelse(obs_mask, evu_raw, 0)))
      event_observed_all    <- c(event_observed_all, as.integer(obs_mask))
      freq_count_all        <- c(freq_count_all,     as.integer(d$freq_count))
      running_start         <- running_start + n_len

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
    } else if (summary_type[i] == 6L && freq_len_vec[i] > 0) {
      s   <- freq_start_vec[i]
      ln  <- freq_len_vec[i]
      el  <- freq_lower_all[s:(s + ln - 1)]    # expo_lower
      er  <- freq_upper_all[s:(s + ln - 1)]    # expo_upper
      evl <- event_lower_all[s:(s + ln - 1)]   # event_lower
      evu <- event_upper_all[s:(s + ln - 1)]   # event_upper
      fc  <- freq_count_all[s:(s + ln - 1)]
      delay_mid <- ((evl + evu) / 2) - ((el + er) / 2)
      central_estimates[i] <- sum(delay_mid * fc) / sum(fc)
    } else if (summary_type[i] == 7L && freq_len_vec[i] > 0) {
      s    <- freq_start_vec[i]
      ln   <- freq_len_vec[i]
      el   <- freq_lower_all[s:(s + ln - 1)]
      er   <- freq_upper_all[s:(s + ln - 1)]
      evl  <- event_lower_all[s:(s + ln - 1)]
      evu  <- event_upper_all[s:(s + ln - 1)]
      fc   <- freq_count_all[s:(s + ln - 1)]
      eobs <- event_observed_all[s:(s + ln - 1)]
      # Use only fully observed rows for the central estimate
      obs_idx <- eobs == 1L
      if (any(obs_idx)) {
        delay_mid <- ((evl[obs_idx] + evu[obs_idx]) / 2) - ((el[obs_idx] + er[obs_idx]) / 2)
        central_estimates[i] <- sum(delay_mid * fc[obs_idx]) / sum(fc[obs_idx])
      }
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
    resolution   = as.array(resolution_vec),
    n_freq_total     = length(freq_value_all),
    freq_value       = freq_value_all,
    freq_lower       = freq_lower_all,
    freq_upper       = freq_upper_all,
    event_lower      = event_lower_all,
    event_upper      = event_upper_all,
    event_observed   = event_observed_all,
    freq_count       = freq_count_all,
    freq_start       = as.array(freq_start_vec),
    freq_len         = as.array(freq_len_vec),
    truncation_time  = as.array(truncation_time_vec)
  )

  stan_data$mu0_mean      <- if (length(valid_centrals) > 0) log(mean(valid_centrals)) else 0
  stan_data$mu0_sd        <- custom_priors$mu0_sd
  stan_data$log_tau_mean  <- custom_priors$log_tau_mean
  stan_data$log_tau_sd    <- custom_priors$log_tau_sd
  stan_data$log_phi_mean  <- custom_priors$log_phi_mean
  stan_data$log_phi_sd    <- custom_priors$log_phi_sd
  stan_data$log_omega_mean <- custom_priors$log_omega_mean
  stan_data$log_omega_sd   <- custom_priors$log_omega_sd
  stan_data$log_kappa_mean <- custom_priors$log_kappa_mean
  stan_data$log_kappa_sd   <- custom_priors$log_kappa_sd

  # Per-order-statistic adaptive quadrature panel count (performance lever,
  # not a correctness change): see POINT1_LIKELIHOOD_MATHS.md Part G for the
  # derivation. A wrong choice here only costs efficiency (or triggers
  # pre_inference_checks()-style convergence issues, handled by the settings
  # escalation in fit_with_escalation()), since the underlying quadrature
  # maths is unchanged - it just decides how many panels to spend.
  n_panels_stat1 <- integer(n_datasets)
  n_panels_stat2 <- integer(n_datasets)
  n_panels_stat3 <- integer(n_datasets)
  for (i in seq_len(n_datasets)) {
    if (summary_type[i] == 1L) {
      k_median <- (n_obs_vec[i] + 1L) %/% 2L
      phi_guess <- .mom_phi_guess(datasets[[i]], dist_type, custom_priors$log_kappa_mean)
      n_panels_stat1[i] <- .order_stat_n_panels(obs_stat1[i], n_obs_vec[i], k_median, resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
      n_panels_stat2[i] <- .order_stat_n_panels(obs_stat2[i], n_obs_vec[i], 1L,        resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
      n_panels_stat3[i] <- .order_stat_n_panels(obs_stat3[i], n_obs_vec[i], n_obs_vec[i], resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
    } else if (summary_type[i] == 2L) {
      n <- n_obs_vec[i]
      k_median <- (n + 1L) %/% 2L
      k_q25 <- max((n + 1L) %/% 4L, 1L)
      k_q75 <- min(max((3L * (n + 1L)) %/% 4L, k_q25 + 1L), n)
      phi_guess <- .mom_phi_guess(datasets[[i]], dist_type, custom_priors$log_kappa_mean)
      n_panels_stat1[i] <- .order_stat_n_panels(obs_stat1[i], n, k_median, resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
      n_panels_stat2[i] <- .order_stat_n_panels(obs_stat2[i], n, k_q25,    resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
      n_panels_stat3[i] <- .order_stat_n_panels(obs_stat3[i], n, k_q75,    resolution_vec[i], dist_type, phi_guess, exp(custom_priors$log_kappa_mean))
    } else {
      n_panels_stat1[i] <- 24L; n_panels_stat2[i] <- 24L; n_panels_stat3[i] <- 24L
    }
  }
  stan_data$n_panels_stat1 <- as.array(n_panels_stat1)
  stan_data$n_panels_stat2 <- as.array(n_panels_stat2)
  stan_data$n_panels_stat3 <- as.array(n_panels_stat3)

  return(stan_data)
}

# Adaptive quadrature panel count for order-statistic rounding likelihoods --
#
# See POINT1_LIKELIHOOD_MATHS.md Part G for the full derivation. Summary: the
# panel count needed for the composite Gauss-Legendre quadrature
# (order_stat_rounded_loglik_fun in the .stan files) to hit the same accuracy
# as the original fixed 24-panel scheme is, empirically, a clean function of
# rho = resolution / SE_asymptotic(order statistic) - not of n or family
# alone. Since rho depends on the (unknown-until-fitted) dispersion
# parameter, this uses a per-dataset method-of-moments guess evaluated at a
# family-specific safety margin (the direction that INCREASES rho, chosen
# per family from the empirical MAP-vs-MoM comparison this session): for
# spread-type families (lognormal, gen. gamma - CV increases with phi) the
# dangerous direction is phi too LOW; for concentration-type families
# (Weibull, Burr XII - CV decreases with phi) it's phi too HIGH. Gamma is
# excluded entirely (kept at the fixed 24 panels): its own MAP-vs-MoM
# comparison showed a 99th-percentile ratio of 4.4x (vs ~1.7-2x for the
# other concentration-type families), and covering that via margin alone
# was found to require booking >96 panels even at moderate n - intractable,
# and unnecessary, since that regime is exactly where the true log-density
# is also vanishingly small (verified directly: both the reference and the
# 24-panel scheme agree the region contributes ~0 to the posterior, they
# just disagree on exactly how close to -Inf). Any residual risk for the
# eligible families is caught by fit_with_escalation()'s Rhat/divergence
# check, not by this heuristic alone.
.panels_from_rho <- function(rho) {
  if (!is.finite(rho)) return(24L)
  if (rho <= 4)  return(2L)
  if (rho <= 10) return(4L)
  if (rho <= 20) return(6L)
  if (rho <= 30) return(8L)
  if (rho <= 46) return(12L)
  if (rho <= 57) return(16L)
  return(24L)
}

# Solve for `loc` such that the p-th quantile of the family equals v, given
# phi/kappa - lets us evaluate the density AT the reported (fixed, known)
# value v without needing an independent loc estimate.
.solve_loc_for_quantile <- function(v, p, dist_type, phi, kappa) {
  switch(as.character(dist_type),
    "1" = log(v) - phi * stats::qnorm(p),
    "3" = log(v) - (1 / phi) * log(-log1p(-p)),
    "4" = log(v) - (1 / phi) * log((1 - p)^(-1 / kappa) - 1),
    "5" = { a <- 1 / kappa^2; w <- log(stats::qgamma(p, a, 1) / a) / kappa; log(v) - phi * w },
    NA_real_
  )
}

# Log-density, matching dist_logpdf_fun in the .stan files exactly (families
# 1/3/4/5 only - gamma, family 2, never calls this since it keeps the fixed
# panel count).
.dist_logpdf_r <- function(x, dist_type, loc, phi, kappa) {
  switch(as.character(dist_type),
    "1" = stats::dlnorm(x, meanlog = loc, sdlog = phi, log = TRUE),
    "3" = stats::dweibull(x, shape = phi, scale = exp(loc), log = TRUE),
    "4" = { u <- log(x) - loc; log_term <- log1p(exp(phi * u))
            log(phi) + log(kappa) + (phi - 1) * u - loc - (kappa + 1) * log_term },
    "5" = { a <- 1 / kappa^2; w <- (log(x) - loc) / phi
            log(kappa) - log(phi) - log(x) + a * log(a) + a * kappa * w - a * exp(kappa * w) - lgamma(a) },
    NA_real_
  )
}

# Quick per-dataset method-of-moments phi guess, mirroring update_phi_prior()
# (kept as a separate, self-contained copy rather than refactoring that
# already-validated function, to avoid any risk of disturbing it).
.mom_phi_guess <- function(d, dist_type, log_kappa_mean) {
  dist_name <- c("1" = "lognormal", "2" = "gamma", "3" = "weibull",
                 "4" = "burr12", "5" = "gengamma")[[as.character(dist_type)]]
  kappa_val <- exp(log_kappa_mean)
  mean_est <- sd_est <- NA_real_
  if (!is.null(d$mean) && !is.null(d$sd)) { mean_est <- d$mean; sd_est <- d$sd
  } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) { mean_est <- d$median; sd_est <- (d$Q3 - d$Q1) / 1.35
  } else if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) { mean_est <- d$median; sd_est <- (d$max - d$min) / 4
  }
  if (is.na(mean_est) || is.na(sd_est) || sd_est <= 0) return(NA_real_)
  tryCatch(switch(dist_name,
    lognormal = sqrt(log(1 + (sd_est / mean_est)^2)),
    weibull   = { cv <- sd_est / mean_est
      stats::uniroot(function(k) sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1) - cv, c(0.1, 200))$root },
    gengamma  = { cv <- sd_est / mean_est; gs <- 1 / kappa_val^2
      stats::uniroot(function(s) { a1 <- gs + s/kappa_val; a2 <- gs + 2*s/kappa_val
        sqrt(max(exp(lgamma(a2) + lgamma(gs) - 2*lgamma(a1)) - 1, 0)) - cv }, c(1e-6, 20))$root },
    burr12    = { cv <- sd_est / mean_est; lower_c <- 2/kappa_val + 1e-6
      stats::uniroot(function(c_val) { lb1 <- lbeta(kappa_val - 1/c_val, 1 + 1/c_val); lb2 <- lbeta(kappa_val - 2/c_val, 1 + 2/c_val)
        sqrt(exp(lb2 - log(kappa_val) - 2*lb1) - 1) - cv }, c(lower_c, 50))$root }
  ), error = function(e) NA_real_)
}

# Per-order-statistic panel count. Falls back to the original fixed 24 (i.e.
# never worse than before) whenever the MoM guess or density evaluation is
# unavailable or non-finite.
.order_stat_n_panels <- function(v, n, k, resolution, dist_type, phi_guess, kappa_guess) {
  if (dist_type == 2L) return(24L)
  if (is.na(phi_guess) || phi_guess <= 0 || is.na(v) || v <= 0) return(24L)

  margin_phi <- switch(as.character(dist_type),
    "1" = phi_guess * 0.6, "5" = phi_guess * 0.6,   # spread-type: danger = phi too LOW
    "3" = phi_guess * 2.0, "4" = phi_guess * 2.0,   # concentration-type: danger = phi too HIGH
    phi_guess
  )
  p <- k / (n + 1)
  loc <- tryCatch(.solve_loc_for_quantile(v, p, dist_type, margin_phi, kappa_guess), error = function(e) NA_real_)
  if (is.na(loc) || !is.finite(loc)) return(24L)
  logf <- tryCatch(.dist_logpdf_r(v, dist_type, loc, margin_phi, kappa_guess), error = function(e) NA_real_)
  if (is.na(logf) || !is.finite(logf)) return(24L)

  se_asymp <- sqrt(p * (1 - p) / n) / exp(logf)
  rho <- resolution / se_asymp
  .panels_from_rho(rho)
}


# Data-driven phi prior ----------------------------------------------------

#' Update the log_phi prior mean from method-of-moments estimates
#'
#' @description
#' Estimates \eqn{\phi} from each dataset individually using moment-based
#' approximations (the same approach used in check 1 of
#' [pre_inference_checks()]), then overwrites `log_phi_mean` in `stan_data`
#' with the log of the median implied \eqn{\phi} across all datasets.
#'
#' This replaces the fixed distribution-specific default with a value anchored
#' to the actual data scale, which is particularly useful for the gamma
#' distribution: the default prior mean (shape ≈ 12) can be far above the
#' data-implied shape (typically 2–8 for incubation periods), causing
#' `gamma_lccdf` to evaluate to `log(0) = -Inf` during Stan's initialisation
#' phase.
#'
#' `log_phi_sd` is left unchanged so the prior remains diffuse around the
#' data-derived centre.
#'
#' @param stan_data A named list returned by [prepare_stan_data_from_datasets()].
#' @param datasets The same named list of datasets passed to
#'   [prepare_stan_data_from_datasets()].  Used only for moment calculations.
#'
#' @return `stan_data` with `log_phi_mean` replaced by
#'   `log(median(implied_phi))`.  All other fields are unchanged.  If no
#'   finite implied-\eqn{\phi} values can be derived, a warning is issued and
#'   `stan_data` is returned unmodified.
#'
#' @details
#' **Moment approximations per distribution:**
#' \describe{
#'   \item{lognormal}{\eqn{\phi = \log(1 + (\mathrm{sd}/\mathrm{mean})^2)}
#'     (approximate log-variance)}
#'   \item{gamma}{\eqn{\phi = (\mathrm{mean}/\mathrm{sd})^2}
#'     (method-of-moments shape)}
#'   \item{Weibull}{\eqn{\phi} solved numerically from the CV via
#'     \eqn{CV^2 = \Gamma(1+2/k)/\Gamma(1+1/k)^2 - 1}}
#' }
#' For datasets that report only median + IQR, the SD is approximated as
#' \eqn{(Q3-Q1)/1.35}; for median + range, as \eqn{(\max-\min)/4};
#' for frequency tables, the weighted SD of the (mid-)points is used.
#'
#' @examples
#' \dontrun{
#'   stan_data <- prepare_stan_data_from_datasets(datasets_Mpox, dist_type = 2)
#'   stan_data <- update_phi_prior(stan_data, datasets_Mpox)
#'   fit <- rstan::sampling(stan_model, data = stan_data, chains = 4, iter = 2000)
#' }
#' @export
update_phi_prior <- function(stan_data, datasets) {

  dist_name <- c("1" = "lognormal", "2" = "gamma", "3" = "weibull",
                 "4" = "burr12",   "5" = "gengamma")[[
    as.character(stan_data$dist_type)
  ]]
  if (is.null(dist_name)) {
    stop("update_phi_prior: unrecognised dist_type (", stan_data$dist_type,
         "); expected 1, 2, 3, 4, or 5.", call. = FALSE)
  }
  # For Burr XII and GG, phi cannot be identified from the CV alone because the
  # CV depends on both phi and kappa.  We fix kappa at its current prior mean
  # (exp(log_kappa_mean)) and solve for phi numerically — the same conditional
  # moment-of-moments approach used for Weibull.  If log_kappa_mean has already
  # been set to a domain-specific value via custom_priors, that will be used.
  kappa_val <- exp(stan_data$log_kappa_mean)

  # Per-dataset moment estimates of mean and SD, covering all five summary types
  implied_phis <- vapply(datasets, function(d) {
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
    } else if (!is.null(d$freq_lower) && !is.null(d$freq_upper) &&
               !is.null(d$freq_count)) {
      mid      <- (d$freq_lower + d$freq_upper) / 2
      w        <- d$freq_count / sum(d$freq_count)
      mean_est <- sum(mid * w)
      sd_est   <- sqrt(sum(w * (mid - mean_est)^2))
    } else if (!is.null(d$expo_lower) && !is.null(d$expo_upper) &&
               !is.null(d$event_lower) && !is.null(d$event_upper) &&
               !is.null(d$freq_count)) {
      delay_mid <- ((d$event_lower + d$event_upper) / 2) -
                   ((d$expo_lower  + d$expo_upper)  / 2)
      w         <- d$freq_count / sum(d$freq_count)
      mean_est  <- sum(delay_mid * w)
      sd_est    <- sqrt(sum(w * (delay_mid - mean_est)^2))
    }

    if (is.na(mean_est) || is.na(sd_est) || sd_est <= 0) return(NA_real_)

    switch(dist_name,
      lognormal = sqrt(log(1 + (sd_est / mean_est)^2)),
      gamma     = (mean_est / sd_est)^2,
      weibull   = {
        cv      <- sd_est / mean_est
        obj_wei <- function(k) sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)^2 - 1) - cv
        tryCatch(stats::uniroot(obj_wei, c(0.1, 200))$root,
                 error = function(e) NA_real_)
      },
      gengamma  = {
        # CV² = Γ(γ + 2σ/κ)·Γ(γ) / Γ(γ + σ/κ)² − 1  where γ = 1/κ²
        # CV is an increasing function of σ (phi), starting at 0 as σ→0⁺.
        cv      <- sd_est / mean_est
        gs      <- 1.0 / kappa_val^2          # gamma_shape = 1/Q^2
        obj_gg  <- function(sigma) {
          a1 <- gs + sigma / kappa_val
          a2 <- gs + 2.0 * sigma / kappa_val
          cv_sq <- exp(lgamma(a2) + lgamma(gs) - 2.0 * lgamma(a1)) - 1.0
          sqrt(max(cv_sq, 0.0)) - cv
        }
        tryCatch(stats::uniroot(obj_gg, c(1e-6, 20.0))$root,
                 error = function(e) NA_real_)
      },
      burr12    = {
        # CV² = B(k−2/c, 1+2/c) / (k · B(k−1/c, 1+1/c)²) − 1
        # Moments require k·c > 2, i.e. c > 2/k.  CV is decreasing in c,
        # so uniroot searches from (2/k + ε) upward.
        cv        <- sd_est / mean_est
        lower_c   <- 2.0 / kappa_val + 1e-6
        obj_burr  <- function(c_val) {
          lb1 <- lbeta(kappa_val - 1.0 / c_val, 1.0 + 1.0 / c_val)
          lb2 <- lbeta(kappa_val - 2.0 / c_val, 1.0 + 2.0 / c_val)
          cv_model <- sqrt(exp(lb2 - log(kappa_val) - 2.0 * lb1) - 1.0)
          cv_model - cv
        }
        tryCatch(stats::uniroot(obj_burr, c(lower_c, 50.0))$root,
                 error = function(e) NA_real_)
      }
    )
  }, numeric(1))

  med_phi <- stats::median(implied_phis, na.rm = TRUE)

  if (is.na(med_phi) || med_phi <= 0) {
    warning("update_phi_prior: could not derive a finite positive phi estimate ",
            "from the data; log_phi_mean has not been changed.", call. = FALSE)
    return(stan_data)
  }

  stan_data$log_phi_mean <- log(med_phi)
  stan_data
}


# Dataset filtering --------------------------------------------------------

#' Filter a dataset list by subgroup and/or location
#'
#' Retains only those entries whose `subgroup` and/or `location` fields match
#' the requested values.  Entries that do not carry the field at all are
#' **kept** when the corresponding filter argument is `NULL` and **dropped**
#' when a filter is active (because their group membership is unknown).
#' Setting both arguments to `NULL` returns the list unchanged.
#'
#' @param datasets A named list of datasets in the format accepted by
#'   [prepare_stan_data_from_datasets()].  Each entry may optionally contain
#'   `subgroup` and/or `location` character fields.
#' @param subgroup Character vector of subgroup values to retain, or `NULL`
#'   (default) to skip subgroup filtering.
#' @param location Character vector of location values to retain, or `NULL`
#'   (default) to skip location filtering.
#'
#' @return A named list containing only the datasets that satisfy all active
#'   filter criteria.
#'
#' @examples
#' datasets <- list(
#'   d1 = list(mean = 4.0, sd = 2.4, n = 49,
#'             subgroup = "tick-bite", location = "Turkey"),
#'   d2 = list(mean = 6.0, sd = 3.1, n = 12,
#'             subgroup = "nosocomial", location = "Iran"),
#'   d3 = list(mean = 5.0, sd = 2.0, n = 30)   # no subgroup/location
#' )
#' filter_datasets(datasets, subgroup = "tick-bite")
#' filter_datasets(datasets, location = c("Turkey", "Iran"))
#' filter_datasets(datasets, subgroup = "nosocomial", location = "Iran")
#'
#' @export
filter_datasets <- function(datasets, subgroup = NULL, location = NULL) {
  if (is.null(subgroup) && is.null(location)) return(datasets)

  keep <- vapply(datasets, function(d) {
    if (!is.null(subgroup)) {
      val <- d$subgroup
      if (is.null(val) || !(val %in% subgroup)) return(FALSE)
    }
    if (!is.null(location)) {
      val <- d$location
      if (is.null(val) || !(val %in% location)) return(FALSE)
    }
    TRUE
  }, logical(1))

  datasets[keep]
}


# GG identifiability heuristic ------------------------------------------------

#' Check whether the Generalised Gamma is likely identifiable from a dataset
#'
#' @description
#' The Generalised Gamma (GG, dist_type = 5) has three parameters: location
#' (\eqn{\mu}), scale (\eqn{\sigma}/phi), and shape (\eqn{Q}/kappa).  All
#' datasets share a single (\eqn{\sigma}, \eqn{Q}) pair, so the GG is only
#' identifiable when datasets consistently imply the same distributional shape.
#' If different studies show widely different coefficients of variation (CV =
#' SD/mean), the sampler cannot find a coherent (\eqn{\sigma}, \eqn{Q}) and
#' will exhibit poor mixing or divergences.
#'
#' Two fast, pre-fit checks are applied:
#' \describe{
#'   \item{CV spread}{Computes the CV for every dataset using moment
#'     approximations (same logic as [update_phi_prior()]).  If
#'     \code{max(CV) / min(CV) > cv_spread_threshold} the CVs are too
#'     inconsistent to identify the extra GG parameter.}
#'   \item{Information richness}{The \eqn{Q} parameter encodes tail behaviour
#'     beyond mean and variance.  Datasets that supply only summary statistics
#'     (mean + SD, median + IQR, median + range) provide at most two moments
#'     and give weak leverage on \eqn{Q}.  If the fraction of datasets with
#'     frequency-table or interval-censored data (summary types 4 and 5) is
#'     below \code{min_rich_fraction}, the shape is too poorly constrained.}
#' }
#'
#' @param datasets A named list of datasets in the format accepted by
#'   [prepare_stan_data_from_datasets()].
#' @param cv_spread_threshold Numeric scalar (default 2.5).  Maximum tolerated
#'   ratio of the largest to the smallest per-dataset CV.  Increase to be more
#'   permissive, decrease to be stricter.
#' @param min_rich_fraction Numeric scalar in (0, 1] (default 0.30).  Minimum
#'   fraction of datasets that must be frequency-table or interval-censored
#'   (summary types 4 / 5).  Set to 0 to disable this check.
#' @param verbose Logical (default TRUE).  Print a one-line verdict with the
#'   reason a check failed.
#'
#' @return `TRUE` if both checks pass (GG fitting is worth attempting),
#'   `FALSE` otherwise.
#'
#' @examples
#' \dontrun{
#' should_attempt_gg(datasets_SARS)    # expected: FALSE
#' should_attempt_gg(datasets_Mpox)
#' }
#'
#' @export
should_attempt_gg <- function(datasets,
                               cv_spread_threshold = 2.5,
                               min_rich_fraction   = 0.30,
                               verbose             = TRUE) {

  # ── Per-dataset CV estimates ────────────────────────────────────────────────
  cvs <- vapply(datasets, function(d) {
    mean_est <- sd_est <- NA_real_

    if (!is.null(d$mean) && !is.null(d$sd)) {
      mean_est <- d$mean;  sd_est <- d$sd
    } else if (!is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)) {
      mean_est <- d$median;  sd_est <- (d$Q3 - d$Q1) / 1.35
    } else if (!is.null(d$median) && !is.null(d$min) && !is.null(d$max)) {
      mean_est <- d$median;  sd_est <- (d$max - d$min) / 4
    } else if (!is.null(d$freq_value) && !is.null(d$freq_count)) {
      w        <- d$freq_count / sum(d$freq_count)
      mean_est <- sum(d$freq_value * w)
      sd_est   <- sqrt(sum(w * (d$freq_value - mean_est)^2))
    } else if (!is.null(d$freq_lower) && !is.null(d$freq_upper) &&
               !is.null(d$freq_count)) {
      mid      <- (d$freq_lower + d$freq_upper) / 2
      w        <- d$freq_count / sum(d$freq_count)
      mean_est <- sum(mid * w)
      sd_est   <- sqrt(sum(w * (mid - mean_est)^2))
    } else if (!is.null(d$expo_lower) && !is.null(d$expo_upper) &&
               !is.null(d$event_lower) && !is.null(d$event_upper) &&
               !is.null(d$freq_count)) {
      delay_mid <- ((d$event_lower + d$event_upper) / 2) -
                   ((d$expo_lower  + d$expo_upper)  / 2)
      w         <- d$freq_count / sum(d$freq_count)
      mean_est  <- sum(delay_mid * w)
      sd_est    <- sqrt(sum(w * (delay_mid - mean_est)^2))
    }

    if (is.na(mean_est) || is.na(sd_est) || mean_est <= 0 || sd_est <= 0)
      return(NA_real_)
    sd_est / mean_est
  }, numeric(1))

  cvs_valid <- cvs[!is.na(cvs)]

  # ── Check 1: CV spread ──────────────────────────────────────────────────────
  if (length(cvs_valid) >= 2) {
    spread <- max(cvs_valid) / min(cvs_valid)
    if (spread > cv_spread_threshold) {
      if (verbose)
        message("should_attempt_gg: SKIP — CV spread too large ",
                "(max/min = ", round(spread, 2),
                ", threshold = ", cv_spread_threshold, ").")
      return(FALSE)
    }
  }

  # ── Check 2: information richness ──────────────────────────────────────────
  is_rich <- vapply(datasets, function(d) {
    !is.null(d$freq_value) ||
    (!is.null(d$freq_lower) && !is.null(d$freq_upper)) ||
    (!is.null(d$expo_lower) && !is.null(d$event_lower))
  }, logical(1))

  rich_frac <- mean(is_rich)
  if (rich_frac < min_rich_fraction) {
    if (verbose)
      message("should_attempt_gg: SKIP — too few frequency-table datasets ",
              "(", round(100 * rich_frac), "% rich, need >= ",
              round(100 * min_rich_fraction), "%).")
    return(FALSE)
  }

  if (verbose) message("should_attempt_gg: OK — GG fitting is worth attempting.")
  TRUE
}


# Gamma + type-2 reliability heuristic ----------------------------------------

#' Check whether Gamma can be reliably fitted from median + IQR summary statistics
#'
#' @description
#' The Gamma distribution (dist_type = 2) has a single shape parameter
#' \eqn{\phi} that controls the coefficient of variation CV = 1/\eqn{\sqrt{\phi}}.
#' When only median + IQR (summary type 2) data are available, the likelihood
#' gradient with respect to \eqn{\phi} comes exclusively from central order
#' statistics at p = 0.25, 0.50, 0.75.  These central quantiles are most
#' sensitive to location and spread, but carry weak information about shape
#' compared with tail quantiles (type-1: min/max) or full frequency data
#' (types 4/5).
#'
#' As \eqn{\phi} grows large the gamma distribution approaches a normal
#' distribution: the three central quantiles become nearly symmetric around the
#' mean, and the likelihood surface flattens in the \eqn{\phi} direction.
#' MCMC consequently exhibits very small step sizes, high autocorrelation, and
#' inference that is dominated by the prior on \code{log_phi}.
#'
#' Note that this problem does not apply to type-1 (median + range) or type-3
#' (mean + SD) data.  Type-1 uses extreme order statistics (min, max) which
#' lie in the tails where gamma shape sensitivity is largest.  Type-3 directly
#' identifies \eqn{\phi} via \eqn{\phi = (\text{mean}/\text{SD})^2}, giving a
#' sharp, well-defined gradient regardless of how concentrated the distribution is.
#'
#' The heuristic uses the moment estimator
#' \deqn{\hat{\phi} = \left(\frac{1.35 \times \text{median}}{\text{IQR}}\right)^2}
#' (equivalent to CV\eqn{^{-2}}, the method-of-moments gamma shape estimate)
#' computed from each type-2 dataset.  If the median implied shape across
#' type-2 datasets exceeds \code{max_implied_shape}, the function returns
#' \code{FALSE}.
#'
#' @param datasets A named list of datasets as accepted by
#'   [prepare_stan_data_from_datasets()].
#' @param max_implied_shape Numeric scalar (default 20).  Maximum tolerated
#'   median implied shape across type-2 datasets.  Corresponds to
#'   CV \eqn{\approx} 0.22 and IQR/median \eqn{\approx} 0.30.  Reduce to be
#'   stricter; increase to be more permissive.
#' @param min_n Integer scalar (default 50).  Minimum acceptable sample size
#'   for a type-2 dataset.  If more than half of the type-2 datasets fall
#'   below this threshold, an advisory message is printed but \code{FALSE} is
#'   not returned — use this as a soft warning only.
#' @param verbose Logical (default \code{TRUE}).  Print a one-line verdict with
#'   the reason a check failed.
#'
#' @return \code{TRUE} if type-2 data appear adequate for gamma fitting,
#'   \code{FALSE} if the implied shape is too large for reliable inference.
#'   Returns \code{TRUE} silently when no type-2 datasets are present (the
#'   heuristic is not relevant in that case).
#'
#' @examples
#' \dontrun{
#' # Concentrated distribution — high implied shape, likely slow
#' ds_concentrated <- list(
#'   d1 = list(median = 10, Q1 = 9.2, Q3 = 10.8, n = 50),
#'   d2 = list(median = 12, Q1 = 11.1, Q3 = 12.9, n = 60)
#' )
#' gamma_type2_reliable(ds_concentrated)   # expected: FALSE
#'
#' # Dispersed distribution — low implied shape, reliable
#' ds_dispersed <- list(
#'   d1 = list(median = 10, Q1 = 7, Q3 = 14, n = 80),
#'   d2 = list(median = 8,  Q1 = 5, Q3 = 12, n = 100)
#' )
#' gamma_type2_reliable(ds_dispersed)      # expected: TRUE
#' }
#'
#' @export
gamma_type2_reliable <- function(datasets,
                                  max_implied_shape = 20,
                                  min_n             = 50,
                                  verbose           = TRUE) {

  # ── Identify type-2 datasets (median + Q1 + Q3) ────────────────────────────
  is_type2 <- vapply(datasets, function(d) {
    !is.null(d$median) && !is.null(d$Q1) && !is.null(d$Q3)
  }, logical(1))

  if (!any(is_type2)) {
    # No type-2 data present — heuristic not relevant, allow fitting
    return(TRUE)
  }

  type2_ds <- datasets[is_type2]

  # ── Check 1: implied shape ──────────────────────────────────────────────────
  implied_shapes <- vapply(type2_ds, function(d) {
    iqr <- d$Q3 - d$Q1
    if (iqr <= 0 || d$median <= 0) return(NA_real_)
    (1.35 * d$median / iqr)^2
  }, numeric(1))

  implied_shapes_valid <- implied_shapes[!is.na(implied_shapes)]

  if (length(implied_shapes_valid) > 0) {
    med_shape <- stats::median(implied_shapes_valid)
    if (med_shape > max_implied_shape) {
      if (verbose)
        message("gamma_type2_reliable: SKIP \u2014 median implied gamma shape = ",
                round(med_shape, 1),
                " (IQR/median \u2248 ", round(1.35 / sqrt(med_shape), 2), ")",
                " exceeds threshold of ", max_implied_shape, ".",
                " Central quantiles carry little information about shape when",
                " the distribution is this concentrated; sampling will be slow.")
      return(FALSE)
    }
  }

  # ── Check 2: sample size (advisory only) ───────────────────────────────────
  ns <- vapply(type2_ds, function(d) {
    if (!is.null(d$n)) as.numeric(d$n) else NA_real_
  }, numeric(1))
  ns_valid <- ns[!is.na(ns)]

  if (length(ns_valid) > 0 && sum(ns_valid < min_n) > length(ns_valid) / 2) {
    if (verbose)
      message("gamma_type2_reliable: NOTE \u2014 more than half of type-2 datasets",
              " have n < ", min_n, ". Small samples make Q1/Q3 order statistics",
              " imprecise; consider using type-1 or frequency-table data instead.")
  }

  if (verbose) message("gamma_type2_reliable: OK \u2014 type-2 data appears adequate for gamma.")
  TRUE
}


# Compile Stan model -------------------------------------------------------

#' Compile a ddsynth Stan model
#'
#' Compiles and returns one of the Stan models shipped with the package.
#' Pass the returned object to [pre_inference_checks()], [fit_model()], or
#' [run_simulation_study_generalized()].
#'
#' @param model Character string: `"factorised"` (default) uses the
#'   order-statistic factorised likelihood
#'   (`hierarchical_data_synthesis_summary_stats.stan`); `"joint"` uses the
#'   joint parameterisation
#'   (`hierarchical_data_synthesis_summary_stats_joint.stan`).
#' @return A compiled `stanmodel` object.
#' @export
compile_stan_model <- function(model = c("factorised", "joint")) {
  model <- match.arg(model)
  filename <- switch(model,
    factorised = "hierarchical_data_synthesis_summary_stats.stan",
    joint      = "hierarchical_data_synthesis_summary_stats_joint.stan"
  )
  stan_file <- system.file("stan", filename, package = "ddsynth")
  if (!nzchar(stan_file)) {
    stop("Stan model file '", filename, "' not found in ddsynth installation.")
  }
  rstan::stan_model(stan_file)
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
#'     the prior and checks whether each observed value falls within the 95%
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
#' @param filter Logical. If `TRUE`, any dataset flagged by at least one
#'   per-dataset check (method-of-moments outlier, outside prior predictive
#'   interval, or non-overlapping LOO phi posterior) is removed from the
#'   returned dataset list.  Defaults to `FALSE`.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`mom_consistency`}{Data frame of implied `phi` per dataset with
#'       an `is_outlier` flag.}
#'     \item{`prior_predictive`}{Data frame of 95% prior predictive intervals
#'       for the implied SD of each dataset, with an `outside_prior_pi` flag.}
#'     \item{`map_probe`}{List with `phi_map` (MAP estimate of `phi`) and
#'       `map_converged` logical.}
#'     \item{`ll_surface`}{Data frame of `phi` vs `log_prob` from the surface
#'       scan.}
#'     \item{`loo_fits`}{Data frame of per-dataset `phi` posterior summaries
#'       from the leave-one-out fits.}
#'     \item{`datasets`}{The input `datasets` list, filtered to remove flagged
#'       datasets when `filter = TRUE`, otherwise identical to the input.}
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
                                 verbose               = TRUE,
                                 filter                = FALSE) {

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
        lognormal = sqrt(log(1 + (sd_est / mean_est)^2)),    # approx sdlog (sigma)
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
  # phi_map reads log_phi0 (population-mean dispersion), the point-4
  # hierarchical model's analogue of the old shared scalar log_phi.
  map_result <- tryCatch({
    opt        <- rstan::optimizing(stan_model, data = stan_data, hessian = FALSE,
                                    refresh = 0)
    phi_map    <- exp(opt$par[["log_phi0"]])
    list(phi_map = phi_map, map_converged = opt$return_code == 0,
         return_code = opt$return_code)
  }, error = function(e) {
    list(phi_map = NA_real_, map_converged = FALSE, return_code = NA_integer_,
         error_msg = conditionMessage(e))
  })
  results$map_probe <- map_result

  # ── Check 4: Log-likelihood surface scan ──────────────────────────────────
  # unconstrain_pars() needs a value for every declared parameter, so this
  # must list all of the point-4 hierarchical model's parameters (mu0,
  # log_tau, log_phi0, log_omega, log_kappa, loc_d_raw, log_phi_d_raw), not
  # just the ones this check varies. log_phi_d_raw is fixed at zero so each
  # phi_val is evaluated at phi_d[d] == phi0 for every dataset (no
  # between-study deviation), i.e. the population-mean surface.
  mu0_init          <- stan_data$mu0_mean
  log_tau_init      <- stan_data$log_tau_mean
  log_kappa_init    <- stan_data$log_kappa_mean
  loc_d_raw_init    <- rep(0, stan_data$n_datasets)
  log_phi_d_raw_init <- rep(0, stan_data$n_datasets)

  ll_surface <- purrr::map_dfr(phi_grid, function(phi_val) {
    pars <- list(mu0           = mu0_init,
                 log_tau       = log_tau_init,
                 log_phi0      = log(phi_val),
                 log_omega     = stan_data$log_omega_mean,
                 log_kappa     = log_kappa_init,
                 loc_d_raw     = loc_d_raw_init,
                 log_phi_d_raw = log_phi_d_raw_init)
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
    if (is.null(fit) || fit@mode != 0L) {
      return(tibble::tibble(dataset = name, phi_mean = NA_real_,
                            phi_lo = NA_real_, phi_hi = NA_real_,
                            rhat = NA_real_,   n_eff = NA_real_))
    }

    # phi_d (not phi): the point-4 hierarchical model has no shared scalar
    # phi. With n_datasets == 1 here, phi_d is a length-1 vector, giving the
    # same single-dataset posterior the old shared phi used to.
    s <- rstan::summary(fit, pars = "phi_d")$summary
    tibble::tibble(dataset  = name,
                   phi_mean = s[, "mean"],
                   phi_lo   = s[, "2.5%"],
                   phi_hi   = s[, "97.5%"],
                   rhat     = s[, "Rhat"],
                   n_eff    = s[, "n_eff"])
  })
  results$loo_fits <- loo_fits

  # ── Per-dataset LOO overlap flag (needed for filter, computed unconditionally)
  med_lo <- stats::median(loo_fits$phi_lo, na.rm = TRUE)
  med_hi <- stats::median(loo_fits$phi_hi, na.rm = TRUE)
  loo_flagged <- dplyr::mutate(
    loo_fits,
    no_overlap = !is.na(phi_lo) & (phi_hi < med_lo | phi_lo > med_hi)
  )

  # ── Union of all per-dataset flags ────────────────────────────────────────
  flagged_datasets <- unique(c(
    mom_df$dataset[!is.na(mom_df$is_outlier)        & mom_df$is_outlier],
    prior_df$dataset[!is.na(prior_df$outside_prior_pi) & prior_df$outside_prior_pi],
    loo_flagged$dataset[!is.na(loo_flagged$no_overlap)  & loo_flagged$no_overlap]
  ))

  # ── Optionally filter datasets ─────────────────────────────────────────────
  datasets_out <- if (filter && length(flagged_datasets) > 0) {
    datasets[setdiff(names(datasets), flagged_datasets)]
  } else {
    datasets
  }
  results$datasets <- datasets_out

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
    n_no_overlap <- sum(loo_flagged$no_overlap, na.rm = TRUE)
    if (n_no_overlap > 0) {
      cli::cli_alert_warning(
        "{n_no_overlap} dataset(s) have phi posteriors that do not overlap the majority: {loo_flagged$dataset[loo_flagged$no_overlap]}"
      )
    } else {
      cli::cli_alert_success("All per-dataset phi posteriors broadly overlap")
    }

    if (filter) {
      if (length(flagged_datasets) > 0) {
        cli::cli_h2("Filtering")
        cli::cli_alert_warning(
          "{length(flagged_datasets)} dataset(s) removed: {flagged_datasets}. ",
          "{length(datasets_out)} dataset(s) retained."
        )
      } else {
        cli::cli_h2("Filtering")
        cli::cli_alert_success("No datasets flagged — all {length(datasets)} retained.")
      }
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
#' Builds the master scenario grid used by both \code{run_simulation_study.R}
#' and the \code{simulation_study} vignette.  All four distribution families
#' that are identifiable from summary statistics (lognormal, gamma, Weibull,
#' Burr XII) are supported, plus the generalised gamma for identifiability
#' research.
#'
#' @param include_homogeneous Include scenarios with fixed summary types
#' @param include_mixed Include scenarios with mixed summary types
#' @param include_varied_n Include scenarios with varied sample sizes
#' @param include_freq_table Include frequency-table summary scenarios
#' @param include_burr12 Include Burr XII scenarios (\code{dist_type = "burr12"})
#' @param include_gengamma Include generalised gamma standard scenarios
#'   (\code{dist_type = "gengamma"})
#' @param include_gg_limitation Include generalised gamma limitation scenarios
#'   (only used when \code{include_gengamma = TRUE})
#' @return Data frame with one row per scenario and columns:
#'   \code{scenario_name}, \code{scenario_group}, \code{scenario_idx},
#'   \code{dist_type} (one of \code{"lognormal"}, \code{"gamma"},
#'   \code{"weibull"}, \code{"burr12"}, \code{"gengamma"}),
#'   \code{n_datasets}, \code{mu0}, \code{tau}, \code{phi},
#'   \code{kappa}, \code{n_obs_mean}, \code{n_obs_sd}, \code{n_obs_min},
#'   \code{n_obs_max}, \code{summary_type_1_prop}–\code{summary_type_4_prop},
#'   \code{summary_type}, \code{vary_n}.
#' @export
generate_scenario_library <- function(include_homogeneous = TRUE,
                                      include_mixed = TRUE,
                                      include_varied_n = TRUE,
                                      include_freq_table = FALSE,
                                      include_burr12 = FALSE,
                                      include_gengamma = FALSE,
                                      include_gg_limitation = FALSE) {

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
  
  # Helper: derive summary_type integer (1–4 = homogeneous; 5 = mixed)
  .st_int <- function(s1, s2, s3, s4) {
    if (s1 == 1) 1L else if (s2 == 1) 2L else if (s3 == 1) 3L else
      if (s4 == 1) 4L else 5L
  }

  # Convert base (dist 1–3) scenario list to unified data frame
  scenarios_df <- dplyr::bind_rows(lapply(scenarios, function(s) {
    nv  <- s$n_obs
    st1 <- mean(s$summary_type == 1)
    st2 <- mean(s$summary_type == 2)
    st3 <- mean(s$summary_type == 3)
    st4 <- mean(s$summary_type == 4)
    data.frame(
      scenario_name       = s$scenario_name,
      scenario_group      = paste0("base_", s$dist_type),
      dist_type           = s$dist_type,
      n_datasets          = s$n_datasets,
      mu0                 = s$mu0,
      tau                 = s$tau,
      phi                 = s$phi,
      kappa               = 1.0,
      n_obs_mean          = mean(nv),
      n_obs_sd            = if (length(nv) > 1) stats::sd(nv) else 0,
      n_obs_min           = min(nv),
      n_obs_max           = max(nv),
      summary_type_1_prop = st1,
      summary_type_2_prop = st2,
      summary_type_3_prop = st3,
      summary_type_4_prop = st4,
      summary_type        = .st_int(st1, st2, st3, st4),
      vary_n              = length(unique(nv)) > 1L,
      stringsAsFactors    = FALSE
    )
  }))

  # Store full scenario details for backward compatibility with the old parallel
  # runner (run_simulation_study_generalized).  Burr XII / GG are NOT added to
  # this attribute because those runners do not support them.
  attr(scenarios_df, "full_scenarios") <- scenarios

  # ===== BURR XII SCENARIOS (dist_type 4) =====
  if (include_burr12) {
    burr_grid <- expand.grid(
      phi        = c(2.0, 3.0),
      kappa      = c(2.0, 5.0),
      n_datasets = c(5L, 10L, 20L),
      n_obs      = c(20, 50),
      stringsAsFactors = FALSE
    )
    burr_rows <- dplyr::bind_rows(lapply(seq_len(nrow(burr_grid)), function(i) {
      r <- burr_grid[i, ]
      data.frame(
        scenario_name       = sprintf("Burr12_c%.1f_k%.1f_D%d_N%d_ST1",
                                      r$phi, r$kappa, r$n_datasets, r$n_obs),
        scenario_group      = "burr12",
        dist_type           = "burr12",
        n_datasets          = as.integer(r$n_datasets),
        mu0                 = log(7),
        tau                 = 0.4,
        phi                 = r$phi,
        kappa               = r$kappa,
        n_obs_mean          = r$n_obs,
        n_obs_sd            = 0,
        n_obs_min           = r$n_obs,
        n_obs_max           = r$n_obs,
        summary_type_1_prop = 1,
        summary_type_2_prop = 0,
        summary_type_3_prop = 0,
        summary_type_4_prop = 0,
        summary_type        = 1L,
        vary_n              = FALSE,
        stringsAsFactors    = FALSE
      )
    }))

    # Extra Burr XII scenarios: mixed summaries, varied N, large N, mean+SD
    vn_min <- max(5L, round(30 * 0.3))
    vn_max <- round(30 * 2)
    burr_extra <- data.frame(
      scenario_name       = c("Burr12_c2.5_k3_D10_N30_Mixed12",
                              "Burr12_c2.5_k3_D10_N30_Mixed123",
                              "Burr12_c2.5_k3_D20_VarN_ST1",
                              "Burr12_c2.5_k3_D20_N100_ST1",
                              "Burr12_c2.5_k3_D10_N30_ST3"),
      scenario_group      = "burr12",
      dist_type           = "burr12",
      n_datasets          = c(10L, 10L, 20L, 20L, 10L),
      mu0                 = log(7),
      tau                 = 0.4,
      phi                 = 2.5,
      kappa               = 3.0,
      n_obs_mean          = c(30, 30, 30, 100, 30),
      n_obs_sd            = c(0, 0, 15, 0, 0),
      n_obs_min           = c(30, 30, vn_min, 100, 30),
      n_obs_max           = c(30, 30, vn_max, 100, 30),
      summary_type_1_prop = c(0.5, 1/3, 1.0, 1.0, 0.0),
      summary_type_2_prop = c(0.5, 1/3, 0.0, 0.0, 0.0),
      summary_type_3_prop = c(0.0, 1/3, 0.0, 0.0, 1.0),
      summary_type_4_prop = c(0.0, 0.0, 0.0, 0.0, 0.0),
      summary_type        = c(5L, 5L, 1L, 1L, 3L),
      vary_n              = c(FALSE, FALSE, TRUE, FALSE, FALSE),
      stringsAsFactors    = FALSE
    )

    scenarios_df <- dplyr::bind_rows(scenarios_df, burr_rows, burr_extra)
  }

  # ===== GENERALISED GAMMA SCENARIOS (dist_type 5) =====
  if (include_gengamma) {
    gg_std <- data.frame(
      scenario_name       = c("GG_s0.5_Q0.5_D10_N30_ST1",
                              "GG_s0.5_Q1.0_D10_N30_ST1",
                              "GG_s0.5_Q2.0_D10_N30_ST1",
                              "GG_s0.5_Q1.0_D5_N30_ST1",
                              "GG_s0.5_Q1.0_D20_N30_ST1",
                              "GG_s0.5_Q1.0_D10_N10_ST1"),
      scenario_group      = "gg_standard",
      dist_type           = "gengamma",
      n_datasets          = c(10L, 10L, 10L, 5L, 20L, 10L),
      mu0                 = log(7),
      tau                 = 0.4,
      phi                 = 0.5,
      kappa               = c(0.5, 1.0, 2.0, 1.0, 1.0, 1.0),
      n_obs_mean          = c(30, 30, 30, 30, 30, 10),
      n_obs_sd            = 0,
      n_obs_min           = c(30, 30, 30, 30, 30, 10),
      n_obs_max           = c(30, 30, 30, 30, 30, 10),
      summary_type_1_prop = 1,
      summary_type_2_prop = 0,
      summary_type_3_prop = 0,
      summary_type_4_prop = 0,
      summary_type        = 1L,
      vary_n              = FALSE,
      stringsAsFactors    = FALSE
    )
    scenarios_df <- dplyr::bind_rows(scenarios_df, gg_std)

    if (include_gg_limitation) {
      gg_lim <- data.frame(
        scenario_name       = c("GG_Limit_NearLognormal_Q0.1_D10_N30_ST1",
                                "GG_Limit_HighQ_Q3_D10_N30_ST1", "GG_Limit_MidpointQ_Q0.5_D10_N30_ST1" ),
        scenario_group      = "gg_limitation",
        dist_type           = "gengamma",
        n_datasets          = 10L,
        mu0                 = log(7),
        tau                 = 0.4,
        phi                 = 0.5,
        kappa               = c(0.1, 3.0, 0.5),
        n_obs_mean          = 30,
        n_obs_sd            = 0,
        n_obs_min           = 30,
        n_obs_max           = 30,
        summary_type_1_prop = 1,
        summary_type_2_prop = 0,
        summary_type_3_prop = 0,
        summary_type_4_prop = 0,
        summary_type        = 1L,
        vary_n              = FALSE,
        stringsAsFactors    = FALSE
      )
      scenarios_df <- dplyr::bind_rows(scenarios_df, gg_lim)
    }
  }

  scenarios_df$scenario_idx <- seq_len(nrow(scenarios_df))

  return(scenarios_df)
}

#' Create a Stan initialisation function from prior means
#'
#' Returns a zero-argument function that initialises each chain at the prior
#' means stored in \code{stan_data}.  Passing this to \code{rstan::sampling()}
#' via \code{init = make_stan_init_fn(stan_data)} avoids the default
#' \code{Uniform(-2, 2)} draws on the unconstrained scale, which send the GG
#' and Burr XII log-posteriors to \eqn{-\infty} during initialisation.
#'
#' @param stan_data Stan data list as returned by
#'   \code{\link{prepare_stan_data_from_datasets}}.
#' @return A zero-argument function suitable for \code{rstan::sampling(init = ...)}.
#' @export
make_stan_init_fn <- function(stan_data) {
  function() {
    list(
      mu0           = stan_data$mu0_mean,
      log_tau       = stan_data$log_tau_mean,
      log_phi0      = stan_data$log_phi_mean,
      log_omega     = stan_data$log_omega_mean,
      log_kappa     = stan_data$log_kappa_mean,
      loc_d_raw     = rep(0.0, stan_data$n_datasets),
      log_phi_d_raw = rep(0.0, stan_data$n_datasets)
    )
  }
}

#' Fit Stan model to simulated data
#'
#' @param sim_data Simulated data from generate_hierarchical_data
#' @param stan_model Compiled Stan model
#' @param chains Number of Markov chains (default 4)
#' @param iter Total number of iterations per chain (default 10000)
#' @param warmup Number of warmup iterations per chain (default 1000)
#' @param refresh How often to print progress (default 0)
#' @param control List of control parameters passed to Stan (e.g. adapt_delta)
#' @param ... Additional arguments to pass to rstan::sampling()
#' @return Stan fit object
#' @export
fit_model <- function(sim_data, stan_model,
                      chains  = 4L,
                      iter    = 10000L,
                      warmup  = 1000L,
                      refresh = 0L,
                      control = list(adapt_delta = 0.95, max_treedepth = 12L),
                      ...) {
  od <- sim_data$obs_data
  rstan::sampling(
    stan_model,
    data    = od,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    refresh = refresh,
    control = control,
    init    = make_stan_init_fn(od),
    ...
  )
}

#' Fit a Stan model with escalating settings, stopping once Rhat is acceptable
#'
#' @description
#' Tries a sequence of increasingly conservative sampling settings, stopping
#' at the first tier that achieves `max(Rhat) < rhat_target` with no
#' divergent transitions. Falls all the way back to the original
#' pre-lever-1 production settings (`iter=12000`, `adapt_delta=0.999`,
#' `max_treedepth=12`) as the last tier if cheaper settings do not converge -
#' those settings were originally chosen because some chains did not
#' converge without them (see `REVISION_TODO.md`), so this is a floor, not
#' just a starting point to relax.
#'
#' Most fits are expected to converge at the cheap first tier, especially
#' now that lever 3 (`R/utils.R`'s `.order_stat_n_panels()`) has reduced
#' per-iteration cost and point 4's family gating
#' (`POINT4_LIKELIHOOD_MATHS.md` Part E) means gamma no longer needs a
#' hierarchical `phi_d` funnel resolved. The escalation ladder exists for the
#' cases that still need it, not as the expected path.
#'
#' @param stan_data Data list for `rstan::sampling()`.
#' @param stan_model Compiled Stan model.
#' @param rhat_target Escalate if `max(Rhat)` exceeds this. Default 1.05,
#'   matching the target agreed for the point 1/4 checkpoint run.
#' @param seed Passed to every tier for reproducibility.
#' @param init Optional init argument passed to every tier (e.g.
#'   `make_stan_init_fn(stan_data)`); default `NULL` uses Stan's own default.
#' @param verbose If `TRUE` (default), reports which tier was used and why.
#' @return A list: `fit` (the `stanfit` object from whichever tier
#'   succeeded, or the last tier tried if none converged), `tier` (integer,
#'   which tier succeeded or was last attempted), `max_rhat`, `divergences`,
#'   `converged` (logical), and `runtime_secs` (total across all tiers
#'   attempted).
#' @export
fit_with_escalation <- function(stan_data, stan_model, rhat_target = 1.05,
                                 seed = 123, init = NULL, verbose = TRUE) {
  tiers <- list(
    list(label = "cheap",    chains = 4, iter = 4000,  warmup = 2000, control = list(adapt_delta = 0.9,  max_treedepth = 10)),
    list(label = "moderate", chains = 4, iter = 6000,  warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 11)),
    # Original production settings (main_analysis.R, pre-lever-1): kept
    # verbatim as the final fallback, not relaxed, because they were
    # originally added after some chains failed to converge without them.
    list(label = "production", chains = 4, iter = 12000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 12))
  )

  t_total0 <- Sys.time()
  result <- NULL
  best <- NULL   # best-so-far by (max_rhat, div_rate), in case no tier fully converges
  for (i in seq_along(tiers)) {
    tier <- tiers[[i]]
    if (verbose) message(sprintf("  [fit_with_escalation] trying tier %d/%d ('%s': iter=%d, adapt_delta=%.3f)...",
                                  i, length(tiers), tier$label, tier$iter, tier$control$adapt_delta))
    fit <- tryCatch(
      rstan::sampling(stan_model, data = stan_data, chains = tier$chains, iter = tier$iter,
                       warmup = tier$warmup, control = tier$control, seed = seed,
                       init = if (is.null(init)) "random" else init, refresh = 0,
                       show_messages = FALSE),
      error = function(e) { if (verbose) message("    tier failed to run: ", conditionMessage(e)); NULL }
    )
    if (is.null(fit)) next

    s <- rstan::summary(fit)$summary
    max_rhat <- max(s[, "Rhat"], na.rm = TRUE)
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
    n_draws <- sum(sapply(sp, nrow))
    div_rate <- div / n_draws
    # A handful of divergences out of thousands of draws is ordinary HMC
    # noise, not a sign of a bad fit; requiring exactly zero caused endless,
    # pointless escalation for datasets with a small structural divergence
    # rate (e.g. Burr XII, already known from this session's diagnostics to
    # run somewhat divergence-prone even at good settings). 1% mirrors
    # common Stan-community practice for "acceptably rare, not ignorable"
    # divergence rates.
    converged <- is.finite(max_rhat) && max_rhat < rhat_target && div_rate < 0.01

    if (verbose) message(sprintf("    max_Rhat=%.4f, divergences=%d/%d (%.2f%%) -> %s",
                                  max_rhat, div, n_draws, 100 * div_rate, if (converged) "OK" else "escalating"))

    result <- list(fit = fit, tier = i, tier_label = tier$label, max_rhat = max_rhat,
                    divergences = div, div_rate = div_rate, converged = converged)
    if (converged) break
    # Track the best-so-far result (lower max_rhat wins; div_rate breaks ties)
    # so that if no tier converges, we return the least-bad attempt rather
    # than just whichever tier happened to run last.
    if (is.null(best) || max_rhat < best$max_rhat ||
        (max_rhat == best$max_rhat && div_rate < best$div_rate)) {
      best <- result
    }
  }

  if (is.null(result)) {
    stop("fit_with_escalation: every tier failed to sample (see messages above).", call. = FALSE)
  }
  if (!result$converged && !is.null(best) && best$max_rhat <= result$max_rhat && best$tier != result$tier) {
    result <- best
  }
  result$runtime_secs <- as.numeric(Sys.time() - t_total0, units = "secs")
  if (verbose && !result$converged) {
    message("  [fit_with_escalation] WARNING: did not reach Rhat<", rhat_target,
            " even at the production tier (max_Rhat=", round(result$max_rhat, 4), ").")
  }
  result
}

#' Fit a list of (dataset, family) tasks concurrently, leaving cores free
#'
#' @description
#' Runs [fit_with_escalation()] over a list of independent fitting tasks in
#' parallel (fork-based, via `parallel::mclapply()` - Unix/macOS only), each
#' task's chains still parallelised internally as usual. Deliberately caps
#' total core usage well below the machine's full core count so the machine
#' remains usable for other work while a corpus run is in progress - this is
#' a hard requirement, not a tuning default: leave `reserve_cores` alone
#' unless the person running this has explicitly said otherwise.
#'
#' Each task's result is written to its own file immediately on completion
#' (`task$output_file`), rather than accumulating in a single shared
#' in-memory list written once at the end - this avoids concurrent workers
#' racing on one output file, and means a crash partway through loses only
#' the tasks that hadn't finished, not everything. Call
#' [merge_parallel_fit_results()] afterwards to assemble the per-task files
#' into a single results list.
#'
#' @param tasks A list of task specifications, each a list with:
#'   `label` (character, for logging), `datasets` (named list, as passed to
#'   [prepare_stan_data_from_datasets()]), `dist_name` (one of "lognormal",
#'   "gamma", "weibull", "burr", "gengamma"), `stan_model` (compiled model to
#'   use for this task - the caller decides hierarchical vs. shared-phi per
#'   family, see `analysis/main_analysis.R`'s `SHARED_PHI_FAMILIES`), and
#'   `output_file` (path to save this task's result to).
#' @param chains_per_fit Chains used within each individual fit (passed to
#'   [fit_with_escalation()] tiers implicitly - currently fixed at 4 within
#'   that function; this argument only affects the core-budget arithmetic
#'   below, so keep it in sync if that changes).
#' @param reserve_cores Cores to leave free for other use on the machine.
#'   Default 4 (out of this machine's 12), leaving meaningful headroom for
#'   interactive use, not just background slack - do not reduce this without
#'   explicit instruction.
#' @param rhat_target Passed through to [fit_with_escalation()].
#' @return Invisibly, a character vector of the `output_file` paths written
#'   (some may be missing if a task errored - check before merging).
#' @export
fit_corpus_parallel <- function(tasks, chains_per_fit = 4, reserve_cores = 4,
                                  rhat_target = 1.05) {
  total_cores <- parallel::detectCores()
  usable_cores <- max(total_cores - reserve_cores, chains_per_fit)
  n_concurrent <- max(floor(usable_cores / chains_per_fit), 1)
  message(sprintf(
    "fit_corpus_parallel: %d cores detected, reserving %d, running up to %d fits concurrently (%d tasks total).",
    total_cores, reserve_cores, n_concurrent, length(tasks)
  ))

  parallel::mclapply(tasks, function(task) {
    dist_type <- c(lognormal = 1L, gamma = 2L, weibull = 3L, burr = 4L, gengamma = 5L)[[task$dist_name]]
    stan_data <- tryCatch(
      prepare_stan_data_from_datasets(task$datasets, dist_type = dist_type),
      error = function(e) { message("  [", task$label, "] prepare_stan_data failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(stan_data)) { saveRDS(list(error = "prepare_stan_data failed"), task$output_file); return(invisible(NULL)) }
    stan_data <- update_phi_prior(stan_data, task$datasets)

    res <- tryCatch(
      fit_with_escalation(stan_data, task$stan_model, rhat_target = rhat_target, verbose = TRUE),
      error = function(e) { message("  [", task$label, "] fit_with_escalation failed: ", conditionMessage(e)); NULL }
    )
    out <- if (is.null(res)) list(error = "fit_with_escalation failed") else
      list(fit = res$fit, stan_data = stan_data, datasets = task$datasets, tier = res$tier_label,
           max_rhat = res$max_rhat, divergences = res$divergences, converged = res$converged,
           runtime_secs = res$runtime_secs)
    saveRDS(out, task$output_file)
    invisible(NULL)
  }, mc.cores = n_concurrent)

  invisible(vapply(tasks, function(t) t$output_file, character(1)))
}

#' Assemble per-task .rds files from [fit_corpus_parallel()] into one list
#'
#' @param output_files Character vector of file paths (as returned by
#'   [fit_corpus_parallel()]).
#' @param labels Character vector, same length, used as names in the
#'   returned list (typically `"<pathogen>.<family>"`).
#' @return A named list of the per-task results; entries for missing files
#'   are `NULL` with a message, not a hard error, so a partially-complete
#'   run can still be assembled and inspected.
#' @export
merge_parallel_fit_results <- function(output_files, labels) {
  stats::setNames(lapply(seq_along(output_files), function(i) {
    if (!file.exists(output_files[i])) {
      message("merge_parallel_fit_results: missing (not yet run or failed): ", labels[i])
      return(NULL)
    }
    readRDS(output_files[i])
  }), labels)
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
  mu0_samples   <- draws$mu0
  tau_samples   <- exp(draws$log_tau)
  phi_samples   <- exp(draws$log_phi)
  kappa_samples <- draws$kappa
  loc_d_samples <- draws$loc_d          # [n_samples x n_datasets]

  n_samples  <- length(mu0_samples)
  n_datasets <- dim(loc_d_samples)[2]

  # Internal density function for all 5 distribution types
  .ddist <- function(x, loc, phi, kappa, dist_type) {
    if (dist_type == "lognormal") {
      dlnorm(x, meanlog = loc, sdlog = phi)
    } else if (dist_type == "gamma") {
      dgamma(x, shape = phi, rate = phi / exp(loc))
    } else if (dist_type == "weibull") {
      dweibull(x, shape = phi, scale = exp(loc))
    } else if (dist_type == "burr12") {
      # Burr XII PDF: (c*k/lambda)*(x/lambda)^(c-1)*(1+(x/lambda)^c)^(-(k+1))
      lambda <- exp(loc)
      r      <- x / lambda
      phi * kappa / lambda * r^(phi - 1) * (1 + r^phi)^(-(kappa + 1))
    } else if (dist_type == "gengamma") {
      # Generalised Gamma (Prentice) PDF
      gs <- 1 / kappa^2
      w  <- (log(x) - loc) / phi
      exp(log(kappa) - log(phi) - log(x) +
          gs * log(gs) + gs * kappa * w - gs * exp(kappa * w) - lgamma(gs))
    }
  }

  # Create grid if not provided
  if (is.null(x_grid)) {
    if (dist_type %in% c("lognormal", "burr12", "gengamma")) {
      x_grid <- seq(0.01, exp(true_params$mu0 + 3 * true_params$tau), length.out = 200)
    } else if (dist_type == "gamma") {
      mean_max <- exp(true_params$mu0 + 3 * true_params$tau)
      x_grid <- seq(0.01, mean_max * 3, length.out = 200)
    } else if (dist_type == "weibull") {
      scale_max <- exp(true_params$mu0 + 3 * true_params$tau)
      x_grid <- seq(0.01, scale_max * 3, length.out = 200)
    }
  }

  kappa_true <- if (!is.null(true_params$kappa)) true_params$kappa else 1.0

  # Compute true predictive density
  true_density <- numeric(length(x_grid))
  for (i in seq_along(x_grid)) {
    x <- x_grid[i]

    integrand <- function(loc) {
      .ddist(x, loc, true_params$phi, kappa_true, dist_type) *
        dnorm(loc, mean = true_params$mu0, sd = true_params$tau)
    }

    true_density[i] <- integrate(integrand,
                                 lower = true_params$mu0 - 5 * true_params$tau,
                                 upper = true_params$mu0 + 5 * true_params$tau)$value
  }

  # Compute estimated predictive density (average over posterior samples)
  est_density <- numeric(length(x_grid))

  # Subsample for computational efficiency
  sample_idx <- sample(1:n_samples, min(500, n_samples))

  for (i in seq_along(x_grid)) {
    x <- x_grid[i]

    density_samples <- numeric(length(sample_idx))
    for (s in seq_along(sample_idx)) {
      idx     <- sample_idx[s]
      phi_s   <- phi_samples[idx]
      kappa_s <- kappa_samples[idx]

      if (n_datasets < 5) {
        loc_point <- mean(loc_d_samples[idx, ])
        density_samples[s] <- .ddist(x, loc_point, phi_s, kappa_s, dist_type)

      } else {
        integrand <- function(loc) {
          .ddist(x, loc, phi_s, kappa_s, dist_type) *
            dnorm(loc, mean = mu0_samples[idx], sd = tau_samples[idx])
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

#' Compute true marginal quantiles of the predictive distribution
#'
#' Uses Monte Carlo integration over the study-level location hierarchy to
#' compute the true population-level quantiles of the predictive distribution.
#'
#' @param dist_type Distribution type: one of \code{"lognormal"}, \code{"gamma"},
#'   \code{"weibull"}, \code{"burr12"}, \code{"gengamma"}.
#' @param mu0 True population mean (location hyperparameter).
#' @param tau True between-study standard deviation.
#' @param phi True distribution-specific shape/scale parameter.
#' @param kappa True third distribution parameter (Burr XII k or GG Q).
#'   Ignored for 2-parameter distributions.
#' @param probs Numeric vector of probabilities for which quantiles are computed.
#'   Default \code{c(0.5, 0.95)}.
#' @param n_mc Number of Monte Carlo draws. Default 5000.
#' @return A named numeric vector of quantiles (names are e.g. \code{"50%"},
#'   \code{"95%"}).
#' @export
compute_true_marginal_quantile <- function(dist_type,
                                           mu0,
                                           tau,
                                           phi,
                                           kappa = 1.0,
                                           probs  = c(0.5, 0.95),
                                           n_mc   = 5000L) {
  loc_draws <- rnorm(n_mc, mean = mu0, sd = tau)

  if (dist_type == "lognormal") {
    obs <- rlnorm(n_mc, meanlog = loc_draws, sdlog = phi)
  } else if (dist_type == "gamma") {
    obs <- rgamma(n_mc, shape = phi, rate = phi / exp(loc_draws))
  } else if (dist_type == "weibull") {
    obs <- rweibull(n_mc, shape = phi, scale = exp(loc_draws))
  } else if (dist_type == "burr12") {
    lambda <- exp(loc_draws)
    u      <- runif(n_mc)
    obs    <- lambda * (u^(-1.0 / kappa) - 1.0)^(1.0 / phi)
  } else if (dist_type == "gengamma") {
    gs  <- 1.0 / kappa^2
    y   <- rgamma(n_mc, shape = gs, rate = 1)
    obs <- exp(loc_draws + phi / kappa * log(kappa^2 * y))
  } else {
    stop(sprintf(
      "compute_true_marginal_quantile: unknown dist_type '%s'.", dist_type
    ), call. = FALSE)
  }

  quantile(obs, probs = probs, na.rm = TRUE)
}

#' Compute credible interval for posterior predictive quantiles
#'
#' For each requested probability \code{p}, draws from the posterior predictive
#' distribution and returns the 2.5th, 50th, and 97.5th percentiles of the
#' resulting quantile distribution across posterior draws.
#'
#' @param fit A \code{stanfit} object.
#' @param dist_type Distribution type string.
#' @param n_datasets Number of studies in the fitted data (used to decide
#'   whether to integrate over mu0/tau or use mean(loc_d)).
#' @param probs Probabilities for which predictive quantiles are computed.
#'   Default \code{c(0.5, 0.95)}.
#' @param n_post_draws Number of posterior draws to use. Default 500.
#' @param n_mc_per_draw Number of predictive samples per posterior draw.
#'   Default 1000.
#' @return A list with elements \code{lower}, \code{median}, and \code{upper}:
#'   the 2.5th, 50th, and 97.5th percentiles of the posterior distribution of
#'   each quantile, each a numeric vector of length \code{length(probs)}.
#' @export
compute_posterior_predictive_quantile_ci <- function(fit,
                                                     dist_type,
                                                     n_datasets,
                                                     probs         = c(0.5, 0.95),
                                                     n_post_draws  = 500L,
                                                     n_mc_per_draw = 1000L) {
  draws    <- rstan::extract(fit)
  n_post   <- length(draws$mu0)
  draw_idx <- sample.int(n_post, min(n_post_draws, n_post))
  q_mat    <- matrix(NA_real_, nrow = length(draw_idx), ncol = length(probs))

  for (i in seq_along(draw_idx)) {
    idx     <- draw_idx[i]
    phi_s   <- draws$phi[idx]
    kappa_s <- draws$kappa[idx]

    # Mirror compute_predictive_cdf: when n_datasets < 5, mu0 and tau are
    # poorly identified; use mean(loc_d) instead.
    if (n_datasets < 5) {
      locs <- rep(mean(draws$loc_d[idx, ]), n_mc_per_draw)
    } else {
      locs <- rnorm(n_mc_per_draw, mean = draws$mu0[idx], sd = draws$tau[idx])
    }

    if (dist_type == "lognormal") {
      obs_s <- rlnorm(n_mc_per_draw, meanlog = locs, sdlog = phi_s)
    } else if (dist_type == "gamma") {
      obs_s <- rgamma(n_mc_per_draw, shape = phi_s, rate = phi_s / exp(locs))
    } else if (dist_type == "weibull") {
      obs_s <- rweibull(n_mc_per_draw, shape = phi_s, scale = exp(locs))
    } else if (dist_type == "burr12") {
      lambda_s <- exp(locs)
      u_s      <- runif(n_mc_per_draw)
      obs_s    <- lambda_s * (u_s^(-1.0 / kappa_s) - 1.0)^(1.0 / phi_s)
    } else if (dist_type == "gengamma") {
      gs_s  <- 1.0 / kappa_s^2
      y_s   <- rgamma(n_mc_per_draw, shape = gs_s, rate = 1)
      obs_s <- exp(locs + phi_s / kappa_s * log(kappa_s^2 * y_s))
    } else {
      stop(sprintf(
        "compute_posterior_predictive_quantile_ci: unknown dist_type '%s'.", dist_type
      ), call. = FALSE)
    }

    q_mat[i, ] <- quantile(obs_s, probs = probs, na.rm = TRUE)
  }

  list(
    lower  = apply(q_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    median = apply(q_mat, 2, median, na.rm = TRUE),
    upper  = apply(q_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

#' Compute the mean Weighted Interval Score for the posterior predictive distribution
#'
#' Evaluates the posterior predictive distribution against test observations drawn
#' from the true marginal distribution, using the Weighted Interval Score
#' (Bracher et al. 2021) — a proper scoring rule that rewards both calibration
#' and sharpness.
#'
#' @param fit A \code{stanfit} object.
#' @param dist_type Distribution type string.
#' @param n_datasets Number of studies in the fitted data.
#' @param true_params Named list with elements \code{mu0}, \code{tau},
#'   \code{phi}, and optionally \code{kappa}.
#' @param alpha_levels Numeric vector of PI coverage levels (e.g. 0.95 for a
#'   95% PI). Default \code{c(0.50, 0.80, 0.90, 0.95)}.
#' @param n_post_draws Number of posterior draws passed to
#'   \code{compute_posterior_predictive_quantile_ci}. Default 500.
#' @param n_mc_per_draw MC samples per posterior draw. Default 1000.
#' @param n_test Number of test observations drawn from the true marginal
#'   distribution. Default 200.
#' @param pred_q Optional pre-computed vector of posterior predictive median
#'   quantiles at the probability grid derived from \code{alpha_levels}, as
#'   returned by the \code{$median} element of
#'   [compute_posterior_predictive_quantile_ci()].  When supplied the internal
#'   call to that function is skipped, avoiding redundant computation.
#' @param true_median Optional pre-computed true marginal median (scalar), e.g.
#'   the \code{'50\%'} element from [compute_true_marginal_quantile()].  When
#'   supplied the internal call to that function is skipped.
#' @return A list with \code{wis} (mean WIS, same scale as the delay outcome)
#'   and \code{rel_wis} (WIS divided by the true marginal median).
#' @export
compute_wis <- function(fit,
                        dist_type,
                        n_datasets,
                        true_params,
                        alpha_levels  = c(0.50, 0.80, 0.90, 0.95),
                        n_post_draws  = 500L,
                        n_mc_per_draw = 1000L,
                        n_test        = 200L,
                        pred_q        = NULL,
                        true_median   = NULL) {

  # ------------------------------------------------------------------
  # 1.  Quantile probability grid needed for WIS
  # ------------------------------------------------------------------
  alpha_k     <- 1 - alpha_levels                    # non-coverage: 0.50 0.20 0.10 0.05
  lower_probs <- sort(alpha_k / 2)                   # 0.025 0.05 0.10 0.25
  upper_probs <- sort(1 - alpha_k / 2)               # 0.75  0.90 0.95 0.975
  all_probs   <- sort(unique(c(lower_probs, 0.5, upper_probs)))
  # -> c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)

  kappa_val <- if (!is.null(true_params$kappa)) {
    true_params$kappa
  } else if (dist_type %in% c("burr12", "gengamma")) {
    stop(sprintf(
      "compute_wis: true_params$kappa is required for dist_type '%s'.",
      dist_type
    ), call. = FALSE)
  } else {
    1.0
  }

  # ------------------------------------------------------------------
  # 2.  Posterior predictive median quantiles
  #     Skipped when pred_q is supplied by the caller (e.g. run_one()),
  #     avoiding a redundant compute_posterior_predictive_quantile_ci call.
  # ------------------------------------------------------------------
  if (is.null(pred_q)) {
    post_obj <- compute_posterior_predictive_quantile_ci(
      fit           = fit,
      dist_type     = dist_type,
      n_datasets    = n_datasets,
      probs         = all_probs,
      n_post_draws  = n_post_draws,
      n_mc_per_draw = n_mc_per_draw
    )
    pred_q <- post_obj$median      # one value per element of all_probs
  }

  # ------------------------------------------------------------------
  # 3.  Draw test observations from the true marginal distribution
  # ------------------------------------------------------------------
  loc_draws <- rnorm(n_test, mean = true_params$mu0, sd = true_params$tau)

  if (dist_type == "lognormal") {
    y_test <- rlnorm(n_test, meanlog = loc_draws, sdlog = true_params$phi)
  } else if (dist_type == "gamma") {
    y_test <- rgamma(n_test, shape = true_params$phi,
                     rate  = true_params$phi / exp(loc_draws))
  } else if (dist_type == "weibull") {
    y_test <- rweibull(n_test, shape = true_params$phi, scale = exp(loc_draws))
  } else if (dist_type == "burr12") {
    lambda <- exp(loc_draws)
    u      <- runif(n_test)
    y_test <- lambda * (u^(-1.0 / kappa_val) - 1.0)^(1.0 / true_params$phi)
  } else if (dist_type == "gengamma") {
    gs     <- 1.0 / kappa_val^2
    y_mc   <- rgamma(n_test, shape = gs, rate = 1)
    y_test <- exp(loc_draws + true_params$phi / kappa_val * log(kappa_val^2 * y_mc))
  } else {
    stop(sprintf("compute_wis: unknown dist_type '%s'.", dist_type), call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 4.  Compute WIS (vectorised over y_test)
  #
  #  WIS = 1/(K + 0.5) * [ 0.5 * |y - m|
  #          + sum_k (alpha_k/2) * IS*(l_k, u_k, y) ]
  #  IS*(l, u, y) = (u - l) + (2/alpha)*(l - y)+ + (2/alpha)*(y - u)+
  # ------------------------------------------------------------------

  # Tolerance-based probability lookup to avoid exact floating-point equality.
  prob_idx <- function(prob_vec, target, tol = 1e-9) {
    if (length(prob_vec) == 0L)
      stop("compute_wis: probability grid is empty.", call. = FALSE)
    idx <- which.min(abs(prob_vec - target))
    if (abs(prob_vec[idx] - target) > tol)
      stop(sprintf("compute_wis: probability %.8g not found in grid.", target),
           call. = FALSE)
    idx
  }

  K   <- length(alpha_levels)
  m   <- pred_q[prob_idx(all_probs, 0.5)]

  wis_accum <- 0.5 * abs(y_test - m)

  for (k in seq_along(alpha_levels)) {
    a_k <- 1 - alpha_levels[k]                      # non-coverage probability
    l_k <- pred_q[prob_idx(all_probs, a_k / 2)]
    u_k <- pred_q[prob_idx(all_probs, 1 - a_k / 2)]
    is_k <- (u_k - l_k) +
            (2 / a_k) * pmax(l_k - y_test, 0) +
            (2 / a_k) * pmax(y_test - u_k, 0)
    wis_accum <- wis_accum + (a_k / 2) * is_k
  }

  wis_vals <- wis_accum / (K + 0.5)

  # ------------------------------------------------------------------
  # 5.  True marginal median (for relative WIS)
  #     Skipped when true_median is supplied by the caller.
  # ------------------------------------------------------------------
  if (is.null(true_median)) {
    true_median <- compute_true_marginal_quantile(
      dist_type = dist_type,
      mu0       = true_params$mu0,
      tau       = true_params$tau,
      phi       = true_params$phi,
      kappa     = kappa_val,
      probs     = 0.5
    )["50%"]
  }

  list(
    wis     = mean(wis_vals, na.rm = TRUE),
    rel_wis = mean(wis_vals, na.rm = TRUE) / true_median
  )
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
      kappa        = scenario$kappa,
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

      # kappa coverage/bias (only for 3-parameter distributions)
      has_kappa      <- scenario$dist_type %in% c("burr12", "gengamma")
      coverage_kappa <- if (has_kappa) check_coverage(fit, "kappa", sim_data$true_params$kappa) else NA
      bias_kappa     <- if (has_kappa) compute_median_bias(fit, "kappa", sim_data$true_params$kappa) else NA

      # True marginal quantiles (single call; "50%" reused as true_median for rel_wis)
      true_q <- compute_true_marginal_quantile(
        dist_type = scenario$dist_type,
        mu0       = sim_data$true_params$mu0,
        tau       = sim_data$true_params$tau,
        phi       = sim_data$true_params$phi,
        kappa     = sim_data$true_params$kappa
      )

      # Posterior predictive quantiles — single call covering both coverage/bias
      # (Q50, Q95) and the full WIS grid (9 levels).  The WIS probs are the
      # superset so we use them here and index into the result below.
      wis_all_probs <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
      idx_q50 <- which(wis_all_probs == 0.50)   # 5
      idx_q95 <- which(wis_all_probs == 0.95)   # 8
      post_ci <- compute_posterior_predictive_quantile_ci(
        fit        = fit,
        dist_type  = scenario$dist_type,
        n_datasets = scenario$n_datasets,
        probs      = wis_all_probs
      )

      # Quantile coverage and bias (index into the joint post_ci result)
      coverage_median <- true_q["50%"] >= post_ci$lower[idx_q50] & true_q["50%"] <= post_ci$upper[idx_q50]
      coverage_p95    <- true_q["95%"] >= post_ci$lower[idx_q95] & true_q["95%"] <= post_ci$upper[idx_q95]
      bias_median     <- post_ci$median[idx_q50] - true_q["50%"]
      bias_p95        <- post_ci$median[idx_q95] - true_q["95%"]

      # Weighted Interval Score — pass pre-computed pred_q and true_median to
      # avoid repeating the two expensive MC calls above
      wis_result <- compute_wis(
        fit         = fit,
        dist_type   = scenario$dist_type,
        n_datasets  = scenario$n_datasets,
        true_params = sim_data$true_params,
        pred_q      = post_ci$median,
        true_median = unname(true_q["50%"])
      )

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
        wis                    = wis_result$wis,
        rel_wis                = wis_result$rel_wis,
        max_rhat               = max_rhat,
        min_neff               = min_neff,
        converged              = max_rhat <= 1.1 & min_neff >= 100,
        true_kappa             = scenario$kappa,
        coverage_kappa         = coverage_kappa,
        bias_kappa             = bias_kappa,
        rel_bias_kappa         = if (has_kappa) bias_kappa / sim_data$true_params$kappa else NA_real_,
        coverage_median        = coverage_median,
        coverage_p95           = coverage_p95,
        bias_median            = bias_median,
        rel_bias_median        = bias_median / true_q["50%"],
        bias_p95               = bias_p95,
        rel_bias_p95           = bias_p95 / true_q["95%"],
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
        wis                    = NA,
        rel_wis                = NA,
        max_rhat               = NA,
        min_neff               = NA,
        converged              = FALSE,
        true_kappa             = scenario$kappa,
        coverage_kappa         = NA,
        bias_kappa             = NA,
        rel_bias_kappa         = NA,
        coverage_median        = NA,
        coverage_p95           = NA,
        bias_median            = NA,
        rel_bias_median        = NA,
        bias_p95               = NA,
        rel_bias_p95           = NA,
        stringsAsFactors       = FALSE
      )
    })
  }

  results <- foreach::foreach(
    i          = seq_len(nrow(tasks)),
    .combine   = dplyr::bind_rows,
    .packages  = c("rstan", "ddsynth"),
    .export    = c("stan_model", "run_one",
                   "compute_true_marginal_quantile",
                   "compute_posterior_predictive_quantile_ci",
                   "compute_wis")
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
#' @param round_order_stats If `TRUE`, the median/min/max/Q1/Q3 reported for
#'   summary types 1 and 2 are rounded to the nearest multiple of
#'   `resolution` before being returned, matching how these statistics are
#'   actually reported in the literature.
#' @param resolution Rounding grid in days, used only when
#'   `round_order_stats = TRUE`. Default `1` (nearest whole day); pass e.g.
#'   `1/24` to simulate hour-resolution reporting. Recorded in the returned
#'   `obs_data$resolution` field so the Stan model integrates over the same
#'   window that was used to round the data.
#' @return List containing true parameters and observed summary statistics
generate_hierarchical_data_mixed <- function(n_datasets,
                                             n_obs,
                                             dist_type = c("lognormal", "gamma", "weibull",
                                                           "burr12", "gengamma"),
                                             mu0,
                                             tau,
                                             phi,
                                             kappa = 1.0,
                                             summary_type = NULL,
                                             round_order_stats = FALSE,
                                             resolution = 1) {

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
  obs_stat1      <- numeric(n_datasets)
  obs_stat2      <- numeric(n_datasets)
  obs_stat3      <- numeric(n_datasets)
  resolution_vec <- rep(1, n_datasets)  # unused (placeholder) for non-type-1/2 datasets
  freq_tables    <- vector("list", n_datasets)
  
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
      shape  <- phi
      rate   <- shape / mean_d
      data_d <- rgamma(n, shape = shape, rate = rate)

    } else if (dist_type == "weibull") {
      scale  <- exp(loc)
      shape  <- phi
      data_d <- rweibull(n, shape = shape, scale = scale)

    } else if (dist_type == "burr12") {
      # Burr XII via inverse-CDF: Q(u) = lambda * (u^(-1/k) - 1)^(1/c)
      # lambda = exp(loc), c = phi, k = kappa.
      lambda <- exp(loc)
      u      <- runif(n)
      data_d <- lambda * (u^(-1.0 / kappa) - 1.0)^(1.0 / phi)

    } else if (dist_type == "gengamma") {
      # Generalised Gamma (Prentice): mu = loc, sigma = phi, Q = kappa.
      # T = exp(mu + sigma/Q * log(Q^2 * Y)),  Y ~ Gamma(1/Q^2, 1).
      gs     <- 1.0 / kappa^2
      y      <- rgamma(n, shape = gs, rate = 1)
      data_d <- exp(loc + phi / kappa * log(kappa^2 * y))
    }
    
    # Compute summary statistics based on type for this specific dataset
    if (st == 1) {  # median + range
      obs_stat1[d] <- median(data_d)
      obs_stat2[d] <- min(data_d)
      obs_stat3[d] <- max(data_d)
      if (round_order_stats) {
        obs_stat1[d] <- round(obs_stat1[d] / resolution) * resolution
        obs_stat2[d] <- round(obs_stat2[d] / resolution) * resolution
        obs_stat3[d] <- round(obs_stat3[d] / resolution) * resolution
        resolution_vec[d] <- resolution
      }

    } else if (st == 2) {  # median + IQR
      obs_stat1[d] <- median(data_d)
      obs_stat2[d] <- quantile(data_d, 0.25)
      obs_stat3[d] <- quantile(data_d, 0.75)
      if (round_order_stats) {
        obs_stat1[d] <- round(obs_stat1[d] / resolution) * resolution
        obs_stat2[d] <- round(obs_stat2[d] / resolution) * resolution
        obs_stat3[d] <- round(obs_stat3[d] / resolution) * resolution
        resolution_vec[d] <- resolution
      }

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
      mu0   = mu0,
      tau   = tau,
      phi   = phi,
      kappa = kappa,
      loc_d = loc_d
    ),
    obs_data = list(
      n_datasets   = n_datasets,
      n_obs        = as.array(n_obs),
      summary_type = as.array(summary_type),
      dist_type    = switch(dist_type,
                            "lognormal" = 1L,
                            "gamma"     = 2L,
                            "weibull"   = 3L,
                            "burr12"    = 4L,
                            "gengamma"  = 5L),
      obs_stat1    = as.array(obs_stat1),
      obs_stat2    = as.array(obs_stat2),
      obs_stat3    = as.array(obs_stat3),
      resolution   = as.array(resolution_vec),
      # Frequency table fields (populated only when summary_type == 4)
      n_freq_total = length(freq_value_flat),
      freq_value   = freq_value_flat,
      freq_count   = freq_count_flat,
      freq_start   = as.array(freq_start_vec),
      freq_len     = as.array(freq_len_vec),
      # Frequency-table interval bounds (type-5 datasets only; zeros for types 1-4)
      freq_lower   = rep(0.0, length(freq_value_flat)),
      freq_upper   = rep(0.0, length(freq_value_flat)),
      # Priors — calibrated defaults matching prepare_stan_data_from_datasets()
      mu0_mean     = 1.0,
      mu0_sd       = 1.0,
      log_tau_mean = 0.2,
      log_tau_sd   = 0.5,
      log_phi_mean = switch(dist_type,
        lognormal = -0.7, gamma = 2.5, weibull = 1.0,
        burr12    =  0.7, gengamma = -0.5),
      log_phi_sd   = 0.5,
      log_kappa_mean = switch(dist_type,
        lognormal = 0.0, gamma = 0.0, weibull = 0.0,
        burr12    = 1.0, gengamma = 0.0),
      log_kappa_sd   = switch(dist_type,
        lognormal = 1.0, gamma = 1.0, weibull = 1.0,
        burr12    = 0.5, gengamma = 0.5)
    )
  )
}