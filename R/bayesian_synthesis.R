#' Bayesian synthesis of delay distributions from summary statistics
#'
#' Combines summary statistics from multiple studies using a Bayesian
#' hierarchical approach to synthesise a posterior predictive delay
#' distribution. Currently implements a conjugate Normal-Gamma model for the
#' log-transformed delay, which corresponds to a Log-Normal marginal likelihood.
#'
#' @param mean Numeric vector. Observed mean delay in each study.
#' @param variance Numeric scalar or vector. Observed variance of the delay.
#'   Provide either `variance` or `sd`, not both.
#' @param sd Numeric scalar or vector. Observed standard deviation of the delay.
#'   Provide either `sd` or `variance`, not both.
#' @param n Integer vector. Sample size of each study.
#' @param prior_mean_meanlog Numeric scalar. Prior mean for the log-scale mean
#'   parameter (meanlog). Defaults to `0`.
#' @param prior_sd_meanlog Numeric scalar. Prior standard deviation for
#'   meanlog. Defaults to `10` (weakly informative).
#' @param prior_shape_precision Numeric scalar. Shape of the Gamma prior on the
#'   log-scale precision (1/sdlog^2). Defaults to `1`.
#' @param prior_rate_precision Numeric scalar. Rate of the Gamma prior on the
#'   log-scale precision. Defaults to `1`.
#' @param n_samples Integer scalar. Number of posterior samples to draw.
#'   Defaults to `10000`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{distribution}{`"lognormal"` — the marginal family of the posterior
#'       predictive.}
#'     \item{posterior_samples}{A data frame of posterior samples with columns
#'       `meanlog` and `sdlog`.}
#'     \item{posterior_mean_meanlog}{Posterior mean of meanlog.}
#'     \item{posterior_sd_meanlog}{Posterior standard deviation of meanlog.}
#'     \item{posterior_mean_sdlog}{Posterior mean of sdlog.}
#'   }
#'
#' @examples
#' result <- bayesian_synthesis(
#'   mean = c(5, 5.2, 4.8),
#'   sd   = c(2, 1.8, 2.1),
#'   n    = c(100, 80, 120)
#' )
#' head(result$posterior_samples)
#'
#' @export
bayesian_synthesis <- function(mean, variance = NULL, sd = NULL, n,
                               prior_mean_meanlog     = 0,
                               prior_sd_meanlog       = 10,
                               prior_shape_precision  = 1,
                               prior_rate_precision   = 1,
                               n_samples              = 10000L) {
  variance <- .resolve_variance(variance, sd)
  stopifnot(length(mean) == length(variance), length(mean) == length(n))
  .check_positive(prior_sd_meanlog, "prior_sd_meanlog")
  .check_positive(prior_shape_precision, "prior_shape_precision")
  .check_positive(prior_rate_precision, "prior_rate_precision")

  # Convert observed moments to log-scale summaries
  sdlog_obs   <- sqrt(log(1 + variance / mean^2))
  meanlog_obs <- log(mean) - sdlog_obs^2 / 2

  # Conjugate Normal-Gamma update for (meanlog, precision)
  # Prior: meanlog | precision ~ N(mu0, 1/(kappa0 * precision))
  #         precision ~ Gamma(alpha0, beta0)
  mu0    <- prior_mean_meanlog
  kappa0 <- 1 / prior_sd_meanlog^2
  alpha0 <- prior_shape_precision
  beta0  <- prior_rate_precision

  k      <- length(mean)
  x_bar  <- sum(n * meanlog_obs) / sum(n)
  n_tot  <- sum(n)

  kappa_n <- kappa0 + n_tot
  mu_n    <- (kappa0 * mu0 + n_tot * x_bar) / kappa_n
  alpha_n <- alpha0 + n_tot / 2
  beta_n  <- beta0 +
    0.5 * sum(n * sdlog_obs^2) +
    0.5 * (kappa0 * n_tot / kappa_n) * (x_bar - mu0)^2

  # Draw posterior samples
  precision_samples <- stats::rgamma(n_samples, shape = alpha_n, rate = beta_n)
  sdlog_samples     <- 1 / sqrt(precision_samples)
  meanlog_samples   <- stats::rnorm(n_samples,
                                    mean = mu_n,
                                    sd   = sdlog_samples / sqrt(kappa_n))

  list(
    distribution          = "lognormal",
    posterior_samples     = data.frame(meanlog = meanlog_samples,
                                       sdlog   = sdlog_samples),
    posterior_mean_meanlog = mean(meanlog_samples),
    posterior_sd_meanlog   = stats::var(meanlog_samples)^0.5,
    posterior_mean_sdlog   = mean(sdlog_samples)
  )
}
