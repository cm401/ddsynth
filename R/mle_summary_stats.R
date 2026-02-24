#' Fit a delay distribution via maximum likelihood from summary statistics
#'
#' Fits a parametric delay distribution by maximising a log-likelihood that is
#' constructed from summary statistics (mean and variance/SD) from one or more
#' previous studies.  A large-sample normal approximation is used: the sample
#' mean is treated as approximately N(μ, σ²/n) and the sample variance as
#' approximately N(σ², 2σ⁴/(n-1)), where μ and σ² are the distributional mean
#' and variance implied by the candidate parameters.
#'
#' @param mean Numeric scalar or vector. Observed mean(s) of the delay.
#' @param variance Numeric scalar or vector. Observed variance(s) of the delay.
#'   Provide either `variance` or `sd`, not both.
#' @param sd Numeric scalar or vector. Observed standard deviation(s) of the
#'   delay.  Provide either `sd` or `variance`, not both.
#' @param n Integer scalar or vector. Sample size(s) associated with each set of
#'   summary statistics.  Used to weight the contribution of each study to the
#'   overall log-likelihood.  Defaults to equal weighting.
#' @param distribution Character string specifying the parametric family.  One
#'   of `"gamma"`, `"lognormal"`, or `"weibull"`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{distribution}{The fitted distribution family.}
#'     \item{parameters}{Named numeric vector of the fitted parameters.}
#'     \item{loglik}{The maximised log-likelihood value.}
#'     \item{convergence}{Convergence code from the optimiser (0 = success).}
#'   }
#'
#' @examples
#' fit <- mle_from_summary(mean = c(5, 5.2), sd = c(2, 1.8), n = c(100, 80),
#'                         distribution = "gamma")
#' fit$parameters
#'
#' @export
mle_from_summary <- function(mean, variance = NULL, sd = NULL,
                             n = NULL, distribution = c("gamma", "lognormal", "weibull")) {
  distribution <- match.arg(distribution)
  variance <- .resolve_variance(variance, sd)
  stopifnot(length(mean) == length(variance))

  if (is.null(n)) {
    n <- rep(1L, length(mean))
  }
  stopifnot(length(n) == length(mean))

  # Use method-of-moments estimates as starting values
  start <- switch(distribution,
    gamma     = {
      mom <- fit_gamma_mom(mean[1], variance = variance[1])
      c(log_shape = log(mom$shape), log_rate = log(mom$rate))
    },
    lognormal = {
      mom <- fit_lognormal_mom(mean[1], variance = variance[1])
      c(meanlog = mom$meanlog, log_sdlog = log(mom$sdlog))
    },
    weibull   = {
      mom <- fit_weibull_mom(mean[1], variance = variance[1])
      c(log_shape = log(mom$shape), log_scale = log(mom$scale))
    }
  )

  # Large-sample normal approximation:
  #   sample mean   X_bar ~ N(mu,   sigma^2 / n)
  #   sample variance S^2 ~ N(sigma^2, 2 * sigma^4 / (n - 1))
  # This approximation is valid for any distribution with finite 4th moment.
  neg_loglik <- switch(distribution,
    gamma = function(par) {
      shape <- exp(par[["log_shape"]])
      rate  <- exp(par[["log_rate"]])
      mu    <- shape / rate
      sig2  <- shape / rate^2
      n_m1  <- pmax(n - 1L, 1L)
      ll_mean <- stats::dnorm(mean,     mean = mu,   sd = sqrt(sig2 / n),        log = TRUE)
      ll_var  <- stats::dnorm(variance, mean = sig2, sd = sqrt(2 * sig2^2 / n_m1), log = TRUE)
      -sum(ll_mean + ll_var)
    },
    lognormal = function(par) {
      ml   <- par[["meanlog"]]
      sl   <- exp(par[["log_sdlog"]])
      mu   <- exp(ml + sl^2 / 2)
      sig2 <- (exp(sl^2) - 1) * exp(2 * ml + sl^2)
      n_m1 <- pmax(n - 1L, 1L)
      ll_mean <- stats::dnorm(mean,     mean = mu,   sd = sqrt(sig2 / n),        log = TRUE)
      ll_var  <- stats::dnorm(variance, mean = sig2, sd = sqrt(2 * sig2^2 / n_m1), log = TRUE)
      -sum(ll_mean + ll_var)
    },
    weibull = function(par) {
      k    <- exp(par[["log_shape"]])
      lam  <- exp(par[["log_scale"]])
      mu   <- lam * gamma(1 + 1 / k)
      sig2 <- lam^2 * (gamma(1 + 2 / k) - gamma(1 + 1 / k)^2)
      n_m1 <- pmax(n - 1L, 1L)
      ll_mean <- stats::dnorm(mean,     mean = mu,   sd = sqrt(sig2 / n),        log = TRUE)
      ll_var  <- stats::dnorm(variance, mean = sig2, sd = sqrt(2 * sig2^2 / n_m1), log = TRUE)
      -sum(ll_mean + ll_var)
    }
  )

  opt <- stats::optim(start, neg_loglik, method = "BFGS",
                      control = list(maxit = 1000))

  params <- switch(distribution,
    gamma     = c(shape = exp(opt$par[["log_shape"]]),
                  rate  = exp(opt$par[["log_rate"]])),
    lognormal = c(meanlog = opt$par[["meanlog"]],
                  sdlog   = exp(opt$par[["log_sdlog"]])),
    weibull   = c(shape = exp(opt$par[["log_shape"]]),
                  scale = exp(opt$par[["log_scale"]]))
  )

  list(
    distribution = distribution,
    parameters   = params,
    loglik       = -opt$value,
    convergence  = opt$convergence
  )
}
