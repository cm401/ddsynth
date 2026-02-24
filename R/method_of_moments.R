#' Fit a Gamma distribution using the method of moments
#'
#' Estimates the shape and rate parameters of a Gamma distribution from the
#' observed mean and variance (or standard deviation) of a delay.
#'
#' @param mean Numeric scalar. The observed mean of the delay.
#' @param variance Numeric scalar. The observed variance of the delay.
#'   Provide either `variance` or `sd`, not both.
#' @param sd Numeric scalar. The observed standard deviation of the delay.
#'   Provide either `sd` or `variance`, not both.
#'
#' @return A list with components:
#'   \describe{
#'     \item{shape}{Estimated shape parameter of the Gamma distribution.}
#'     \item{rate}{Estimated rate parameter of the Gamma distribution.}
#'     \item{scale}{Estimated scale parameter (reciprocal of rate).}
#'     \item{mean}{Mean implied by the fitted parameters.}
#'     \item{variance}{Variance implied by the fitted parameters.}
#'   }
#'
#' @examples
#' fit <- fit_gamma_mom(mean = 5, variance = 4)
#' fit$shape
#' fit$rate
#'
#' @export
fit_gamma_mom <- function(mean, variance = NULL, sd = NULL) {
  variance <- .resolve_variance(variance, sd)
  .check_positive(mean, "mean")
  .check_positive(variance, "variance")

  shape <- mean^2 / variance
  rate  <- mean / variance

  list(
    shape    = shape,
    rate     = rate,
    scale    = 1 / rate,
    mean     = shape / rate,
    variance = shape / rate^2
  )
}


#' Fit a Log-Normal distribution using the method of moments
#'
#' Estimates the meanlog and sdlog parameters of a Log-Normal distribution from
#' the observed mean and variance (or standard deviation) of a delay.
#'
#' @inheritParams fit_gamma_mom
#'
#' @return A list with components:
#'   \describe{
#'     \item{meanlog}{Estimated mean of the log of the delay.}
#'     \item{sdlog}{Estimated standard deviation of the log of the delay.}
#'     \item{mean}{Mean implied by the fitted parameters.}
#'     \item{variance}{Variance implied by the fitted parameters.}
#'   }
#'
#' @examples
#' fit <- fit_lognormal_mom(mean = 5, sd = 2)
#' fit$meanlog
#' fit$sdlog
#'
#' @export
fit_lognormal_mom <- function(mean, variance = NULL, sd = NULL) {
  variance <- .resolve_variance(variance, sd)
  .check_positive(mean, "mean")
  .check_positive(variance, "variance")

  sdlog   <- sqrt(log(1 + variance / mean^2))
  meanlog <- log(mean) - sdlog^2 / 2

  list(
    meanlog  = meanlog,
    sdlog    = sdlog,
    mean     = exp(meanlog + sdlog^2 / 2),
    variance = (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
  )
}


#' Fit a Weibull distribution using the method of moments
#'
#' Estimates the shape and scale parameters of a Weibull distribution from the
#' observed mean and variance (or standard deviation) of a delay.
#'
#' The shape parameter is solved numerically, so convergence is not guaranteed
#' for all inputs.
#'
#' @inheritParams fit_gamma_mom
#'
#' @return A list with components:
#'   \describe{
#'     \item{shape}{Estimated shape parameter of the Weibull distribution.}
#'     \item{scale}{Estimated scale parameter of the Weibull distribution.}
#'     \item{mean}{Mean implied by the fitted parameters.}
#'     \item{variance}{Variance implied by the fitted parameters.}
#'   }
#'
#' @examples
#' fit <- fit_weibull_mom(mean = 5, variance = 4)
#' fit$shape
#' fit$scale
#'
#' @export
fit_weibull_mom <- function(mean, variance = NULL, sd = NULL) {
  variance <- .resolve_variance(variance, sd)
  .check_positive(mean, "mean")
  .check_positive(variance, "variance")

  # cv^2 = Gamma(1 + 2/k) / Gamma(1 + 1/k)^2 - 1
  cv_sq <- variance / mean^2
  obj <- function(k) {
    gamma(1 + 2 / k) / gamma(1 + 1 / k)^2 - 1 - cv_sq
  }
  fit  <- uniroot(obj, interval = c(1e-4, 1e4), tol = 1e-8)
  k    <- fit$root
  lam  <- mean / gamma(1 + 1 / k)

  list(
    shape    = k,
    scale    = lam,
    mean     = lam * gamma(1 + 1 / k),
    variance = lam^2 * (gamma(1 + 2 / k) - gamma(1 + 1 / k)^2)
  )
}
