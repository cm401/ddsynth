# Internal helper utilities for the delaydistributions package.
# These functions are NOT exported.

#' Resolve variance from variance/sd arguments
#'
#' @param variance Numeric or NULL.
#' @param sd Numeric or NULL.
#' @return Numeric variance.
#' @noRd
.resolve_variance <- function(variance, sd) {
  if (!is.null(variance) && !is.null(sd)) {
    stop("Provide either 'variance' or 'sd', not both.", call. = FALSE)
  }
  if (is.null(variance) && is.null(sd)) {
    stop("One of 'variance' or 'sd' must be supplied.", call. = FALSE)
  }
  if (!is.null(sd)) {
    return(sd^2)
  }
  variance
}

#' Check that a value is strictly positive
#'
#' @param x Numeric scalar.
#' @param name Character scalar. Name of the argument (for error messages).
#' @noRd
.check_positive <- function(x, name) {
  if (any(x <= 0)) {
    stop(sprintf("'%s' must be strictly positive.", name), call. = FALSE)
  }
  invisible(x)
}
