# Tests for extract_quantiles() in R/utils.R

# ---------------------------------------------------------------------------
# Local fixture: toy CDF summary + matrix from LN(meanlog=2, sdlog=0.5)
# ---------------------------------------------------------------------------
.make_cdf_objects <- function(n_draws = 50) {
  x_seq <- seq(0.01, 20, length.out = 200)
  set.seed(42)
  cdf_mat <- t(replicate(n_draws, {
    mu_d <- rnorm(1, 2, 0.1)
    s_d  <- exp(rnorm(1, log(0.5), 0.05))
    plnorm(x_seq, meanlog = mu_d, sdlog = s_d)
  }))
  cdf_summary <- data.frame(
    x      = x_seq,
    median = apply(cdf_mat, 2, median),
    mean   = apply(cdf_mat, 2, mean),
    low    = apply(cdf_mat, 2, quantile, 0.025),
    high   = apply(cdf_mat, 2, quantile, 0.975)
  )
  list(summary = cdf_summary, cdf_mat = cdf_mat)
}

# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

test_that("returns a data frame with columns quantile, quantile_label, x_low, x_high", {
  objs   <- .make_cdf_objects()
  result <- suppressWarnings(
    extract_quantiles(objs$summary, probs = c(0.5, 0.95), cdf_mat = objs$cdf_mat)
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("quantile", "quantile_label", "x_low", "x_high") %in% names(result)))
})

test_that("returns one row per requested probability (cdf_mat path)", {
  objs   <- .make_cdf_objects()
  result <- suppressWarnings(
    extract_quantiles(objs$summary, probs = c(0.25, 0.5, 0.75, 0.95), cdf_mat = objs$cdf_mat)
  )
  expect_equal(nrow(result), 4L)
  expect_equal(sort(result$quantile), c(0.25, 0.50, 0.75, 0.95))
})

test_that("quantile_label is 'Q' followed by the probability as a percentage", {
  objs   <- .make_cdf_objects()
  result <- suppressWarnings(
    extract_quantiles(objs$summary, probs = c(0.5, 0.95), cdf_mat = objs$cdf_mat)
  )
  expect_equal(sort(result$quantile_label), c("Q50", "Q95"))
})

# ---------------------------------------------------------------------------
# Bounds ordering
# ---------------------------------------------------------------------------

test_that("x_low <= x_high for every probability", {
  objs   <- .make_cdf_objects()
  result <- suppressWarnings(
    extract_quantiles(objs$summary, probs = c(0.5, 0.95), cdf_mat = objs$cdf_mat)
  )
  expect_true(all(result$x_low <= result$x_high))
})

# ---------------------------------------------------------------------------
# Accuracy
# ---------------------------------------------------------------------------

test_that("quantile midpoints are plausible for LN(meanlog=2, sdlog=0.5)", {
  # True median = exp(2) ≈ 7.4; true 95th pct ≈ 18.5
  objs   <- .make_cdf_objects()
  result <- suppressWarnings(
    extract_quantiles(objs$summary, probs = c(0.5, 0.95), cdf_mat = objs$cdf_mat)
  )
  q50_mid <- mean(c(result[result$quantile == 0.50, "x_low"],
                    result[result$quantile == 0.50, "x_high"]))
  q95_mid <- mean(c(result[result$quantile == 0.95, "x_low"],
                    result[result$quantile == 0.95, "x_high"]))
  expect_gt(q50_mid, 5);  expect_lt(q50_mid, 11)
  expect_gt(q95_mid, 12); expect_lt(q95_mid, 25)
})

# ---------------------------------------------------------------------------
# Fallback path (cdf_mat = NULL)
# ---------------------------------------------------------------------------

test_that("fallback path (cdf_mat = NULL) returns same column structure", {
  objs   <- .make_cdf_objects()
  result <- extract_quantiles(objs$summary, probs = c(0.5, 0.95), cdf_mat = NULL)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("quantile", "quantile_label", "x_low", "x_high") %in% names(result)))
  expect_equal(nrow(result), 2L)
})

# ---------------------------------------------------------------------------
# Warning when x_seq is too narrow
# ---------------------------------------------------------------------------

test_that("warns when more than 5% of draws do not reach the requested probability", {
  x_seq      <- seq(0, 5, length.out = 100)   # too narrow for p = 0.9999
  cdf_mat    <- matrix(rep(plnorm(x_seq, 2, 0.5), 30), nrow = 30, byrow = TRUE)
  summary_df <- data.frame(x = x_seq,
                           median = plnorm(x_seq, 2, 0.5),
                           mean   = plnorm(x_seq, 2, 0.5),
                           low    = plnorm(x_seq, 2, 0.5),
                           high   = plnorm(x_seq, 2, 0.5))
  expect_warning(
    extract_quantiles(summary_df, probs = 0.9999, cdf_mat = cdf_mat),
    regexp = "did not reach"
  )
})
