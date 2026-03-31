# Tests for summarise_parameters() in R/bsl_data_synthesis.R

# Local fixture: a posterior sample list centred at known values
.post_list <- function(n = 400, mu0 = 2, log_tau = log(0.5), log_phi = log(3)) {
  set.seed(7)
  post <- cbind(rnorm(n, mu0,     0.1),
                rnorm(n, log_tau, 0.1),
                rnorm(n, log_phi, 0.1))
  list(lognormal = post)
}

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("returns a data frame with the expected columns", {
  result <- ddsynth:::summarise_parameters(.post_list())
  expect_s3_class(result, "data.frame")
  expect_true(all(c("model", "parameter", "mean", "median",
                    "ci_95", "low_2.5", "high_97.5") %in% names(result)))
})

test_that("returns exactly 3 rows (mu0, tau, phi) for a single model", {
  result <- ddsynth:::summarise_parameters(.post_list())
  expect_equal(nrow(result), 3L)
  expect_equal(as.character(result$parameter), c("mu0", "tau", "phi"))
})

test_that("returns 3 rows per model when multiple models are supplied", {
  post      <- cbind(rnorm(200, 2, 0.1), rnorm(200, 0, 0.1), rnorm(200, 1, 0.1))
  post_list <- list(lognormal = post, gamma = post)
  result    <- ddsynth:::summarise_parameters(post_list)
  expect_equal(nrow(result), 6L)
  expect_equal(unique(result$model), c("lognormal", "gamma"))
})

# ---------------------------------------------------------------------------
# Scale: tau and phi must be on the natural (exp) scale
# ---------------------------------------------------------------------------

test_that("tau and phi summaries are on the natural (exp) scale", {
  # log_tau centred at 0 → tau ≈ 1; log_phi centred at log(3) → phi ≈ 3
  result  <- ddsynth:::summarise_parameters(.post_list(n = 2000, log_tau = 0, log_phi = log(3)))
  tau_row <- result[result$parameter == "tau", ]
  phi_row <- result[result$parameter == "phi", ]
  expect_equal(tau_row$median, 1.0, tolerance = 0.05)
  expect_equal(phi_row$median, 3.0, tolerance = 0.10)
})

# ---------------------------------------------------------------------------
# Credible interval properties
# ---------------------------------------------------------------------------

test_that("95% CrI brackets the true mu0 value", {
  result  <- ddsynth:::summarise_parameters(.post_list(n = 2000, mu0 = 2))
  mu0_row <- result[result$parameter == "mu0", ]
  expect_lte(mu0_row$low_2.5,   2.0)
  expect_gte(mu0_row$high_97.5, 2.0)
})

test_that("low_2.5 <= median <= high_97.5 for every parameter", {
  result <- ddsynth:::summarise_parameters(.post_list())
  expect_true(all(result$low_2.5 <= result$median))
  expect_true(all(result$median  <= result$high_97.5))
})

# ---------------------------------------------------------------------------
# ci_95 formatting
# ---------------------------------------------------------------------------

test_that("ci_95 is a character column containing an em dash", {
  result <- ddsynth:::summarise_parameters(.post_list())
  expect_type(result$ci_95, "character")
  expect_true(all(grepl("\u2014", result$ci_95)))
})

# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

test_that("errors when the posterior matrix has fewer than 3 columns", {
  bad_post <- list(bad = matrix(rnorm(200), ncol = 2))
  expect_error(ddsynth:::summarise_parameters(bad_post), regexp = "3 columns")
})
