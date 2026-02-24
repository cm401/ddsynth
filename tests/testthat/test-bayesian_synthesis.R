test_that("bayesian_synthesis returns expected structure", {
  set.seed(42)
  result <- bayesian_synthesis(
    mean = c(5, 5.2, 4.8),
    sd   = c(2, 1.8, 2.1),
    n    = c(100L, 80L, 120L)
  )
  expect_type(result, "list")
  expect_equal(result$distribution, "lognormal")
  expect_s3_class(result$posterior_samples, "data.frame")
  expect_named(result$posterior_samples, c("meanlog", "sdlog"))
  expect_equal(nrow(result$posterior_samples), 10000L)
  expect_true(all(result$posterior_samples$sdlog > 0))
})

test_that("bayesian_synthesis posterior mean is near MoM estimate", {
  set.seed(123)
  mu   <- 5
  sig2 <- 4
  mom  <- fit_lognormal_mom(mean = mu, variance = sig2)
  result <- bayesian_synthesis(
    mean = rep(mu, 5),
    sd   = rep(sqrt(sig2), 5),
    n    = rep(200L, 5)
  )
  # With lots of data the posterior should be close to the MoM estimate
  expect_equal(result$posterior_mean_meanlog, mom$meanlog, tolerance = 0.05)
})

test_that("bayesian_synthesis errors if lengths are inconsistent", {
  result <- tryCatch(
    bayesian_synthesis(mean = c(5, 5.2), sd = c(2, 1.8), n = c(100L)),
    error = function(e) "error"
  )
  expect_equal(result, "error")
})

test_that("bayesian_synthesis errors if both variance and sd are given", {
  expect_error(
    bayesian_synthesis(mean = 5, variance = 4, sd = 2, n = 100L),
    "not both"
  )
})

test_that("bayesian_synthesis n_samples controls output size", {
  set.seed(1)
  result <- bayesian_synthesis(mean = 5, sd = 2, n = 50L, n_samples = 500L)
  expect_equal(nrow(result$posterior_samples), 500L)
})
