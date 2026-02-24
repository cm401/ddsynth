test_that("mle_from_summary returns a list with expected elements for gamma", {
  fit <- mle_from_summary(mean = c(5, 5.2), sd = c(2, 1.8),
                          n = c(100, 80), distribution = "gamma")
  expect_type(fit, "list")
  expect_named(fit, c("distribution", "parameters", "loglik", "convergence"),
               ignore.order = TRUE)
  expect_equal(fit$distribution, "gamma")
  expect_named(fit$parameters, c("shape", "rate"))
  expect_true(all(fit$parameters > 0))
})

test_that("mle_from_summary returns a list with expected elements for lognormal", {
  fit <- mle_from_summary(mean = c(5, 5.2), sd = c(2, 1.8),
                          n = c(100, 80), distribution = "lognormal")
  expect_equal(fit$distribution, "lognormal")
  expect_named(fit$parameters, c("meanlog", "sdlog"))
  expect_true(fit$parameters[["sdlog"]] > 0)
})

test_that("mle_from_summary returns a list with expected elements for weibull", {
  fit <- mle_from_summary(mean = c(5, 5.2), sd = c(2, 1.8),
                          n = c(100, 80), distribution = "weibull")
  expect_equal(fit$distribution, "weibull")
  expect_named(fit$parameters, c("shape", "scale"))
  expect_true(all(fit$parameters > 0))
})

test_that("mle_from_summary errors if mean and variance lengths differ", {
  result <- tryCatch(
    mle_from_summary(mean = c(5, 5.2), sd = c(2), n = c(100, 80)),
    error = function(e) "error"
  )
  expect_equal(result, "error")
})

test_that("mle_from_summary uses equal weights when n is not supplied", {
  fit_no_n  <- mle_from_summary(mean = 5, sd = 2, distribution = "gamma")
  fit_n_one <- mle_from_summary(mean = 5, sd = 2, n = 1L, distribution = "gamma")
  expect_equal(fit_no_n$parameters, fit_n_one$parameters, tolerance = 1e-6)
})
