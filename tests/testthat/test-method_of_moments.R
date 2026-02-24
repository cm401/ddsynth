test_that("fit_gamma_mom returns correct parameters", {
  # For Gamma: shape = mean^2/var, rate = mean/var
  fit <- fit_gamma_mom(mean = 4, variance = 2)
  expect_equal(fit$shape, 8, tolerance = 1e-8)
  expect_equal(fit$rate, 2, tolerance = 1e-8)
  expect_equal(fit$scale, 0.5, tolerance = 1e-8)
  expect_equal(fit$mean, 4, tolerance = 1e-8)
  expect_equal(fit$variance, 2, tolerance = 1e-8)
})

test_that("fit_gamma_mom accepts sd argument", {
  fit_v <- fit_gamma_mom(mean = 4, variance = 4)
  fit_s <- fit_gamma_mom(mean = 4, sd = 2)
  expect_equal(fit_v$shape, fit_s$shape)
  expect_equal(fit_v$rate, fit_s$rate)
})

test_that("fit_gamma_mom errors on invalid inputs", {
  expect_error(fit_gamma_mom(mean = -1, variance = 1), "'mean' must be strictly positive")
  expect_error(fit_gamma_mom(mean = 1, variance = -1), "'variance' must be strictly positive")
  expect_error(fit_gamma_mom(mean = 1, variance = 1, sd = 1), "not both")
  expect_error(fit_gamma_mom(mean = 1), "must be supplied")
})

test_that("fit_lognormal_mom returns correct parameters", {
  mu   <- 5
  sig2 <- 4
  fit  <- fit_lognormal_mom(mean = mu, variance = sig2)
  # Recover implied mean and variance
  expect_equal(fit$mean, mu, tolerance = 1e-6)
  expect_equal(fit$variance, sig2, tolerance = 1e-6)
  expect_true(fit$sdlog > 0)
})

test_that("fit_lognormal_mom accepts sd argument", {
  fit_v <- fit_lognormal_mom(mean = 5, variance = 4)
  fit_s <- fit_lognormal_mom(mean = 5, sd = 2)
  expect_equal(fit_v$meanlog, fit_s$meanlog)
  expect_equal(fit_v$sdlog, fit_s$sdlog)
})

test_that("fit_weibull_mom returns correct parameters", {
  mu   <- 5
  sig2 <- 4
  fit  <- fit_weibull_mom(mean = mu, variance = sig2)
  expect_equal(fit$mean, mu, tolerance = 1e-5)
  expect_equal(fit$variance, sig2, tolerance = 1e-5)
  expect_true(fit$shape > 0)
  expect_true(fit$scale > 0)
})

test_that("fit_weibull_mom accepts sd argument", {
  fit_v <- fit_weibull_mom(mean = 5, variance = 4)
  fit_s <- fit_weibull_mom(mean = 5, sd = 2)
  expect_equal(fit_v$shape, fit_s$shape, tolerance = 1e-8)
  expect_equal(fit_v$scale, fit_s$scale, tolerance = 1e-8)
})
