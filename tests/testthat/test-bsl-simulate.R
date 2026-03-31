# Tests for bsl_simulate_dataset() and bsl_bridge_estimate() in R/bsl_data_synthesis.R

# ---------------------------------------------------------------------------
# bsl_simulate_dataset()
# ---------------------------------------------------------------------------

test_that("returns a numeric vector of the requested length for all distributions", {
  for (dist in c("lognormal", "gamma", "weibull")) {
    x <- ddsynth:::bsl_simulate_dataset(dist, loc = log(5), phi = 1.5, n = 100)
    expect_type(x, "double")
    expect_length(x, 100L)
  }
})

test_that("all returned values are strictly positive", {
  for (dist in c("lognormal", "gamma", "weibull")) {
    x <- ddsynth:::bsl_simulate_dataset(dist, loc = log(5), phi = 1.5, n = 200)
    expect_true(all(x > 0))
  }
})

test_that("lognormal: log-scale mean and SD match the supplied parameters", {
  set.seed(101)
  mu <- log(8); phi <- 0.4
  x  <- ddsynth:::bsl_simulate_dataset("lognormal", loc = mu, phi = phi, n = 5000)
  expect_equal(mean(log(x)), mu,  tolerance = 0.05)
  expect_equal(sd(log(x)),   phi, tolerance = 0.05)
})

test_that("gamma: natural-scale mean matches exp(loc) and shape ≈ phi", {
  set.seed(202)
  loc <- log(7); phi <- 3.0
  x   <- ddsynth:::bsl_simulate_dataset("gamma", loc = loc, phi = phi, n = 5000)
  expect_equal(mean(x), exp(loc), tolerance = 0.10)
  # MoM shape estimate: (mean/sd)^2
  expect_equal((mean(x) / sd(x))^2, phi, tolerance = 0.30)
})

test_that("weibull: median matches scale * (log 2)^(1/shape)", {
  set.seed(303)
  loc <- log(6); phi <- 2.0
  x   <- ddsynth:::bsl_simulate_dataset("weibull", loc = loc, phi = phi, n = 5000)
  expected_median <- exp(loc) * (log(2))^(1 / phi)
  expect_equal(median(x), expected_median, tolerance = 0.10)
})

test_that("output contains no non-finite values under normal parameters", {
  x <- ddsynth:::bsl_simulate_dataset("lognormal", loc = 0, phi = 1, n = 500)
  expect_true(all(is.finite(x)))
})

test_that("errors with an informative message for an unknown distribution name", {
  expect_error(
    ddsynth:::bsl_simulate_dataset("pareto", loc = 0, phi = 1, n = 100),
    regexp = "Unknown dist"
  )
})

# ---------------------------------------------------------------------------
# bsl_bridge_estimate()
# ---------------------------------------------------------------------------

test_that("returns a single finite numeric scalar", {
  loglik <- rnorm(200, mean = -100, sd = 5)
  result <- bsl_bridge_estimate(loglik)
  expect_type(result, "double")
  expect_length(result, 1L)
  expect_true(is.finite(result))
})

test_that("log-mean-exp of a constant vector equals that constant", {
  expect_equal(bsl_bridge_estimate(rep(-50.0, 100)), -50.0, tolerance = 1e-10)
})

test_that("result matches the manual log-mean-exp formula", {
  set.seed(42)
  loglik <- rnorm(300, -80, 3)
  m      <- max(loglik)
  manual <- m + log(mean(exp(loglik - m)))
  expect_equal(bsl_bridge_estimate(loglik), manual, tolerance = 1e-12)
})

test_that("works correctly for a length-1 input vector", {
  expect_equal(bsl_bridge_estimate(-42), -42, tolerance = 1e-10)
})
