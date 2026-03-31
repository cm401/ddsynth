# Tests for Burr XII (dist_type = 4) and Generalised Gamma (dist_type = 5)
# extensions to the hierarchical data synthesis model.
#
# All tests operate at the R level — no Stan compilation or sampling is
# performed.  Coverage:
#   1. Burr XII CDF formula correctness
#   2. Generalised Gamma CDF formula correctness (via pgamma identity)
#   3. prepare_stan_data_from_datasets() produces correct Stan data lists
#      for dist_type 4 and 5
#   4. update_phi_prior() warns gracefully and returns data unchanged for
#      dist_type 4 and 5

# ---------------------------------------------------------------------------
# Inline CDF helpers — mirror what compute_predictive_cdf() does internally
# ---------------------------------------------------------------------------

.pburr12 <- function(x, lambda, c, k) {
  1 - (1 + (x / lambda)^c)^(-k)
}

.pgengamma <- function(x, mu, sigma, Q) {
  gamma_shape <- 1 / Q^2
  w <- (log(x) - mu) / sigma
  pgamma(gamma_shape * exp(Q * w), shape = gamma_shape, rate = 1)
}

# ---------------------------------------------------------------------------
# 1. Burr XII CDF
# ---------------------------------------------------------------------------

test_that("Burr XII CDF is 0 at x=0 and approaches 1 for large x", {
  expect_equal(.pburr12(0,   lambda = 1, c = 2, k = 3), 0)
  expect_gt(   .pburr12(1e6, lambda = 1, c = 2, k = 3), 0.9999)
})

test_that("Burr XII CDF is strictly increasing and in [0, 1]", {
  x    <- c(0.1, 0.5, 1, 2, 5, 10)
  vals <- .pburr12(x, lambda = 1, c = 2, k = 3)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_true(all(diff(vals) > 0))
})

test_that("Burr XII CDF matches closed-form at known evaluation points", {
  # F(x; lambda=1, c=2, k=3) = 1 - (1 + x^2)^(-3)
  expect_equal(.pburr12(1,   lambda = 1, c = 2, k = 3),
               1 - (1 + 1^2)^(-3),    tolerance = 1e-12)
  expect_equal(.pburr12(0.5, lambda = 1, c = 2, k = 3),
               1 - (1 + 0.5^2)^(-3),  tolerance = 1e-12)
  expect_equal(.pburr12(2,   lambda = 1, c = 2, k = 3),
               1 - (1 + 2^2)^(-3),    tolerance = 1e-12)
})

test_that("Burr XII CDF respects scale parameter", {
  # Scaling x by lambda is equivalent to using lambda=1 with x/lambda
  lambda <- 3; x <- 2; c_s <- 1.5; k_s <- 2
  expect_equal(.pburr12(x,          lambda = lambda, c = c_s, k = k_s),
               .pburr12(x / lambda, lambda = 1,      c = c_s, k = k_s),
               tolerance = 1e-12)
})

test_that("Burr XII CDF increases with k (heavier tail = slower rise)", {
  # Larger k => lighter right tail => CDF rises faster at moderate x
  x <- 2
  cdf_small_k <- .pburr12(x, lambda = 1, c = 2, k = 1)
  cdf_large_k <- .pburr12(x, lambda = 1, c = 2, k = 5)
  expect_gt(cdf_large_k, cdf_small_k)
})

# ---------------------------------------------------------------------------
# 2. Generalised Gamma CDF
# ---------------------------------------------------------------------------

test_that("GG CDF is close to 0 near 0 and approaches 1 for large x", {
  expect_lt(.pgengamma(1e-6, mu = 0, sigma = 1, Q = 1), 1e-3)
  expect_gt(.pgengamma(1e6,  mu = 0, sigma = 1, Q = 1), 0.9999)
})

test_that("GG CDF is strictly increasing and in [0, 1]", {
  x    <- c(0.5, 1, 2, 5, 10)
  vals <- .pgengamma(x, mu = 0, sigma = 1, Q = 1)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_true(all(diff(vals) > 0))
})

test_that("GG with Q=1, sigma=1 matches Exponential(rate=1)", {
  # GG(mu=0, sigma=1, Q=1) = Weibull(shape=1/sigma=1, scale=exp(0)=1)
  # = Exponential(rate=1)
  x_vals <- c(0.5, 1, 2, 3, 5)
  expect_equal(.pgengamma(x_vals, mu = 0, sigma = 1, Q = 1),
               pexp(x_vals, rate = 1),
               tolerance = 1e-10)
})

test_that("GG mu parameter acts as log-scale shift", {
  # GG(mu=log(2), sigma=1, Q=1) = Exp(rate=1/exp(log(2))) = Exp(rate=0.5)
  mu     <- log(2)
  x_vals <- c(0.5, 1, 2, 4)
  expect_equal(.pgengamma(x_vals, mu = mu, sigma = 1, Q = 1),
               pexp(x_vals, rate = 1 / exp(mu)),
               tolerance = 1e-10)
})

test_that("GG with very small Q approximates Lognormal", {
  # As Q -> 0+, GG(mu, sigma, Q) -> Lognormal(meanlog=mu, sdlog=sigma)
  x_vals <- c(1, 2, 5, 10, 20)
  expect_equal(.pgengamma(x_vals, mu = 1, sigma = 0.5, Q = 0.001),
               plnorm(x_vals, meanlog = 1, sdlog = 0.5),
               tolerance = 1e-3)
})

test_that("GG CDF uses pgamma identity correctly", {
  # F(x; mu, sigma, Q) = pgamma(gamma_shape * exp(Q*w), gamma_shape, 1)
  # where gamma_shape = 1/Q^2, w = (log(x)-mu)/sigma
  mu <- 1; sigma <- 0.5; Q <- 2; x <- 3
  gamma_shape <- 1 / Q^2
  w   <- (log(x) - mu) / sigma
  expected <- pgamma(gamma_shape * exp(Q * w), shape = gamma_shape, rate = 1)
  expect_equal(.pgengamma(x, mu = mu, sigma = sigma, Q = Q),
               expected, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3. prepare_stan_data_from_datasets() for dist_type 4 and 5
# ---------------------------------------------------------------------------

# Small but valid dataset list (2 studies with median+IQR)
.make_datasets <- function() {
  list(
    list(median = 5, Q1 = 3, Q3 = 8, n = 50),
    list(median = 6, Q1 = 4, Q3 = 9, n = 40)
  )
}

test_that("prepare_stan_data includes log_kappa fields for all dist_types 1-5", {
  for (dt in 1:5) {
    sd <- suppressWarnings(
      prepare_stan_data_from_datasets(.make_datasets(), dist_type = dt)
    )
    expect_true(!is.null(sd$log_kappa_mean),
                info = paste("log_kappa_mean missing for dist_type", dt))
    expect_true(!is.null(sd$log_kappa_sd),
                info = paste("log_kappa_sd missing for dist_type", dt))
    expect_true(is.numeric(sd$log_kappa_mean),
                info = paste("log_kappa_mean not numeric for dist_type", dt))
    expect_gt(sd$log_kappa_sd, 0,
              label = paste("log_kappa_sd for dist_type", dt))
  }
})

test_that("prepare_stan_data uses wide uninformative kappa prior for dist_type 1-3", {
  for (dt in 1:3) {
    sd <- suppressWarnings(
      prepare_stan_data_from_datasets(.make_datasets(), dist_type = dt)
    )
    expect_equal(sd$log_kappa_sd, 1.0,
                 info = paste("expected sd=1 for dist_type", dt))
  }
})

test_that("prepare_stan_data uses correct defaults for dist_type 4 (Burr XII)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 4)
  )
  expect_equal(sd$dist_type,      4L)
  expect_equal(sd$log_kappa_mean, 1.0)
  expect_equal(sd$log_kappa_sd,   0.5)
  expect_equal(sd$log_phi_mean,   0.7)   # c ~ 2, centre of typical range
})

test_that("prepare_stan_data uses correct defaults for dist_type 5 (GG)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 5)
  )
  expect_equal(sd$dist_type,      5L)
  expect_equal(sd$log_kappa_mean, 0.0)   # Q ~ 1, near Weibull
  expect_equal(sd$log_kappa_sd,   0.5)
  expect_equal(sd$log_phi_mean,  -0.5)   # sigma ~ 0.6
})

test_that("prepare_stan_data errors for dist_type outside 1-5", {
  expect_error(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 6),
    regexp = "dist_type"
  )
  expect_error(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 0),
    regexp = "dist_type"
  )
})

test_that("prepare_stan_data respects custom kappa priors via custom_priors", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      .make_datasets(), dist_type = 4,
      custom_priors = list(log_kappa_mean = 2.5, log_kappa_sd = 0.3)
    )
  )
  expect_equal(sd$log_kappa_mean, 2.5)
  expect_equal(sd$log_kappa_sd,   0.3)
})

test_that("prepare_stan_data dist_type 4 and 5 return valid Stan data structure", {
  for (dt in 4:5) {
    sd <- suppressWarnings(
      prepare_stan_data_from_datasets(.make_datasets(), dist_type = dt)
    )
    # Required fields present
    for (field in c("n_datasets", "n_obs", "summary_type", "dist_type",
                    "mu0_mean", "mu0_sd", "log_tau_mean", "log_tau_sd",
                    "log_phi_mean", "log_phi_sd",
                    "log_kappa_mean", "log_kappa_sd")) {
      expect_true(!is.null(sd[[field]]),
                  info = paste("field", field, "missing for dist_type", dt))
    }
    expect_equal(sd$n_datasets, 2L)
  }
})

# ---------------------------------------------------------------------------
# 4. update_phi_prior() works for dist_type 4/5 (no warning, updates log_phi_mean)
# ---------------------------------------------------------------------------

test_that("update_phi_prior emits no warning for dist_type 4 (Burr XII)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 4)
  )
  expect_no_warning(update_phi_prior(sd, .make_datasets()))
})

test_that("update_phi_prior changes log_phi_mean for dist_type 4 (Burr XII)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 4)
  )
  sd2 <- update_phi_prior(sd, .make_datasets())
  # log_phi_mean must change; log_kappa_mean is never touched
  expect_false(isTRUE(all.equal(sd2$log_phi_mean, sd$log_phi_mean)))
  expect_equal(sd2$log_kappa_mean, sd$log_kappa_mean)
})

test_that("update_phi_prior emits no warning for dist_type 5 (GG)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 5)
  )
  expect_no_warning(update_phi_prior(sd, .make_datasets()))
})

test_that("update_phi_prior changes log_phi_mean for dist_type 5 (GG)", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(.make_datasets(), dist_type = 5)
  )
  sd2 <- update_phi_prior(sd, .make_datasets())
  expect_false(isTRUE(all.equal(sd2$log_phi_mean, sd$log_phi_mean)))
})
