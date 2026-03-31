# Tests for update_phi_prior() in R/utils.R
# Fixtures (.make_sd, .five_meansd, .ds_freq4, .ds_freq5) are in helper-fixtures.R.

# ---------------------------------------------------------------------------
# Lognormal (dist_type = 1): phi = sqrt(log(1 + (sd/mean)^2))  (sdlog, not variance)
# ---------------------------------------------------------------------------

test_that("lognormal: log_phi_mean = log(median(sqrt(log(1 + (sd/mean)^2))))", {
  means  <- c(5, 6, 7, 8, 9)
  sds    <- c(2, 3, 2, 4, 3)
  ds     <- .five_meansd(means, sds)
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  implied  <- sqrt(log(1 + (sds / means)^2))
  expect_equal(result$log_phi_mean, log(median(implied)), tolerance = 1e-10)
})

test_that("lognormal: only log_phi_mean is changed; other fields are unmodified", {
  ds     <- .five_meansd(c(5,6,7,8,9), c(2,3,2,4,3))
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  expect_equal(result$log_phi_sd,   sd_obj$log_phi_sd)
  expect_equal(result$log_tau_mean, sd_obj$log_tau_mean)
  expect_equal(result$n_datasets,   sd_obj$n_datasets)
})

# ---------------------------------------------------------------------------
# Gamma (dist_type = 2): phi = (mean/sd)^2
# ---------------------------------------------------------------------------

test_that("gamma: log_phi_mean = log(median((mean/sd)^2))", {
  means  <- c(5, 6, 7, 8, 9)
  sds    <- c(2, 3, 2, 4, 3)
  ds     <- .five_meansd(means, sds)
  sd_obj <- .make_sd(ds, dist_type = 2)
  result <- update_phi_prior(sd_obj, ds)

  implied <- (means / sds)^2
  expect_equal(result$log_phi_mean, log(median(implied)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Weibull (dist_type = 3): phi solved from CV via uniroot
# ---------------------------------------------------------------------------

test_that("weibull: recovers shape k from data with CV matching known k=2", {
  # CV for Weibull(shape=k): sqrt(gamma(1+2/k)/gamma(1+1/k)^2 - 1)
  k_true <- 2
  cv     <- sqrt(gamma(1 + 2/k_true) / gamma(1 + 1/k_true)^2 - 1)
  ds     <- .five_meansd(rep(5, 5), rep(5 * cv, 5))
  sd_obj <- .make_sd(ds, dist_type = 3)
  result <- update_phi_prior(sd_obj, ds)

  expect_equal(exp(result$log_phi_mean), k_true, tolerance = 0.01)
})

# ---------------------------------------------------------------------------
# All five summary-type inputs
# ---------------------------------------------------------------------------

test_that("works with IQR datasets: sd approximated as (Q3-Q1)/1.35", {
  ds     <- lapply(1:5, function(i) list(median = 7, Q1 = 5, Q3 = 9, n = 30))
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  sd_est  <- (9 - 5) / 1.35
  implied <- sqrt(log(1 + (sd_est / 7)^2))
  expect_equal(result$log_phi_mean, log(median(rep(implied, 5))), tolerance = 1e-10)
})

test_that("works with range datasets: sd approximated as (max-min)/4", {
  ds     <- lapply(1:5, function(i) list(median = 6, min = 1, max = 13, n = 25))
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  sd_est  <- (13 - 1) / 4
  implied <- sqrt(log(1 + (sd_est / 6)^2))
  expect_equal(result$log_phi_mean, log(median(rep(implied, 5))), tolerance = 1e-10)
})

test_that("works with freq_value datasets: uses weighted SD", {
  fv     <- c(3, 4, 5, 6, 7)
  fc     <- c(4L, 8L, 12L, 8L, 4L)
  ds     <- lapply(1:5, function(i) list(freq_value = fv, freq_count = fc))
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  w        <- fc / sum(fc)
  mean_est <- sum(fv * w)
  sd_est   <- sqrt(sum(w * (fv - mean_est)^2))
  implied  <- sqrt(log(1 + (sd_est / mean_est)^2))
  expect_equal(result$log_phi_mean, log(median(rep(implied, 5))), tolerance = 1e-10)
})

test_that("works with interval-censored datasets: uses midpoint weighted SD", {
  fl     <- c(0, 5, 10); fu <- c(5, 10, 15); fc <- c(5L, 10L, 5L)
  ds     <- lapply(1:5, function(i) list(freq_lower = fl, freq_upper = fu, freq_count = fc))
  sd_obj <- .make_sd(ds, dist_type = 1)
  result <- update_phi_prior(sd_obj, ds)

  mid      <- (fl + fu) / 2
  w        <- fc / sum(fc)
  mean_est <- sum(mid * w)
  sd_est   <- sqrt(sum(w * (mid - mean_est)^2))
  implied  <- sqrt(log(1 + (sd_est / mean_est)^2))
  expect_equal(result$log_phi_mean, log(implied), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Generalised Gamma (dist_type = 5): phi = sigma solved from CV with kappa fixed
# ---------------------------------------------------------------------------

test_that("gengamma: no warning is emitted", {
  ds     <- .five_meansd(c(5, 6, 7, 8, 9), c(2, 3, 2, 4, 3))
  sd_obj <- .make_sd(ds, dist_type = 5)
  expect_no_warning(update_phi_prior(sd_obj, ds))
})

test_that("gengamma: recovers sigma from CV for known (sigma, kappa) at kappa=1", {
  # With kappa=1 (default GG prior), gamma_shape=1, and the CV formula reduces to
  # CV^2 = Gamma(1+2*sigma)/Gamma(1+sigma)^2 - 1.
  # For sigma=0.5: CV = sqrt(Gamma(2)/Gamma(1.5)^2 - 1)
  sigma_true <- 0.5
  kappa_val  <- 1.0                   # exp(log_kappa_mean=0)
  gs         <- 1.0 / kappa_val^2     # = 1
  cv_true    <- sqrt(exp(lgamma(gs + 2*sigma_true/kappa_val) +
                         lgamma(gs) -
                         2*lgamma(gs + sigma_true/kappa_val)) - 1)
  ds     <- .five_meansd(rep(5, 5), rep(5 * cv_true, 5))
  sd_obj <- .make_sd(ds, dist_type = 5)   # log_kappa_mean = 0 → kappa=1
  result <- update_phi_prior(sd_obj, ds)

  expect_equal(exp(result$log_phi_mean), sigma_true, tolerance = 0.01)
})

test_that("gengamma: log_phi_mean increases as CV increases", {
  # Higher CV → higher sigma needed
  ds_low  <- .five_meansd(rep(5, 5), rep(5 * 0.3, 5))
  ds_high <- .five_meansd(rep(5, 5), rep(5 * 0.8, 5))
  sd_obj  <- .make_sd(ds_low, dist_type = 5)
  r_low   <- update_phi_prior(sd_obj,                    ds_low)
  r_high  <- update_phi_prior(.make_sd(ds_high, 5), ds_high)

  expect_lt(r_low$log_phi_mean, r_high$log_phi_mean)
})

# ---------------------------------------------------------------------------
# Burr XII (dist_type = 4): phi = c solved from CV with kappa fixed
# ---------------------------------------------------------------------------

test_that("burr12: no warning is emitted", {
  ds     <- .five_meansd(c(5, 6, 7, 8, 9), c(2, 3, 2, 4, 3))
  sd_obj <- .make_sd(ds, dist_type = 4)
  expect_no_warning(update_phi_prior(sd_obj, ds))
})

test_that("burr12: recovers c from CV for known (c, k) at k=exp(1)", {
  # With kappa=exp(log_kappa_mean=1)=e, and known c (phi), compute the expected CV.
  c_true    <- 3.0
  kappa_val <- exp(1)               # default Burr XII prior mean
  lb1 <- lbeta(kappa_val - 1/c_true, 1 + 1/c_true)
  lb2 <- lbeta(kappa_val - 2/c_true, 1 + 2/c_true)
  cv_true <- sqrt(exp(lb2 - log(kappa_val) - 2*lb1) - 1)

  ds     <- .five_meansd(rep(5, 5), rep(5 * cv_true, 5))
  sd_obj <- .make_sd(ds, dist_type = 4)   # log_kappa_mean = 1 → kappa=e
  result <- update_phi_prior(sd_obj, ds)

  expect_equal(exp(result$log_phi_mean), c_true, tolerance = 0.01)
})

test_that("burr12: log_phi_mean decreases as CV increases (CV and c are inversely related)", {
  ds_low  <- .five_meansd(rep(5, 5), rep(5 * 0.3, 5))
  ds_high <- .five_meansd(rep(5, 5), rep(5 * 0.8, 5))
  r_low   <- update_phi_prior(.make_sd(ds_low,  4), ds_low)
  r_high  <- update_phi_prior(.make_sd(ds_high, 4), ds_high)

  expect_gt(r_low$log_phi_mean, r_high$log_phi_mean)
})

# ---------------------------------------------------------------------------
# Edge case: no finite phi can be derived
# ---------------------------------------------------------------------------

test_that("warns and returns stan_data unchanged when no finite phi is derivable", {
  # Datasets with no recognised format → vapply returns all NAs → median = NA
  ds_unrecognised <- lapply(1:5, function(i) list(foo = 5, n = 10))
  sd_obj <- .make_sd(
    lapply(1:5, function(i) list(median = 5, min = 1, max = 10, n = 20)),
    dist_type = 1
  )
  phi_before <- sd_obj$log_phi_mean

  expect_warning(
    update_phi_prior(sd_obj, ds_unrecognised),
    regexp = "could not derive"
  )
  result <- suppressWarnings(update_phi_prior(sd_obj, ds_unrecognised))
  expect_equal(result$log_phi_mean, phi_before)
})
