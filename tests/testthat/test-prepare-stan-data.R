# Tests for prepare_stan_data_from_datasets() in R/utils.R
# Fixtures (.ds_range, .ds_iqr, .ds_meansd, .ds_freq4, .ds_freq5, .make_sd)
# are defined in helper-fixtures.R and loaded automatically by testthat.

# ---------------------------------------------------------------------------
# Summary type detection
# ---------------------------------------------------------------------------

test_that("type-1 (median + range) datasets are stored in obs_stat1/2/3", {
  sd <- .make_sd(.ds_range())
  expect_equal(as.integer(sd$summary_type), c(1L, 1L))
  expect_equal(as.vector(sd$obs_stat1), c(5, 7))    # medians
  expect_equal(as.vector(sd$obs_stat2), c(1, 2))    # mins
  expect_equal(as.vector(sd$obs_stat3), c(10, 15))  # maxs
})

test_that("type-2 (median + IQR) datasets are stored in obs_stat1/2/3", {
  sd <- .make_sd(.ds_iqr())
  expect_equal(as.integer(sd$summary_type), c(2L, 2L))
  expect_equal(as.vector(sd$obs_stat1), c(5, 6))
  expect_equal(as.vector(sd$obs_stat2), c(3, 4))
  expect_equal(as.vector(sd$obs_stat3), c(8, 9))
})

test_that("type-3 (mean + sd) datasets are stored in obs_stat1/2", {
  sd <- .make_sd(.ds_meansd())
  expect_equal(as.integer(sd$summary_type), c(3L, 3L))
  expect_equal(as.vector(sd$obs_stat1), c(5.0, 7.0))
  expect_equal(as.vector(sd$obs_stat2), c(2.0, 3.0))
})

test_that("type-4 (freq_value/freq_count) flat arrays and indexing are correct", {
  sd <- .make_sd(.ds_freq4())
  expect_equal(as.integer(sd$summary_type), c(4L, 4L))
  expect_equal(sd$n_freq_total,           9L)   # 5 + 4 entries
  expect_equal(as.integer(sd$freq_start), c(1L, 6L))
  expect_equal(as.integer(sd$freq_len),   c(5L, 4L))
  expect_equal(as.vector(sd$freq_value[1:5]), c(1, 2, 3, 4, 5))
  expect_equal(as.vector(sd$freq_value[6:9]), c(3, 4, 5, 6))
  # n: first dataset omits n → sum(freq_count) = 20; second provides n = 20
  expect_equal(as.integer(sd$n_obs), c(20L, 20L))
})

test_that("type-5 (interval-censored freq) flat arrays and indexing are correct", {
  sd <- .make_sd(.ds_freq5())
  expect_equal(as.integer(sd$summary_type), c(5L, 5L))
  expect_equal(sd$n_freq_total,           7L)   # 4 + 3 entries
  expect_equal(as.integer(sd$freq_start), c(1L, 5L))
  expect_equal(as.integer(sd$freq_len),   c(4L, 3L))
  expect_equal(as.vector(sd$freq_lower[1:4]), c(0, 3, 6, 9))
  expect_equal(as.vector(sd$freq_upper[1:4]), c(3, 6, 9, 12))
  # n: first dataset omits n → sum(freq_count) = 26; second provides n = 16
  expect_equal(as.integer(sd$n_obs), c(26L, 16L))
})

# ---------------------------------------------------------------------------
# freq_lower / freq_upper always present (backward compatibility)
# ---------------------------------------------------------------------------

test_that("freq_lower and freq_upper are present even when no freq-table datasets exist", {
  sd <- .make_sd(.ds_range())
  expect_false(is.null(sd$freq_lower))
  expect_false(is.null(sd$freq_upper))
  expect_equal(sd$n_freq_total,      0L)
  expect_equal(length(sd$freq_lower), 0L)
  expect_equal(length(sd$freq_upper), 0L)
})

# ---------------------------------------------------------------------------
# mu0 prior computation
# ---------------------------------------------------------------------------

test_that("mu0_mean is log of the mean of the medians (type 1)", {
  sd <- .make_sd(.ds_range())
  expect_equal(sd$mu0_mean, log(mean(c(5, 7))), tolerance = 1e-10)
})

test_that("mu0_mean uses the weighted mean of freq_value for type-4 datasets", {
  sd  <- .make_sd(.ds_freq4())
  wm1 <- sum(c(1,2,3,4,5) * c(2,5,8,4,1)) / sum(c(2,5,8,4,1))
  wm2 <- sum(c(3,4,5,6)   * c(3,6,5,2))   / sum(c(3,6,5,2))
  expect_equal(sd$mu0_mean, log(mean(c(wm1, wm2))), tolerance = 1e-10)
})

test_that("mu0_mean uses the interval midpoint mean for type-5 datasets", {
  sd   <- .make_sd(.ds_freq5())
  mid1 <- c(1.5, 4.5, 7.5, 10.5); cnt1 <- c(5, 10, 8, 3)
  mid2 <- c(2.5, 5.5, 8.5);        cnt2 <- c(4, 9,  3)
  wm1  <- sum(mid1 * cnt1) / sum(cnt1)
  wm2  <- sum(mid2 * cnt2) / sum(cnt2)
  expect_equal(sd$mu0_mean, log(mean(c(wm1, wm2))), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

test_that("errors when a dataset has no recognised summary statistic format", {
  ds <- list(d1 = list(foo = 5, bar = 2, n = 10))
  expect_error(
    suppressWarnings(prepare_stan_data_from_datasets(ds, dist_type = 1)),
    regexp = "does not have recognized"
  )
})

test_that("errors when freq_lower, freq_upper, freq_count have different lengths", {
  ds <- list(d1 = list(freq_lower = c(0, 3),
                       freq_upper = c(3, 6, 9),  # one extra
                       freq_count = c(5, 10)))
  expect_error(
    suppressWarnings(prepare_stan_data_from_datasets(ds, dist_type = 1)),
    regexp = "same length"
  )
})

test_that("errors when any freq_lower > its paired freq_upper", {
  ds <- list(d1 = list(freq_lower = c(0, 6),
                       freq_upper = c(3, 4),   # 6 > 4
                       freq_count = c(5, 10)))
  expect_error(
    suppressWarnings(prepare_stan_data_from_datasets(ds, dist_type = 1)),
    regexp = "freq_lower values must be <="
  )
})

# ---------------------------------------------------------------------------
# Mixed summary types
# ---------------------------------------------------------------------------

test_that("all five summary types can be mixed in a single call", {
  ds <- list(
    d1 = list(median = 5, min = 1, max = 10, n = 30),
    d2 = list(median = 6, Q1 = 4, Q3 = 9,   n = 40),
    d3 = list(mean   = 7, sd  = 2,           n = 25),
    d4 = list(freq_value = c(4, 5, 6), freq_count = c(3L, 7L, 2L)),
    d5 = list(freq_lower = c(0, 5), freq_upper = c(5, 10),
              freq_count = c(8L, 4L))
  )
  sd <- suppressWarnings(prepare_stan_data_from_datasets(ds, dist_type = 1))
  expect_equal(as.integer(sd$summary_type), c(1L, 2L, 3L, 4L, 5L))
  expect_equal(sd$n_datasets,   5L)
  expect_equal(sd$n_freq_total, 5L)   # 3 (d4) + 2 (d5)
})

# ---------------------------------------------------------------------------
# n_datasets < 5 warning
# ---------------------------------------------------------------------------

test_that("warns when n_datasets < 5", {
  expect_warning(
    prepare_stan_data_from_datasets(.ds_range(), dist_type = 1),
    regexp = "n_datasets"
  )
})

test_that("does not warn when n_datasets >= 5", {
  ds <- c(.ds_range(), .ds_iqr(), list(d5 = list(mean = 6, sd = 2, n = 30)))
  expect_silent(prepare_stan_data_from_datasets(ds, dist_type = 1))
})

# ---------------------------------------------------------------------------
# Default priors per dist_type
# ---------------------------------------------------------------------------

test_that("lognormal (dist_type=1) default: log_phi_mean=-0.7, log_phi_sd=0.5", {
  sd <- .make_sd(.ds_range(), dist_type = 1)
  expect_equal(sd$log_phi_mean, -0.7)
  expect_equal(sd$log_phi_sd,    0.5)
})

test_that("gamma (dist_type=2) default: log_phi_mean=2.5, log_phi_sd=0.5", {
  sd <- .make_sd(.ds_range(), dist_type = 2)
  expect_equal(sd$log_phi_mean, 2.5)
  expect_equal(sd$log_phi_sd,   0.5)
})

test_that("Weibull (dist_type=3) default: log_phi_mean=1.0, log_phi_sd=0.5", {
  sd <- .make_sd(.ds_range(), dist_type = 3)
  expect_equal(sd$log_phi_mean, 1.0)
  expect_equal(sd$log_phi_sd,   0.5)
})

# ---------------------------------------------------------------------------
# custom_priors override
# ---------------------------------------------------------------------------

test_that("custom_priors overrides individual fields while leaving others at defaults", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      .ds_range(), dist_type = 1,
      custom_priors = list(log_phi_mean = -0.3, mu0_sd = 2.0)
    )
  )
  expect_equal(sd$log_phi_mean, -0.3)
  expect_equal(sd$mu0_sd,        2.0)
  expect_equal(sd$log_phi_sd,    0.5)   # default unchanged
})
