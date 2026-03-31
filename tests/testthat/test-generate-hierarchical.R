# Tests for generate_hierarchical_data_mixed() in R/utils.R

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("returns a list with true_params and obs_data sub-lists", {
  set.seed(1)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 30, dist_type = "lognormal",
    mu0 = log(7), tau = 0.3, phi = 0.5, summary_type = 1L
  )
  expect_type(out, "list")
  expect_false(is.null(out$true_params))
  expect_false(is.null(out$obs_data))
})

test_that("true_params stores the supplied parameter values exactly", {
  set.seed(2)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 30, dist_type = "lognormal",
    mu0 = log(7), tau = 0.3, phi = 0.5, summary_type = 1L
  )
  expect_equal(out$true_params$mu0, log(7))
  expect_equal(out$true_params$tau, 0.3)
  expect_equal(out$true_params$phi, 0.5)
})

# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------

test_that("all three distributions produce positive obs_stat1 values (type-1)", {
  set.seed(11)
  for (dist in c("lognormal", "gamma", "weibull")) {
    out <- ddsynth:::generate_hierarchical_data_mixed(
      n_datasets = 5, n_obs = 20, dist_type = dist,
      mu0 = log(5), tau = 0.2, phi = 1.5, summary_type = 1L
    )
    expect_true(all(as.vector(out$obs_data$obs_stat1) > 0),
                info = paste("Non-positive median for", dist))
  }
})

# ---------------------------------------------------------------------------
# n_obs handling
# ---------------------------------------------------------------------------

test_that("scalar n_obs is replicated to length n_datasets", {
  set.seed(22)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 6, n_obs = 50, dist_type = "lognormal",
    mu0 = log(6), tau = 0.2, phi = 0.4, summary_type = 1L
  )
  n_obs_vec <- as.vector(out$obs_data$n_obs)
  expect_length(n_obs_vec, 6L)
  expect_true(all(n_obs_vec == 50L))
})

test_that("errors when n_obs length does not match n_datasets", {
  expect_error(
    ddsynth:::generate_hierarchical_data_mixed(
      n_datasets = 5, n_obs = c(30, 40),
      dist_type = "lognormal", mu0 = 0, tau = 0.2, phi = 0.5
    ),
    regexp = "n_obs"
  )
})

# ---------------------------------------------------------------------------
# Summary type output invariants
# ---------------------------------------------------------------------------

test_that("type-1: min <= median <= max for every dataset", {
  set.seed(33)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 100, dist_type = "lognormal",
    mu0 = log(5), tau = 0.1, phi = 0.3, summary_type = 1L
  )
  s1 <- as.vector(out$obs_data$obs_stat1)
  s2 <- as.vector(out$obs_data$obs_stat2)
  s3 <- as.vector(out$obs_data$obs_stat3)
  expect_true(all(s2 <= s1))
  expect_true(all(s1 <= s3))
})

test_that("type-2: Q1 <= median <= Q3 for every dataset", {
  set.seed(44)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 100, dist_type = "lognormal",
    mu0 = log(5), tau = 0.1, phi = 0.3, summary_type = 2L
  )
  s1 <- as.vector(out$obs_data$obs_stat1)
  s2 <- as.vector(out$obs_data$obs_stat2)
  s3 <- as.vector(out$obs_data$obs_stat3)
  expect_true(all(s2 <= s1))
  expect_true(all(s1 <= s3))
})

test_that("type-3: mean and sd are both positive for every dataset", {
  set.seed(55)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 100, dist_type = "gamma",
    mu0 = log(5), tau = 0.1, phi = 2.0, summary_type = 3L
  )
  s1 <- as.vector(out$obs_data$obs_stat1)
  s2 <- as.vector(out$obs_data$obs_stat2)
  expect_true(all(s1 > 0))
  expect_true(all(s2 > 0))
})

test_that("type-4: frequency tables are populated and n_freq_total is positive", {
  set.seed(66)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 5, n_obs = 50, dist_type = "lognormal",
    mu0 = log(5), tau = 0.2, phi = 0.4, summary_type = 4L
  )
  expect_true(out$obs_data$n_freq_total > 0L)
  expect_true(out$obs_data$n_freq_total <= 5 * 50)
})

# ---------------------------------------------------------------------------
# Probability-vector summary_type dispatch
# ---------------------------------------------------------------------------

test_that("probability-vector summary_type assigns types 1–4 across datasets", {
  set.seed(77)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 20, n_obs = 50, dist_type = "lognormal",
    mu0 = log(5), tau = 0.2, phi = 0.4,
    summary_type = c(0.25, 0.25, 0.25, 0.25)
  )
  st <- as.vector(out$obs_data$summary_type)
  expect_length(st, 20L)
  expect_true(all(st %in% 1:4))
})

# ---------------------------------------------------------------------------
# Accuracy
# ---------------------------------------------------------------------------

test_that("obs_stat1 (type-3 mean) is centred near exp(mu0) for large samples", {
  set.seed(88)
  out <- ddsynth:::generate_hierarchical_data_mixed(
    n_datasets = 20, n_obs = 500, dist_type = "lognormal",
    mu0 = log(7), tau = 0.05, phi = 0.3, summary_type = 3L
  )
  expect_equal(median(as.vector(out$obs_data$obs_stat1)), 7.0, tolerance = 0.5)
})

# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

test_that("errors with 'arg' message for an unknown distribution name", {
  expect_error(
    ddsynth:::generate_hierarchical_data_mixed(
      n_datasets = 5, n_obs = 30, dist_type = "pareto",
      mu0 = 0, tau = 0.2, phi = 0.5
    ),
    regexp = "arg"
  )
})
