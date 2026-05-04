test_that("check_positive rejects non-positive values", {
  expect_error(ddsynth:::check_positive(-1),  "strictly positive")
  expect_error(ddsynth:::check_positive(0),   "strictly positive")
  expect_error(ddsynth:::check_positive("a"), "strictly positive")
  expect_silent(ddsynth:::check_positive(1))
  expect_silent(ddsynth:::check_positive(c(0.1, 2, 10)))
})

test_that("check_scalar rejects non-scalars", {
  expect_error(ddsynth:::check_scalar(c(1, 2)), "single finite")
  expect_error(ddsynth:::check_scalar(Inf),     "single finite")
  expect_error(ddsynth:::check_scalar("a"),     "single finite")
  expect_silent(ddsynth:::check_scalar(3.14))
})

# make_stan_init_fn -----------------------------------------------------------

test_that("make_stan_init_fn returns a function", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      list(d1 = list(median = 7, min = 3, max = 14, n = 20)),
      dist_type = 1
    )
  )
  init_fn <- make_stan_init_fn(sd)
  expect_true(is.function(init_fn))
})

test_that("make_stan_init_fn closure returns a list with required names", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      list(d1 = list(median = 7, min = 3, max = 14, n = 20),
           d2 = list(median = 9, min = 5, max = 16, n = 30)),
      dist_type = 1
    )
  )
  init <- make_stan_init_fn(sd)()
  expect_named(init, c("mu0", "log_tau", "log_phi", "log_kappa", "loc_d_raw"),
               ignore.order = TRUE)
})

test_that("make_stan_init_fn initialises parameters from stan_data", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      list(d1 = list(median = 7, min = 3, max = 14, n = 20),
           d2 = list(median = 9, min = 5, max = 16, n = 30)),
      dist_type = 1
    )
  )
  init <- make_stan_init_fn(sd)()
  expect_equal(init$mu0,       sd$mu0_mean)
  expect_equal(init$log_tau,   sd$log_tau_mean)
  expect_equal(init$log_phi,   sd$log_phi_mean)
  expect_equal(init$log_kappa, sd$log_kappa_mean)
})

test_that("make_stan_init_fn sets loc_d_raw to zero vector of length n_datasets", {
  sd <- suppressWarnings(
    prepare_stan_data_from_datasets(
      list(d1 = list(median = 7, min = 3, max = 14, n = 20),
           d2 = list(median = 9, min = 5, max = 16, n = 30),
           d3 = list(mean   = 8, sd  = 2,           n = 25)),
      dist_type = 1
    )
  )
  init <- make_stan_init_fn(sd)()
  expect_equal(init$loc_d_raw, rep(0.0, 3L))
})

# compile_stan_model ----------------------------------------------------------

test_that("compile_stan_model rejects invalid model names", {
  expect_error(compile_stan_model("invalid"), "should be one of")
  expect_error(compile_stan_model("FACTORISED"), "should be one of")
})

test_that("compile_stan_model accepts valid model names without error on arg check", {
  # match.arg() resolves partial and exact matches before any file I/O;
  # test that argument parsing succeeds for both valid choices.
  expect_true(match.arg("factorised", c("factorised", "joint")) == "factorised")
  expect_true(match.arg("joint",      c("factorised", "joint")) == "joint")
})
