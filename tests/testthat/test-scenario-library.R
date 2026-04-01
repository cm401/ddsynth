# Tests for generate_scenario_library() in R/utils.R
#
# Coverage:
#   1. Returned schema — all expected columns present with correct types
#   2. scenario_idx is unique, sequential, and starts at 1
#   3. Base (dist 1–3) defaults: kappa = 1, correct scenario_group
#   4. include_burr12 adds the right rows and column values
#   5. include_gengamma adds the right rows; include_gg_limitation is guarded
#   6. Row counts match expectations for each flag combination
#   7. summary_type integer is consistent with the proportion columns
#   8. vary_n flag is set correctly for varied-N scenarios

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.base_only <- function(...) {
  generate_scenario_library(
    include_homogeneous = TRUE,
    include_mixed       = TRUE,
    include_varied_n    = TRUE,
    include_freq_table  = FALSE,
    include_burr12      = FALSE,
    include_gengamma    = FALSE,
    ...
  )
}

# ---------------------------------------------------------------------------
# 1. Schema
# ---------------------------------------------------------------------------

test_that("generate_scenario_library returns a data frame with all required columns", {
  sc <- .base_only()
  required <- c("scenario_name", "scenario_group", "dist_type", "n_datasets",
                "mu0", "tau", "phi", "kappa",
                "n_obs_mean", "n_obs_sd", "n_obs_min", "n_obs_max",
                "summary_type_1_prop", "summary_type_2_prop",
                "summary_type_3_prop", "summary_type_4_prop",
                "summary_type", "vary_n", "scenario_idx")
  for (col in required) {
    expect_true(col %in% names(sc),
                info = paste("column missing:", col))
  }
})

test_that("column types are correct", {
  sc <- .base_only()
  expect_true(is.character(sc$scenario_name))
  expect_true(is.character(sc$scenario_group))
  expect_true(is.character(sc$dist_type))
  expect_true(is.numeric(sc$mu0))
  expect_true(is.numeric(sc$kappa))
  expect_true(is.numeric(sc$n_obs_mean))
  expect_true(is.integer(sc$summary_type))
  expect_true(is.logical(sc$vary_n))
  expect_true(is.integer(sc$scenario_idx))
})

# ---------------------------------------------------------------------------
# 2. scenario_idx integrity
# ---------------------------------------------------------------------------

test_that("scenario_idx is sequential starting at 1 with no gaps", {
  sc <- generate_scenario_library(
    include_burr12 = TRUE, include_gengamma = TRUE, include_gg_limitation = TRUE
  )
  expect_equal(sc$scenario_idx, seq_len(nrow(sc)))
})

test_that("scenario_idx is unique across all scenarios", {
  sc <- generate_scenario_library(include_burr12 = TRUE, include_gengamma = TRUE)
  expect_equal(length(unique(sc$scenario_idx)), nrow(sc))
})

test_that("scenario_name is unique across all scenarios", {
  sc <- generate_scenario_library(include_burr12 = TRUE, include_gengamma = TRUE)
  expect_equal(length(unique(sc$scenario_name)), nrow(sc))
})

# ---------------------------------------------------------------------------
# 3. Base scenarios (dist_type 1–3)
# ---------------------------------------------------------------------------

test_that("base scenarios have kappa = 1.0", {
  sc <- .base_only()
  expect_true(all(sc$kappa == 1.0))
})

test_that("base scenarios have correct scenario_group values", {
  sc <- .base_only()
  expected_groups <- paste0("base_", c("lognormal", "gamma", "weibull"))
  expect_true(all(sc$scenario_group %in% expected_groups))
})

test_that("dist_type values for base scenarios are lognormal, gamma, or weibull", {
  sc <- .base_only()
  expect_true(all(sc$dist_type %in% c("lognormal", "gamma", "weibull")))
})

# ---------------------------------------------------------------------------
# 4. Burr XII scenarios
# ---------------------------------------------------------------------------

test_that("include_burr12 = FALSE adds no Burr XII rows", {
  sc <- generate_scenario_library(include_burr12 = FALSE)
  expect_equal(sum(sc$scenario_group == "burr12"), 0L)
})

test_that("include_burr12 = TRUE adds Burr XII rows", {
  sc <- generate_scenario_library(include_burr12 = TRUE)
  n_burr <- sum(sc$scenario_group == "burr12")
  # Grid: 2 phi × 2 kappa × 3 D × 2 N = 24, plus 5 extras = 29
  expect_equal(n_burr, 29L)
})

test_that("Burr XII scenarios have dist_type = 'burr12'", {
  sc <- generate_scenario_library(include_burr12 = TRUE)
  burr <- sc[sc$scenario_group == "burr12", ]
  expect_true(all(burr$dist_type == "burr12"))
})

test_that("Burr XII scenarios have kappa > 1 (shape2 parameter)", {
  sc <- generate_scenario_library(include_burr12 = TRUE)
  burr <- sc[sc$scenario_group == "burr12", ]
  expect_true(all(burr$kappa >= 2.0))
})

test_that("Burr XII grid scenarios all use summary_type 1", {
  sc   <- generate_scenario_library(include_burr12 = TRUE)
  grid <- sc[grepl("^Burr12_c[0-9.]+_k[0-9.]+_D[0-9]+_N[0-9]+_ST1$",
                   sc$scenario_name), ]
  expect_true(all(grid$summary_type == 1L))
})

# ---------------------------------------------------------------------------
# 5. Generalised Gamma scenarios
# ---------------------------------------------------------------------------

test_that("include_gengamma = FALSE adds no GG rows", {
  sc <- generate_scenario_library(include_gengamma = FALSE)
  expect_equal(sum(sc$dist_type == "gengamma"), 0L)
})

test_that("include_gengamma = TRUE adds 6 standard GG rows", {
  sc <- generate_scenario_library(include_gengamma = TRUE,
                                  include_gg_limitation = FALSE)
  expect_equal(sum(sc$scenario_group == "gg_standard"),   6L)
  expect_equal(sum(sc$scenario_group == "gg_limitation"), 0L)
})

test_that("include_gg_limitation = TRUE requires include_gengamma = TRUE", {
  sc <- generate_scenario_library(include_gengamma     = TRUE,
                                  include_gg_limitation = TRUE)
  expect_equal(sum(sc$scenario_group == "gg_limitation"), 2L)
})

test_that("include_gg_limitation = TRUE without gengamma adds no limitation rows", {
  sc <- generate_scenario_library(include_gengamma     = FALSE,
                                  include_gg_limitation = TRUE)
  expect_equal(sum(sc$scenario_group == "gg_limitation"), 0L)
})

test_that("GG standard scenarios have dist_type = 'gengamma'", {
  sc <- generate_scenario_library(include_gengamma = TRUE)
  gg <- sc[sc$scenario_group == "gg_standard", ]
  expect_true(all(gg$dist_type == "gengamma"))
})

test_that("GG kappa values span the expected range (0.5 to 2.0 for standard)", {
  sc    <- generate_scenario_library(include_gengamma = TRUE)
  kvals <- sc$kappa[sc$scenario_group == "gg_standard"]
  expect_true(all(kvals %in% c(0.5, 1.0, 2.0)))
})

# ---------------------------------------------------------------------------
# 6. Row counts
# ---------------------------------------------------------------------------

test_that("full scenario library has 193 rows", {
  sc <- generate_scenario_library(
    include_homogeneous   = TRUE,
    include_mixed         = TRUE,
    include_varied_n      = TRUE,
    include_freq_table    = TRUE,
    include_burr12        = TRUE,
    include_gengamma      = TRUE,
    include_gg_limitation = TRUE
  )
  expect_equal(nrow(sc), 193L)
})

test_that("base-only library (no freq table) has consistent row count across dists", {
  sc <- generate_scenario_library(
    include_homogeneous = TRUE, include_mixed = TRUE,
    include_varied_n = TRUE, include_freq_table = FALSE
  )
  n_ln  <- sum(sc$dist_type == "lognormal")
  n_gam <- sum(sc$dist_type == "gamma")
  n_wei <- sum(sc$dist_type == "weibull")
  expect_equal(n_ln, n_gam)
  expect_equal(n_ln, n_wei)
})

# ---------------------------------------------------------------------------
# 7. summary_type consistency
# ---------------------------------------------------------------------------

test_that("summary_type = 1 iff summary_type_1_prop = 1", {
  sc <- generate_scenario_library(include_burr12 = TRUE, include_gengamma = TRUE)
  st1_rows <- sc[sc$summary_type == 1L, ]
  expect_true(all(st1_rows$summary_type_1_prop == 1))
  non_st1  <- sc[sc$summary_type != 1L, ]
  expect_true(all(non_st1$summary_type_1_prop < 1))
})

test_that("summary_type = 5 (mixed) iff no single prop equals 1", {
  sc   <- generate_scenario_library(include_mixed = TRUE)
  mixed <- sc[sc$summary_type == 5L, ]
  for (col in c("summary_type_1_prop", "summary_type_2_prop",
                "summary_type_3_prop", "summary_type_4_prop")) {
    expect_true(all(mixed[[col]] < 1),
                info = paste("mixed scenario has prop=1 in", col))
  }
})

test_that("summary_type proportions sum to 1 for every scenario", {
  sc <- generate_scenario_library(include_burr12 = TRUE, include_gengamma = TRUE)
  prop_sums <- rowSums(sc[, c("summary_type_1_prop", "summary_type_2_prop",
                               "summary_type_3_prop", "summary_type_4_prop")])
  expect_true(all(abs(prop_sums - 1) < 1e-10))
})

# ---------------------------------------------------------------------------
# 8. vary_n flag
# ---------------------------------------------------------------------------

test_that("vary_n = TRUE only for VarN scenarios", {
  sc <- generate_scenario_library(include_varied_n = TRUE)
  vary_rows   <- sc[sc$vary_n, ]
  stable_rows <- sc[!sc$vary_n, ]
  expect_true(all(grepl("VarN", vary_rows$scenario_name)))
  expect_true(all(!grepl("VarN", stable_rows$scenario_name)))
})

test_that("vary_n scenarios have finite n_obs_min and n_obs_max", {
  sc        <- generate_scenario_library(include_varied_n = TRUE,
                                         include_burr12   = TRUE)
  vary_rows <- sc[sc$vary_n, ]
  expect_true(all(is.finite(vary_rows$n_obs_min)))
  expect_true(all(is.finite(vary_rows$n_obs_max)))
})

test_that("non-vary_n scenarios have n_obs_min = n_obs_max = n_obs_mean for all dist types", {
  sc          <- generate_scenario_library(include_burr12 = TRUE,
                                           include_gengamma = TRUE,
                                           include_gg_limitation = TRUE)
  stable_rows <- sc[!sc$vary_n, ]
  expect_true(all(stable_rows$n_obs_min == stable_rows$n_obs_mean),
              info = "n_obs_min != n_obs_mean for some fixed-n scenario")
  expect_true(all(stable_rows$n_obs_max == stable_rows$n_obs_mean),
              info = "n_obs_max != n_obs_mean for some fixed-n scenario")
})
