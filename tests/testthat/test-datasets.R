# Smoke tests for built-in datasets: datasets_Nipah, datasets_MVD, datasets_EVD
# Checks structure, recognised formats, and round-trip through prepare_stan_data_from_datasets().

# ---------------------------------------------------------------------------
# Helper: verify every entry has at least one recognised summary-stat combination
# ---------------------------------------------------------------------------
.check_recognised <- function(ds, ds_name) {
  for (nm in names(ds)) {
    d <- ds[[nm]]
    has_range <- !is.null(d$median) && !is.null(d$min) && !is.null(d$max)
    has_iqr   <- !is.null(d$median) && !is.null(d$Q1)  && !is.null(d$Q3)
    has_mean  <- !is.null(d$mean)   && !is.null(d$sd)
    has_freq4 <- !is.null(d$freq_value) && !is.null(d$freq_count)
    has_freq5 <- !is.null(d$freq_lower) && !is.null(d$freq_upper) && !is.null(d$freq_count)
    expect_true(has_range || has_iqr || has_mean || has_freq4 || has_freq5,
                info = paste(ds_name, nm, "has no recognised format"))
  }
}

# ---------------------------------------------------------------------------
# datasets_Nipah
# ---------------------------------------------------------------------------

test_that("datasets_Nipah is a named list with 11 entries", {
  expect_type(datasets_Nipah, "list")
  expect_equal(length(datasets_Nipah), 11L)
})

test_that("every datasets_Nipah entry has a recognised summary-stat format", {
  .check_recognised(datasets_Nipah, "datasets_Nipah")
})

test_that("datasets_Nipah round-trips through prepare_stan_data_from_datasets", {
  sd <- prepare_stan_data_from_datasets(datasets_Nipah, dist_type = 1)
  expect_equal(sd$n_datasets, 11L)
  expect_true(all(sd$n_obs > 0))
})

# ---------------------------------------------------------------------------
# datasets_MVD
# ---------------------------------------------------------------------------

test_that("datasets_MVD is a named list with 2 entries", {
  expect_type(datasets_MVD, "list")
  expect_equal(length(datasets_MVD), 2L)
})

test_that("every datasets_MVD entry has a recognised summary-stat format", {
  .check_recognised(datasets_MVD, "datasets_MVD")
})

test_that("datasets_MVD round-trips through prepare_stan_data_from_datasets", {
  sd <- suppressWarnings(prepare_stan_data_from_datasets(datasets_MVD, dist_type = 1))
  expect_equal(sd$n_datasets, 2L)
  expect_true(all(sd$n_obs > 0))
})

# ---------------------------------------------------------------------------
# datasets_EVD
# ---------------------------------------------------------------------------

test_that("datasets_EVD is a named list with 13 entries", {
  expect_type(datasets_EVD, "list")
  expect_equal(length(datasets_EVD), 11L)
})

test_that("every datasets_EVD entry has a recognised summary-stat format", {
  .check_recognised(datasets_EVD, "datasets_EVD")
})

test_that("datasets_EVD round-trips through prepare_stan_data_from_datasets", {
  sd <- prepare_stan_data_from_datasets(datasets_EVD, dist_type = 1)
  expect_equal(sd$n_datasets, 11L)
  expect_true(all(sd$n_obs > 0))
})

# ---------------------------------------------------------------------------
# source fields (all datasets)
# ---------------------------------------------------------------------------

test_that("every source field that is present is a non-empty character string", {
  all_ds <- list(Nipah = datasets_Nipah, MVD = datasets_MVD, EVD = datasets_EVD)
  for (ds_name in names(all_ds)) {
    for (nm in names(all_ds[[ds_name]])) {
      src <- all_ds[[ds_name]][[nm]]$source
      if (!is.null(src)) {
        expect_type(src, "character")
        expect_gt(nchar(src), 0)
      }
    }
  }
})
