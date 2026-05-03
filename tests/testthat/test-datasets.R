# Smoke tests for all built-in datasets.
# Checks: correct type, expected entry count, all recognised summary formats,
# and a round-trip through prepare_stan_data_from_datasets().

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.check_recognised <- function(ds, ds_name) {
  for (nm in names(ds)) {
    d <- ds[[nm]]
    has_range <- !is.null(d$median) && !is.null(d$min)  && !is.null(d$max)
    has_iqr   <- !is.null(d$median) && !is.null(d$Q1)   && !is.null(d$Q3)
    has_mean  <- !is.null(d$mean)   && !is.null(d$sd)
    has_freq4 <- !is.null(d$freq_value) && !is.null(d$freq_count)
    has_freq5 <- !is.null(d$freq_lower) && !is.null(d$freq_upper) &&
                 !is.null(d$freq_count)
    expect_true(has_range || has_iqr || has_mean || has_freq4 || has_freq5,
                info = paste(ds_name, nm, "has no recognised format"))
  }
}

.round_trip <- function(ds, expected_n, ds_name) {
  sd <- suppressWarnings(prepare_stan_data_from_datasets(ds, dist_type = 1))
  expect_equal(sd$n_datasets, expected_n,
               info = paste(ds_name, "n_datasets mismatch"))
  expect_true(all(sd$n_obs > 0),
              info = paste(ds_name, "has zero-n entry"))
}

# ---------------------------------------------------------------------------
# datasets_Nipah (11 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Nipah is a named list with 11 entries", {
  expect_type(datasets_Nipah, "list")
  expect_equal(length(datasets_Nipah), 11L)
})

test_that("every datasets_Nipah entry has a recognised summary-stat format", {
  .check_recognised(datasets_Nipah, "datasets_Nipah")
})

test_that("datasets_Nipah round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Nipah, 11L, "datasets_Nipah")
})

# ---------------------------------------------------------------------------
# datasets_MVD (2 entries)
# ---------------------------------------------------------------------------

test_that("datasets_MVD is a named list with 2 entries", {
  expect_type(datasets_MVD, "list")
  expect_equal(length(datasets_MVD), 2L)
})

test_that("every datasets_MVD entry has a recognised summary-stat format", {
  .check_recognised(datasets_MVD, "datasets_MVD")
})

test_that("datasets_MVD round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_MVD, 2L, "datasets_MVD")
})

# ---------------------------------------------------------------------------
# datasets_EVD (11 entries)
# ---------------------------------------------------------------------------

test_that("datasets_EVD is a named list with 11 entries", {
  expect_type(datasets_EVD, "list")
  expect_equal(length(datasets_EVD), 11L)
})

test_that("every datasets_EVD entry has a recognised summary-stat format", {
  .check_recognised(datasets_EVD, "datasets_EVD")
})

test_that("datasets_EVD round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_EVD, 11L, "datasets_EVD")
})

# ---------------------------------------------------------------------------
# datasets_SARS (22 entries)
# ---------------------------------------------------------------------------

test_that("datasets_SARS is a named list with 22 entries", {
  expect_type(datasets_SARS, "list")
  expect_equal(length(datasets_SARS), 22L)
})

test_that("every datasets_SARS entry has a recognised summary-stat format", {
  .check_recognised(datasets_SARS, "datasets_SARS")
})

test_that("datasets_SARS round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_SARS, 22L, "datasets_SARS")
})

# ---------------------------------------------------------------------------
# datasets_MERS (11 entries)
# ---------------------------------------------------------------------------

test_that("datasets_MERS is a named list with 11 entries", {
  expect_type(datasets_MERS, "list")
  expect_equal(length(datasets_MERS), 11L)
})

test_that("every datasets_MERS entry has a recognised summary-stat format", {
  .check_recognised(datasets_MERS, "datasets_MERS")
})

test_that("datasets_MERS round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_MERS, 11L, "datasets_MERS")
})

# ---------------------------------------------------------------------------
# datasets_Lassa (1 entry)
# ---------------------------------------------------------------------------

test_that("datasets_Lassa is a named list with 1 entry", {
  expect_type(datasets_Lassa, "list")
  expect_equal(length(datasets_Lassa), 1L)
})

test_that("every datasets_Lassa entry has a recognised summary-stat format", {
  .check_recognised(datasets_Lassa, "datasets_Lassa")
})

test_that("datasets_Lassa round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Lassa, 1L, "datasets_Lassa")
})

# ---------------------------------------------------------------------------
# datasets_Measles (12 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Measles is a named list with 12 entries", {
  expect_type(datasets_Measles, "list")
  expect_equal(length(datasets_Measles), 12L)
})

test_that("every datasets_Measles entry has a recognised summary-stat format", {
  .check_recognised(datasets_Measles, "datasets_Measles")
})

test_that("datasets_Measles round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Measles, 12L, "datasets_Measles")
})

# ---------------------------------------------------------------------------
# datasets_Mpox (16 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Mpox is a named list with 16 entries", {
  expect_type(datasets_Mpox, "list")
  expect_equal(length(datasets_Mpox), 16L)
})

test_that("every datasets_Mpox entry has a recognised summary-stat format", {
  .check_recognised(datasets_Mpox, "datasets_Mpox")
})

test_that("datasets_Mpox round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Mpox, 16L, "datasets_Mpox")
})

# ---------------------------------------------------------------------------
# datasets_Cholera (16 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Cholera is a named list with 16 entries", {
  expect_type(datasets_Cholera, "list")
  expect_equal(length(datasets_Cholera), 16L)
})

test_that("every datasets_Cholera entry has a recognised summary-stat format", {
  .check_recognised(datasets_Cholera, "datasets_Cholera")
})

test_that("datasets_Cholera round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Cholera, 16L, "datasets_Cholera")
})

# ---------------------------------------------------------------------------
# datasets_RVF (2 entries)
# ---------------------------------------------------------------------------

test_that("datasets_RVF is a named list with 2 entries", {
  expect_type(datasets_RVF, "list")
  expect_equal(length(datasets_RVF), 2L)
})

test_that("every datasets_RVF entry has a recognised summary-stat format", {
  .check_recognised(datasets_RVF, "datasets_RVF")
})

test_that("datasets_RVF round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_RVF, 2L, "datasets_RVF")
})

# ---------------------------------------------------------------------------
# datasets_CCHF (6 entries)
# ---------------------------------------------------------------------------

test_that("datasets_CCHF is a named list with 6 entries", {
  expect_type(datasets_CCHF, "list")
  expect_equal(length(datasets_CCHF), 6L)
})

test_that("every datasets_CCHF entry has a recognised summary-stat format", {
  .check_recognised(datasets_CCHF, "datasets_CCHF")
})

test_that("datasets_CCHF round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_CCHF, 6L, "datasets_CCHF")
})

# ---------------------------------------------------------------------------
# datasets_COVID_19 (51 entries)
# ---------------------------------------------------------------------------

test_that("datasets_COVID_19 is a named list with 51 entries", {
  expect_type(datasets_COVID_19, "list")
  expect_equal(length(datasets_COVID_19), 51L)
})

test_that("every datasets_COVID_19 entry has a recognised summary-stat format", {
  .check_recognised(datasets_COVID_19, "datasets_COVID_19")
})

test_that("datasets_COVID_19 round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_COVID_19, 51L, "datasets_COVID_19")
})

# ---------------------------------------------------------------------------
# datasets_Dengue (14 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Dengue is a named list with 14 entries", {
  expect_type(datasets_Dengue, "list")
  expect_equal(length(datasets_Dengue), 14L)
})

test_that("every datasets_Dengue entry has a recognised summary-stat format", {
  .check_recognised(datasets_Dengue, "datasets_Dengue")
})

test_that("datasets_Dengue round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Dengue, 14L, "datasets_Dengue")
})

# ---------------------------------------------------------------------------
# datasets_YFV (3 entries)
# ---------------------------------------------------------------------------

test_that("datasets_YFV is a named list with 3 entries", {
  expect_type(datasets_YFV, "list")
  expect_equal(length(datasets_YFV), 3L)
})

test_that("every datasets_YFV entry has a recognised summary-stat format", {
  .check_recognised(datasets_YFV, "datasets_YFV")
})

test_that("datasets_YFV round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_YFV, 3L, "datasets_YFV")
})

# ---------------------------------------------------------------------------
# datasets_flu (14 entries)
# ---------------------------------------------------------------------------

test_that("datasets_flu is a named list with 14 entries", {
  expect_type(datasets_flu, "list")
  expect_equal(length(datasets_flu), 14L)
})

test_that("every datasets_flu entry has a recognised summary-stat format", {
  .check_recognised(datasets_flu, "datasets_flu")
})

test_that("datasets_flu round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_flu, 14L, "datasets_flu")
})

# ---------------------------------------------------------------------------
# datasets_typhoid (22 entries)
# ---------------------------------------------------------------------------

test_that("datasets_typhoid is a named list with 22 entries", {
  expect_type(datasets_typhoid, "list")
  expect_equal(length(datasets_typhoid), 22L)
})

test_that("every datasets_typhoid entry has a recognised summary-stat format", {
  .check_recognised(datasets_typhoid, "datasets_typhoid")
})

test_that("datasets_typhoid round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_typhoid, 22L, "datasets_typhoid")
})

# ---------------------------------------------------------------------------
# datasets_Smallpox (4 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Smallpox is a named list with 4 entries", {
  expect_type(datasets_Smallpox, "list")
  expect_equal(length(datasets_Smallpox), 4L)
})

test_that("every datasets_Smallpox entry has a recognised summary-stat format", {
  .check_recognised(datasets_Smallpox, "datasets_Smallpox")
})

test_that("datasets_Smallpox round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Smallpox, 4L, "datasets_Smallpox")
})

# ---------------------------------------------------------------------------
# datasets_Zika (2 entries)
# ---------------------------------------------------------------------------

test_that("datasets_Zika is a named list with 2 entries", {
  expect_type(datasets_Zika, "list")
  expect_equal(length(datasets_Zika), 2L)
})

test_that("every datasets_Zika entry has a recognised summary-stat format", {
  .check_recognised(datasets_Zika, "datasets_Zika")
})

test_that("datasets_Zika round-trips through prepare_stan_data_from_datasets", {
  .round_trip(datasets_Zika, 2L, "datasets_Zika")
})

# ---------------------------------------------------------------------------
# Source fields: all datasets
# ---------------------------------------------------------------------------

test_that("every source field that is present is a non-empty character string", {
  all_ds <- list(
    Nipah    = datasets_Nipah,
    MVD      = datasets_MVD,
    EVD      = datasets_EVD,
    SARS     = datasets_SARS,
    MERS     = datasets_MERS,
    Lassa    = datasets_Lassa,
    Measles  = datasets_Measles,
    Mpox     = datasets_Mpox,
    Cholera  = datasets_Cholera,
    RVF      = datasets_RVF,
    CCHF     = datasets_CCHF,
    COVID_19 = datasets_COVID_19,
    Dengue   = datasets_Dengue,
    YFV      = datasets_YFV,
    flu      = datasets_flu,
    typhoid  = datasets_typhoid,
    Smallpox = datasets_Smallpox,
    Zika     = datasets_Zika
  )
  for (ds_name in names(all_ds)) {
    for (nm in names(all_ds[[ds_name]])) {
      src <- all_ds[[ds_name]][[nm]]$source
      if (!is.null(src)) {
        expect_type(src, "character")
        expect_gt(nchar(src), 0,
                  label = paste(ds_name, nm, "source"))
      }
    }
  }
})
