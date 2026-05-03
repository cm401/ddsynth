# Tests for should_attempt_gg() -----------------------------------------------

# Helpers ---------------------------------------------------------------------

# Frequency-table dataset (type-4): CV ≈ 0.24, "rich"
.gg_freq <- function(n = 30) {
  list(freq_value = c(5, 7, 9, 11, 13),
       freq_count = c(3, 8, 12, 8, 3), n = n)
}

# Interval-censored dataset (type-5): "rich"
.gg_icens <- function(n = 30) {
  list(freq_lower = c(4, 6, 8, 10, 12),
       freq_upper = c(6, 8, 10, 12, 14),
       freq_count = c(3, 8, 12, 8, 3), n = n)
}

# Range dataset (type-1): CV ≈ 0.22, not "rich"
.gg_range <- function(n = 30) {
  list(median = 9, min = 5, max = 13, n = n)
}

# IQR dataset (type-2): CV ≈ 0.30, not "rich"
.gg_iqr <- function(n = 30) {
  list(median = 9, Q1 = 6.3, Q3 = 11.7, n = n)
}

# Mean+SD dataset (type-3): not "rich"
.gg_meansd <- function(cv = 0.30, mean = 9, n = 30) {
  list(mean = mean, sd = mean * cv, n = n)
}

# ── 1. Returns TRUE when richness and spread are both acceptable ──────────────

test_that("returns TRUE with sufficient rich datasets and low CV spread", {
  ds <- list(d1 = .gg_freq(), d2 = .gg_freq(), d3 = .gg_range())
  # rich_frac = 2/3 ≈ 0.67 >= 0.30; all CVs similar => spread << 2.5
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

test_that("returns TRUE with only frequency-table datasets", {
  ds <- list(d1 = .gg_freq(), d2 = .gg_freq())
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

test_that("type-5 interval-censored datasets count as rich", {
  ds <- list(d1 = .gg_icens(), d2 = .gg_range())
  # rich_frac = 0.5 >= 0.30
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

test_that("returns TRUE for a single dataset (no spread check possible)", {
  ds <- list(d1 = .gg_freq())
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

# ── 2. Returns FALSE when rich fraction is too low ───────────────────────────

test_that("returns FALSE when no datasets are rich (all type-1/2/3)", {
  ds <- list(d1 = .gg_range(), d2 = .gg_range(), d3 = .gg_iqr(), d4 = .gg_meansd())
  # rich_frac = 0 < 0.30
  expect_false(should_attempt_gg(ds, verbose = FALSE))
})

test_that("returns FALSE when rich fraction is below custom threshold", {
  ds <- list(d1 = .gg_freq(), d2 = .gg_range(), d3 = .gg_range(), d4 = .gg_range())
  # rich_frac = 0.25; default min_rich_fraction = 0.30 → FALSE
  expect_false(should_attempt_gg(ds, verbose = FALSE))
})

test_that("returns TRUE when rich fraction exactly meets the threshold", {
  # 1 of 3 datasets is rich: 1/3 ≈ 0.33 >= 0.30
  ds <- list(d1 = .gg_freq(), d2 = .gg_range(), d3 = .gg_range())
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

test_that("custom min_rich_fraction threshold is respected", {
  ds <- list(d1 = .gg_freq(), d2 = .gg_range())  # rich_frac = 0.5
  expect_true( should_attempt_gg(ds, min_rich_fraction = 0.40, verbose = FALSE))
  expect_false(should_attempt_gg(ds, min_rich_fraction = 0.60, verbose = FALSE))
})

# ── 3. Returns FALSE when CV spread is too large ─────────────────────────────

test_that("returns FALSE when CV spread exceeds threshold", {
  # d1 CV ≈ 0.05 (concentrated), d2 CV ≈ 0.69 (dispersed) → spread ≈ 13.9
  d_conc <- list(freq_value = c(8, 9, 10), freq_count = c(1, 8, 1), n = 10)
  d_disp <- list(freq_value = c(2, 9, 20), freq_count = c(3, 4, 3), n = 10)
  ds <- list(d1 = d_conc, d2 = d_disp)
  expect_false(should_attempt_gg(ds, verbose = FALSE))
})

test_that("custom cv_spread_threshold is respected", {
  # spread ≈ 2; just straddles default of 2.5
  d1 <- list(freq_value = c(7, 9, 11), freq_count = c(2, 6, 2), n = 10)  # CV ≈ 0.19
  d2 <- list(freq_value = c(4, 9, 14), freq_count = c(2, 6, 2), n = 10)  # CV ≈ 0.38
  ds <- list(d1 = d1, d2 = d2)
  expect_true( should_attempt_gg(ds, cv_spread_threshold = 3.0, verbose = FALSE))
  expect_false(should_attempt_gg(ds, cv_spread_threshold = 1.5, verbose = FALSE))
})

# ── 4. CV spread check is skipped for a single valid CV ───────────────────────

test_that("spread check is skipped when only one dataset has a valid CV", {
  # Only one CV computable; no spread can be formed → proceed to richness check
  ds <- list(
    d1 = .gg_freq(),
    d2 = list(n = 10)  # no stat fields → NA CV
  )
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

# ── 5. Verbose messaging ──────────────────────────────────────────────────────

test_that("verbose=TRUE emits SKIP when richness is too low", {
  ds <- list(d1 = .gg_range(), d2 = .gg_range())
  expect_message(should_attempt_gg(ds, verbose = TRUE), "SKIP")
})

test_that("verbose=TRUE emits SKIP when CV spread is too large", {
  d_conc <- list(freq_value = c(8, 9, 10), freq_count = c(1, 8, 1), n = 10)
  d_disp <- list(freq_value = c(2, 9, 20), freq_count = c(3, 4, 3), n = 10)
  ds <- list(d1 = d_conc, d2 = d_disp)
  expect_message(should_attempt_gg(ds, verbose = TRUE), "SKIP")
})

test_that("verbose=TRUE emits OK when all checks pass", {
  ds <- list(d1 = .gg_freq(), d2 = .gg_freq())
  expect_message(should_attempt_gg(ds, verbose = TRUE), "OK")
})

test_that("verbose=FALSE suppresses all messages", {
  ds <- list(d1 = .gg_range(), d2 = .gg_range())
  expect_no_message(should_attempt_gg(ds, verbose = FALSE))
})

# ── 6. NA handling ────────────────────────────────────────────────────────────

test_that("datasets with no recognisable format produce NA CV without error", {
  ds <- list(
    d1 = .gg_freq(),
    d2 = .gg_freq(),
    d3 = list(n = 10)  # no stat fields
  )
  expect_true(should_attempt_gg(ds, verbose = FALSE))
})

test_that("all-NA CVs bypass the spread check and proceed to richness", {
  # All entries have unrecognised format → no valid CV → no spread check
  ds <- list(d1 = list(n = 10), d2 = list(n = 15))
  # rich_frac = 0 < 0.30 → FALSE (fails richness, not spread)
  expect_false(should_attempt_gg(ds, verbose = FALSE))
})
