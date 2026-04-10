# Tests for gamma_type2_reliable() -------------------------------------------

# Helper: build a type-2 dataset list
.t2 <- function(median, Q1, Q3, n = 100) {
  list(median = median, Q1 = Q1, Q3 = Q3, n = n)
}

# Helper: build a type-1 dataset (median + range)
.t1 <- function(median, min, max, n = 50) {
  list(median = median, min = min, max = max, n = n)
}

# Helper: build a type-3 dataset (mean + sd)
.t3 <- function(mean, sd, n = 50) {
  list(mean = mean, sd = sd, n = n)
}

# ── 1. No type-2 datasets ─────────────────────────────────────────────────────

test_that("returns TRUE silently when no type-2 datasets present", {
  ds <- list(d1 = .t1(10, 2, 20), d2 = .t3(10, 3))
  expect_true(gamma_type2_reliable(ds, verbose = FALSE))
})

# ── 2. High implied shape → FALSE ─────────────────────────────────────────────

test_that("returns FALSE when median implied shape exceeds max_implied_shape", {
  # IQR/median = 1.6/10 = 0.16 → phi_impl = (1.35/0.16)^2 ≈ 71
  ds <- list(
    d1 = .t2(median = 10, Q1 = 9.2, Q3 = 10.8),
    d2 = .t2(median = 12, Q1 = 11.2, Q3 = 12.8)
  )
  expect_false(gamma_type2_reliable(ds, verbose = FALSE))
})

test_that("returns FALSE with default threshold for concentrated datasets", {
  # phi_impl ≈ (1.35 * 5 / 0.5)^2 = (13.5)^2 = 182
  ds <- list(d1 = .t2(median = 5, Q1 = 4.75, Q3 = 5.25))
  expect_false(gamma_type2_reliable(ds, verbose = FALSE))
})

# ── 3. Low implied shape → TRUE ───────────────────────────────────────────────

test_that("returns TRUE when implied shape is below threshold", {
  # IQR/median = 8/10 = 0.8 → phi_impl = (1.35/0.8)^2 ≈ 2.85
  ds <- list(
    d1 = .t2(median = 10, Q1 = 6, Q3 = 14),
    d2 = .t2(median = 8,  Q1 = 4, Q3 = 12)
  )
  expect_true(gamma_type2_reliable(ds, verbose = FALSE))
})

test_that("returns TRUE at exactly the boundary (implied shape = max_implied_shape)", {
  # phi_impl = 20 → IQR/median = 1.35/sqrt(20) ≈ 0.302
  # IQR ≈ 0.302 * 10 = 3.02
  ds <- list(d1 = .t2(median = 10, Q1 = 10 - 1.51, Q3 = 10 + 1.51))
  # Should be TRUE (just at threshold, not exceeding)
  result <- gamma_type2_reliable(ds, max_implied_shape = 20, verbose = FALSE)
  expect_true(result)
})

# ── 4. max_implied_shape threshold respected ──────────────────────────────────

test_that("custom max_implied_shape changes the decision", {
  # phi_impl ≈ 8 (moderate)
  ds <- list(d1 = .t2(median = 10, Q1 = 7.6, Q3 = 12.4))
  expect_true( gamma_type2_reliable(ds, max_implied_shape = 20, verbose = FALSE))
  expect_false(gamma_type2_reliable(ds, max_implied_shape = 5,  verbose = FALSE))
})

# ── 5. Verbose messages ───────────────────────────────────────────────────────

test_that("verbose=TRUE emits a SKIP message for high implied shape", {
  ds <- list(d1 = .t2(median = 10, Q1 = 9.5, Q3 = 10.5))
  expect_message(
    gamma_type2_reliable(ds, verbose = TRUE),
    "SKIP"
  )
})

test_that("verbose=TRUE emits an OK message when check passes", {
  ds <- list(d1 = .t2(median = 10, Q1 = 6, Q3 = 14))
  expect_message(
    gamma_type2_reliable(ds, verbose = TRUE),
    "OK"
  )
})

test_that("verbose=FALSE suppresses all messages", {
  ds <- list(d1 = .t2(median = 10, Q1 = 9.5, Q3 = 10.5))
  expect_no_message(gamma_type2_reliable(ds, verbose = FALSE))
})

# ── 6. Small n advisory ───────────────────────────────────────────────────────

test_that("advisory NOTE fires when most type-2 datasets have n < min_n", {
  ds <- list(
    d1 = .t2(median = 10, Q1 = 6, Q3 = 14, n = 20),
    d2 = .t2(median = 8,  Q1 = 5, Q3 = 12, n = 15)
  )
  expect_message(
    gamma_type2_reliable(ds, min_n = 50, verbose = TRUE),
    "NOTE"
  )
})

test_that("advisory NOTE does not cause FALSE return", {
  ds <- list(
    d1 = .t2(median = 10, Q1 = 6, Q3 = 14, n = 20),
    d2 = .t2(median = 8,  Q1 = 5, Q3 = 12, n = 15)
  )
  # Implied shape is low (OK), but n is small — should still return TRUE
  result <- suppressMessages(gamma_type2_reliable(ds, min_n = 50, verbose = TRUE))
  expect_true(result)
})

# ── 7. NA handling ────────────────────────────────────────────────────────────

test_that("datasets with IQR = 0 are skipped without error", {
  ds <- list(
    d1 = .t2(median = 10, Q1 = 10, Q3 = 10),   # IQR = 0 → NA
    d2 = .t2(median = 10, Q1 = 6,  Q3 = 14)    # phi_impl ≈ 2.85
  )
  expect_true(gamma_type2_reliable(ds, verbose = FALSE))
})

test_that("all-NA implied shapes do not cause error and return TRUE", {
  # All have IQR = 0 — no valid implied shape to compare
  ds <- list(
    d1 = .t2(median = 10, Q1 = 10, Q3 = 10),
    d2 = .t2(median = 8,  Q1 = 8,  Q3 = 8)
  )
  expect_true(gamma_type2_reliable(ds, verbose = FALSE))
})

# ── 8. Mixed summary types ────────────────────────────────────────────────────

test_that("only type-2 datasets contribute to implied shape check", {
  # Type-1 and type-3 datasets mixed with a low-shape type-2 dataset
  ds <- list(
    d1 = .t1(10, 2, 25),           # type-1 — ignored by heuristic
    d2 = .t3(10, 4),               # type-3 — ignored by heuristic
    d3 = .t2(10, 6, 14)            # type-2 — phi_impl ≈ 2.85
  )
  expect_true(gamma_type2_reliable(ds, verbose = FALSE))
})
