# Tests for filter_datasets() in R/utils.R

test_that("NULL args returns the input list unchanged", {
  ds <- list(
    a = list(mean = 5, sd = 2, n = 10, subgroup = "adults"),
    b = list(mean = 6, sd = 3, n = 12)
  )
  expect_identical(filter_datasets(ds), ds)
  expect_identical(filter_datasets(ds, subgroup = NULL, location = NULL), ds)
})

test_that("drops entries that lack the queried field entirely", {
  ds <- list(
    a = list(mean = 5, sd = 2, n = 10, subgroup = "adults"),
    b = list(mean = 6, sd = 3, n = 12)               # no subgroup field
  )
  res <- filter_datasets(ds, subgroup = "adults")
  expect_equal(length(res), 1L)
  expect_named(res, "a")
})

test_that("keeps only entries whose subgroup matches", {
  ds <- list(
    a = list(mean = 5, sd = 1, n = 10, subgroup = "adults"),
    b = list(mean = 6, sd = 2, n = 12, subgroup = "children"),
    c = list(mean = 7, sd = 3, n = 15, subgroup = "adults")
  )
  res <- filter_datasets(ds, subgroup = "adults")
  expect_equal(length(res), 2L)
  expect_true(all(names(res) %in% c("a", "c")))
})

test_that("accepts a vector of subgroup values", {
  ds <- list(
    a = list(mean = 5, sd = 1, n = 10, subgroup = "adults"),
    b = list(mean = 6, sd = 2, n = 12, subgroup = "children"),
    c = list(mean = 7, sd = 3, n = 15, subgroup = "elderly")
  )
  res <- filter_datasets(ds, subgroup = c("adults", "elderly"))
  expect_equal(length(res), 2L)
  expect_true(all(names(res) %in% c("a", "c")))
})

test_that("filters by location", {
  ds <- list(
    a = list(mean = 5, sd = 1, n = 10, location = "Turkey"),
    b = list(mean = 6, sd = 2, n = 12, location = "Iran"),
    c = list(mean = 7, sd = 3, n = 15, location = "Turkey")
  )
  res <- filter_datasets(ds, location = "Turkey")
  expect_equal(length(res), 2L)
  expect_true(all(names(res) %in% c("a", "c")))
})

test_that("applies subgroup AND location as an intersection", {
  ds <- list(
    a = list(mean = 5, sd = 1, n = 10, subgroup = "adults",   location = "Turkey"),
    b = list(mean = 6, sd = 2, n = 12, subgroup = "children", location = "Turkey"),
    c = list(mean = 7, sd = 3, n = 15, subgroup = "adults",   location = "Iran")
  )
  res <- filter_datasets(ds, subgroup = "adults", location = "Turkey")
  expect_equal(length(res), 1L)
  expect_named(res, "a")
})

test_that("returns an empty list when nothing matches", {
  ds <- list(a = list(mean = 5, sd = 1, n = 10, subgroup = "adults"))
  res <- filter_datasets(ds, subgroup = "children")
  expect_equal(length(res), 0L)
})

test_that("preserves the names of retained entries", {
  ds <- list(
    foo = list(mean = 5, sd = 1, n = 10, subgroup = "adults"),
    bar = list(mean = 6, sd = 2, n = 12, subgroup = "children")
  )
  res <- filter_datasets(ds, subgroup = "adults")
  expect_named(res, "foo")
})
