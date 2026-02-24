test_that("check_positive rejects non-positive values", {
  expect_error(delaydistribution:::check_positive(-1),  "strictly positive")
  expect_error(delaydistribution:::check_positive(0),   "strictly positive")
  expect_error(delaydistribution:::check_positive("a"), "strictly positive")
  expect_silent(delaydistribution:::check_positive(1))
  expect_silent(delaydistribution:::check_positive(c(0.1, 2, 10)))
})

test_that("check_scalar rejects non-scalars", {
  expect_error(delaydistribution:::check_scalar(c(1, 2)), "single finite")
  expect_error(delaydistribution:::check_scalar(Inf),     "single finite")
  expect_error(delaydistribution:::check_scalar("a"),     "single finite")
  expect_silent(delaydistribution:::check_scalar(3.14))
})
