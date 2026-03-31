# Shared dataset fixtures used across multiple test files.
# testthat automatically sources files beginning with "helper" before any tests run.

.ds_range <- function() {
  list(
    d1 = list(median = 5, min = 1, max = 10, n = 30),
    d2 = list(median = 7, min = 2, max = 15, n = 25)
  )
}

.ds_iqr <- function() {
  list(
    d1 = list(median = 5, Q1 = 3, Q3 = 8, n = 40),
    d2 = list(median = 6, Q1 = 4, Q3 = 9, n = 35)
  )
}

.ds_meansd <- function() {
  list(
    d1 = list(mean = 5.0, sd = 2.0, n = 50),
    d2 = list(mean = 7.0, sd = 3.0, n = 40)
  )
}

# summary type 4 — exact frequency table (freq_value / freq_count)
.ds_freq4 <- function() {
  list(
    d1 = list(freq_value = c(1, 2, 3, 4, 5),
              freq_count = c(2L, 5L, 8L, 4L, 1L)),  # n implicit = 20
    d2 = list(freq_value = c(3, 4, 5, 6),
              freq_count = c(3L, 6L, 5L, 2L), n = 20L)
  )
}

# summary type 5 — interval-censored frequency table
.ds_freq5 <- function() {
  list(
    d1 = list(freq_lower = c(0, 3, 6, 9),
              freq_upper = c(3, 6, 9, 12),
              freq_count = c(5L, 10L, 8L, 3L)),   # n implicit = 26
    d2 = list(freq_lower = c(1, 4, 7),
              freq_upper = c(4, 7, 10),
              freq_count = c(4L, 9L, 3L), n = 16L)
  )
}

# Build a Stan data list from datasets, suppressing the n < 5 warning
.make_sd <- function(datasets, dist_type = 1) {
  suppressWarnings(prepare_stan_data_from_datasets(datasets, dist_type = dist_type))
}

# Build five mean+sd datasets from parallel vectors — avoids n < 5 warning
.five_meansd <- function(mean_vals, sd_vals) {
  mapply(function(m, s) list(mean = m, sd = s, n = 50),
         mean_vals, sd_vals, SIMPLIFY = FALSE)
}
