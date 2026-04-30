# analysis/run_meta_analysis.R
# =============================================================================
# Classical meta-analysis of incubation period central estimates
# using metamean() from the {meta} package.
#
# For each pathogen with >= 5 usable datasets the script:
#   - Pools the mean via random effects, log-transformed (sm = "MLN")
#   - Approximates mean/SD from median + IQR or median + range (Luo/Shi)
#   - Computes mean/SD from frequency tables and interval-censored data
#   - Includes subgroup analysis when >= 2 distinct subgroups exist
#   - Saves a forest plot to results/meta_analysis/<pathogen>_forest.pdf
#
# Exclusion rule for interval-censored datasets:
#   Skip any entry where ALL lower bounds equal 0.1 (upper-bound-only data
#   carrying no lower-bound information; e.g. COVID-19 d40).
# =============================================================================

suppressPackageStartupMessages({
  library(meta)
  library(here)
  library(rstan)
})
devtools::load_all(here(), quiet = TRUE)

MIN_STUDIES <- 5L
OUTDIR      <- here("results", "meta_analysis")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# Helper functions
# =============================================================================

# TRUE when all freq_lower values are 0.1 (upper-bound-only interval censoring)
.is_upper_only <- function(d) {
  !is.null(d$freq_lower) && all(abs(d$freq_lower - 0.1) < 1e-9)
}

# "Surname (year)" from source string; falls back to a truncated source
.short_cite <- function(src) {
  if (is.null(src) || !nzchar(src)) return("Unknown")
  m <- regmatches(src, regexpr("^[^(]+\\(\\d{4}\\)", src))
  if (length(m) == 1L && nzchar(m)) return(trimws(m))
  if (nchar(src) > 35L) paste0(substr(src, 1L, 32L), "...") else src
}

# Extract mean, SD, n from a single dataset entry.
# For interval-censored data the SD combines between-midpoint variance and
# within-bin variance (uniform assumption): Var_total = Var_between + E[Var_within].
.extract_central <- function(d) {
  n <- if (!is.null(d$n)) d$n else sum(d$freq_count)

  if (!is.null(d$mean)) {
    list(mean = d$mean, sd = d$sd, n = n,
         median = NA_real_, q1 = NA_real_, q3 = NA_real_,
         min_val = NA_real_, max_val = NA_real_)

  } else if (!is.null(d$median) && !is.null(d$Q1)) {
    list(mean = NA_real_, sd = NA_real_, n = n,
         median = d$median, q1 = d$Q1, q3 = d$Q3,
         min_val = NA_real_, max_val = NA_real_)

  } else if (!is.null(d$median)) {
    list(mean = NA_real_, sd = NA_real_, n = n,
         median = d$median, q1 = NA_real_, q3 = NA_real_,
         min_val = d$min, max_val = d$max)

  } else if (!is.null(d$freq_value)) {
    vals <- rep(d$freq_value, d$freq_count)
    list(mean    = mean(vals),
         sd      = if (length(vals) > 1L) sd(vals) else NA_real_,
         n       = n,
         median  = NA_real_, q1 = NA_real_, q3 = NA_real_,
         min_val = NA_real_, max_val = NA_real_)

  } else if (!is.null(d$freq_lower)) {
    mids        <- (d$freq_lower + d$freq_upper) / 2
    vals        <- rep(mids, d$freq_count)
    between_var <- if (length(vals) > 1L) var(vals) else 0
    within_var  <- mean(rep((d$freq_upper - d$freq_lower)^2 / 12, d$freq_count))
    list(mean    = mean(vals),
         sd      = sqrt(between_var + within_var),
         n       = n,
         median  = NA_real_, q1 = NA_real_, q3 = NA_real_,
         min_val = NA_real_, max_val = NA_real_)

  } else {
    NULL
  }
}

# Build the data-frame fed to metamean() from a pathogen dataset list.
.build_meta_df <- function(datasets) {
  nms  <- names(datasets)
  rows <- lapply(seq_along(datasets), function(i) {
    d <- datasets[[i]]
    if (.is_upper_only(d)) return(NULL)
    ex <- .extract_central(d)
    if (is.null(ex)) return(NULL)
    data.frame(
      study    = .short_cite(d$source),
      n        = ex$n,
      mean     = ex$mean,
      sd       = ex$sd,
      median   = ex$median,
      q1       = ex$q1,
      q3       = ex$q3,
      min_val  = ex$min_val,
      max_val  = ex$max_val,
      subgroup = if (!is.null(d$subgroup)) d$subgroup else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df
}

# Fit metamean() with random effects and log transformation.
# Passes subgroup only when >= 2 distinct non-NA subgroups are present.
.run_metamean <- function(df) {
  sgs          <- unique(na.omit(df$subgroup))
  use_subgroup <- length(sgs) >= 2L

  # Studies with NA subgroup get labelled "Unclassified" when subgroup
  # analysis is active, so they appear in the forest plot correctly.
  sg_vec <- NULL
  if (use_subgroup) {
    sg_vec <- ifelse(is.na(df$subgroup), "Unclassified", df$subgroup)
  }

  metamean(
    n        = df$n,
    mean     = df$mean,
    sd       = df$sd,
    median   = df$median,
    q1       = df$q1,
    q3       = df$q3,
    min      = df$min_val,
    max      = df$max_val,
    studlab  = df$study,
    subgroup = sg_vec,
    method.mean = "Luo",
    method.sd   = "Shi",
    sm          = "MLN",
    common      = FALSE,
    random      = TRUE
  )
}

# Extract tau posterior summary from the Stan lognormal fit.
# Prefers the "all" slot (matching the unfiltered meta-analysis);
# falls back to "filtered" if "all" is unavailable.
.get_stan_tau <- function(main_res, pathogen_name) {
  pg <- main_res[[pathogen_name]]
  if (is.null(pg)) return(NULL)

  slot_nm <- if (!is.null(pg[["all"]][["lognormal"]])) "all"
             else if (!is.null(pg[["filtered"]][["lognormal"]])) "filtered"
             else return(NULL)

  sfit <- pg[[slot_nm]][["lognormal"]][["fit"]]
  if (is.null(sfit)) return(NULL)

  s <- tryCatch(
    rstan::summary(sfit, pars = "tau")$summary,
    error = function(e) NULL
  )
  if (is.null(s)) return(NULL)

  list(
    median = s[1, "50%"],
    lo     = s[1, "2.5%"],
    hi     = s[1, "97.5%"],
    slot   = slot_nm
  )
}

# Save a forest plot for one metamean result.
.save_forest <- function(m, df, label, path) {
  n_studies  <- nrow(df)
  n_subgroups <- if (!is.null(m$subgroup))
    length(unique(na.omit(m$subgroup))) else 0L
  # Allow extra lines for subgroup headers + within-subgroup pooled rows
  n_rows  <- n_studies + n_subgroups * 3L + 5L
  fig_h   <- max(7, 0.22 * n_rows + 3)

  pdf(path, width = 16, height = fig_h)
  on.exit(dev.off(), add = TRUE)

  forest(
    m,
    sortvar             = TE,
    print.tau2          = TRUE,
    print.I2            = TRUE,
    print.Q             = FALSE,
    digits              = 2L,
    xlab                = "Mean incubation period (days)",
    smlab               = "",
    col.diamond.random  = "#2166AC",
    col.square          = "grey40",
    leftcols            = c("studlab", "n"),
    leftlabs            = c("Study [dataset]", "N"),
    rightcols           = c("effect", "ci"),
    rightlabs           = c("Mean", "95% CI"),
    print.subgroup.name = FALSE,
    header.line         = "both",
    spacing             = 0.8,
    main                = paste0(label, ": random-effects meta-analysis (MLN)")
  )
}


# =============================================================================
# Pathogen registry
# =============================================================================
# Each entry: name (used for file names), data (dataset list), label (for plot)

registry <- list(
  list(name = "Nipah",    data = datasets_Nipah,    label = "Nipah virus disease"),
  list(name = "EVD",      data = datasets_EVD,      label = "Ebola virus disease (EVD)"),
  list(name = "SARS",     data = datasets_SARS,     label = "SARS (SARS-CoV-1)"),
  list(name = "MERS",     data = datasets_MERS,     label = "MERS (MERS-CoV)"),
  list(name = "Measles",  data = datasets_Measles,  label = "Measles"),
  list(name = "Mpox",     data = datasets_Mpox,     label = "Mpox"),
  list(name = "Cholera",  data = datasets_Cholera,  label = "Cholera"),
  list(name = "CCHF",     data = datasets_CCHF,     label = "Crimean-Congo haemorrhagic fever (CCHF)"),
  list(name = "COVID_19", data = datasets_COVID_19, label = "COVID-19"),
  list(name = "Dengue",   data = datasets_Dengue,   label = "Dengue"),
  list(name = "Flu",      data = datasets_flu,      label = "Influenza"),
  list(name = "Typhoid",  data = datasets_typhoid,  label = "Typhoid fever")
)


# =============================================================================
# Run meta-analyses
# =============================================================================

results <- list()

for (pg in registry) {
  df <- .build_meta_df(pg$data)

  n_usable <- if (is.null(df)) 0L else nrow(df)
  if (n_usable < MIN_STUDIES) {
    message("Skipping ", pg$name, ": only ", n_usable, " usable datasets (< ", MIN_STUDIES, ")")
    next
  }

  message("Running: ", pg$name, " (", n_usable, " datasets)")

  m <- tryCatch(
    withCallingHandlers(
      .run_metamean(df),
      warning = function(w) {
        message("  WARN:  ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )

  if (is.null(m)) next

  results[[pg$name]] <- list(meta = m, df = df, label = pg$label)
  message("  Done. Pooled mean = ",
          round(exp(m$TE.random), 2), " days  [",
          round(exp(m$lower.random), 2), ", ",
          round(exp(m$upper.random), 2), "]  I2 = ",
          round(m$I2 * 100, 1), "%")
}

message("\n", length(results), " pathogens analysed.")


# =============================================================================
# Generate and save forest plots
# =============================================================================

# Load Stan results to annotate forest plots with Bayesian tau estimates.
stan_res_path <- here("results", "main_results.rds")
stan_res      <- if (file.exists(stan_res_path)) {
  message("Loading Stan results from ", stan_res_path)
  readRDS(stan_res_path)
} else {
  message("main_results.rds not found — forest plots will not include Stan tau")
  NULL
}

for (nm in names(results)) {
  res  <- results[[nm]]
  path <- file.path(OUTDIR, paste0(nm, "_forest.pdf"))

  tryCatch(
    .save_forest(res$meta, res$df, res$label, path),
    error = function(e) message("  Forest plot failed for ", nm, ": ", conditionMessage(e))
  )

  message("Saved: ", path)
}


# =============================================================================
# Save pooled summary table + tau comparison
# =============================================================================

summary_rows <- lapply(names(results), function(nm) {
  res <- results[[nm]]
  m   <- res$meta

  stan_tau <- if (!is.null(stan_res)) .get_stan_tau(stan_res, nm) else NULL

  data.frame(
    pathogen            = nm,
    label               = res$label,
    k                   = m$k,
    pooled_mean         = round(exp(m$TE.random),    2),
    ci_lower            = round(exp(m$lower.random), 2),
    ci_upper            = round(exp(m$upper.random), 2),
    I2_pct              = round(m$I2 * 100,          1),
    metamean_tau        = round(sqrt(m$tau2),         3),
    stan_tau_median     = if (!is.null(stan_tau)) round(stan_tau$median, 3) else NA_real_,
    stan_tau_lo         = if (!is.null(stan_tau)) round(stan_tau$lo,     3) else NA_real_,
    stan_tau_hi         = if (!is.null(stan_tau)) round(stan_tau$hi,     3) else NA_real_,
    stan_slot           = if (!is.null(stan_tau)) stan_tau$slot             else NA_character_,
    stringsAsFactors    = FALSE
  )
})
summary_tbl <- do.call(rbind, summary_rows)

write.csv(summary_tbl,
          file.path(OUTDIR, "pooled_estimates.csv"),
          row.names = FALSE, quote = FALSE)

saveRDS(results, file.path(OUTDIR, "meta_analysis_results.rds"))

message("\nSummary table (metamean τ vs Stan log-normal τ):\n")
print(summary_tbl[, c("pathogen", "k", "I2_pct", "metamean_tau",
                       "stan_tau_median", "stan_tau_lo", "stan_tau_hi")],
      row.names = FALSE)
message("\nAll outputs written to: ", OUTDIR)

# =============================================================================
# Explanatory note (printed to console for reference)
# =============================================================================
message("
── Interpreting tau comparisons ──────────────────────────────────────────────
Both tau estimates are on the log scale and measure between-study variability
in the incubation period location parameter:

  metamean tau: SD of log(mean_i) across studies, estimated by REML under
    a random-effects model using the inverse-variance method.

  Stan log-normal tau: posterior SD of loc_d_i (= log median_i for log-normal)
    across studies, estimated from the joint likelihood of all data types.

For a log-normal distribution: log(mean_i) = loc_d_i + 0.5 * sigma_i^2,
where sigma = 1/phi (within-study log-scale SD). Because the Stan model uses
a shared phi across studies, the two tau estimates differ only by a constant
within-study variance term and are directly comparable.

Discrepancies can arise from:
  (a) Different data used — metamean includes approximated means from
      median/range/IQR entries; Stan uses the full likelihood.
  (b) Filtering — Stan results shown are from the 'all' slot unless
      only 'filtered' was available.
  (c) Prior regularisation — Stan's half-normal prior on tau pulls estimates
      towards zero relative to the REML estimator.
──────────────────────────────────────────────────────────────────────────────
")
