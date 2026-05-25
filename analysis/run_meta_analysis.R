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
    total_n  <- sum(d$freq_count)
    mean_val <- sum(d$freq_value * d$freq_count) / total_n
    sd_val   <- if (total_n > 1L) {
      sqrt(sum(d$freq_count * (d$freq_value - mean_val)^2) / (total_n - 1L))
    } else {
      NA_real_
    }
    list(mean    = mean_val,
         sd      = sd_val,
         n       = n,
         median  = NA_real_, q1 = NA_real_, q3 = NA_real_,
         min_val = NA_real_, max_val = NA_real_)

  } else if (!is.null(d$freq_lower)) {
    mids        <- (d$freq_lower + d$freq_upper) / 2
    total_n     <- sum(d$freq_count)
    mean_val    <- sum(mids * d$freq_count) / total_n
    between_var <- if (total_n > 1L) {
      sum(d$freq_count * (mids - mean_val)^2) / (total_n - 1L)
    } else {
      0
    }
    within_var  <- sum(d$freq_count * (d$freq_upper - d$freq_lower)^2 / 12) / total_n
    list(mean    = mean_val,
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

# Extract tau posterior summary from the best-fitting Stan model.
# Best distribution is taken from model_weights (highest stacking weight).
# Prefers the "all" slot; falls back to "filtered" if "all" is unavailable.
# Returns NULL when < 5 datasets were used (tau not reliably estimated).
.get_stan_tau <- function(main_res, model_weights, pathogen_name) {
  pg <- main_res[[pathogen_name]]
  if (is.null(pg)) return(NULL)

  # Identify best-fitting distribution for this pathogen
  wts      <- model_weights[[pathogen_name]]
  best_dist <- if (!is.null(wts) && length(wts) > 0L) names(which.max(wts)) else "lognormal"

  slot_nm <- if (!is.null(pg[["all"]][[best_dist]])) "all"
             else if (!is.null(pg[["filtered"]][[best_dist]])) "filtered"
             else return(NULL)

  fit_obj <- pg[[slot_nm]][[best_dist]]
  if (length(fit_obj$datasets) < 5L) return(NULL)
  sfit <- fit_obj$fit
  if (is.null(sfit)) return(NULL)

  s <- tryCatch(
    rstan::summary(sfit, pars = "tau")$summary,
    error = function(e) NULL
  )
  if (is.null(s)) return(NULL)

  list(
    median    = s[1, "50%"],
    lo        = s[1, "2.5%"],
    hi        = s[1, "97.5%"],
    slot      = slot_nm,
    best_dist = best_dist
  )
}

# Save a forest plot for one metamean result.
# meta 8.x uses grid graphics and auto-calculates the correct figure height
# when file= is supplied, so we let it open and close the PDF device itself.
.save_forest <- function(m, df, label, path) {
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
    print.subgroup.name           = FALSE,
    header.line                   = "both",
    addrows.below.overall         = 2L,
    fontsize                      = 11,
    colgap.forest.left            = "4cm",
    file                          = path,
    width                         = 16
  )
}


# =============================================================================
# Pathogen registry
# =============================================================================
# Each entry: name (used for file names), data (dataset list), label (for plot)

registry <- list(
  list(name = "Nipah",    data = datasets_Nipah,          label = "Nipah virus disease"),
  list(name = "EVD",      data = datasets_EVD,            label = "Ebola virus disease (EVD)"),
  list(name = "SARS",     data = datasets_SARS,           label = "SARS (SARS-CoV-1)"),
  list(name = "MERS",     data = datasets_MERS,           label = "MERS (MERS-CoV)"),
  list(name = "Measles",  data = datasets_Measles,        label = "Measles"),
  list(name = "Mpox",     data = datasets_Mpox,           label = "Mpox"),
  list(name = "Cholera",  data = datasets_Cholera,        label = "Cholera"),
  list(name = "CCHF",     data = datasets_CCHF_extended,  label = "Crimean-Congo haemorrhagic fever (CCHF)"),
  list(name = "COVID_19", data = datasets_COVID_19,       label = "COVID-19"),
  list(name = "Dengue",   data = datasets_Dengue,         label = "Dengue"),
  list(name = "Flu",      data = datasets_flu,            label = "Influenza"),
  list(name = "Typhoid",  data = datasets_typhoid,        label = "Typhoid fever")
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

# Load Stan results and model weights for the tau comparison table.
stan_res_path <- here("results", "main_results.rds")
weights_path  <- here("results", "model_weights.rds")

stan_res <- if (file.exists(stan_res_path)) {
  message("Loading Stan results from ", stan_res_path)
  readRDS(stan_res_path)
} else {
  message("main_results.rds not found — Stan tau will not be reported")
  NULL
}

model_weights <- if (file.exists(weights_path)) {
  message("Loading model weights from ", weights_path)
  readRDS(weights_path)
} else {
  message("model_weights.rds not found — defaulting to lognormal for Stan tau")
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

  stan_tau <- if (!is.null(stan_res)) .get_stan_tau(stan_res, model_weights, nm) else NULL

  data.frame(
    pathogen            = nm,
    label               = res$label,
    k                   = m$k,
    pooled_mean         = round(exp(m$TE.random),    2),
    ci_lower            = round(exp(m$lower.random), 2),
    ci_upper            = round(exp(m$upper.random), 2),
    I2_pct              = round(m$I2 * 100,          1),
    metamean_tau        = round(sqrt(m$tau2),         3),
    stan_tau_median     = if (!is.null(stan_tau)) round(stan_tau$median,    3) else NA_real_,
    stan_tau_lo         = if (!is.null(stan_tau)) round(stan_tau$lo,        3) else NA_real_,
    stan_tau_hi         = if (!is.null(stan_tau)) round(stan_tau$hi,        3) else NA_real_,
    stan_best_dist      = if (!is.null(stan_tau)) stan_tau$best_dist         else NA_character_,
    stan_slot           = if (!is.null(stan_tau)) stan_tau$slot              else NA_character_,
    stringsAsFactors    = FALSE
  )
})
summary_tbl <- do.call(rbind, summary_rows)

write.csv(summary_tbl,
          file.path(OUTDIR, "pooled_estimates.csv"),
          row.names = FALSE, quote = FALSE)

saveRDS(results, file.path(OUTDIR, "meta_analysis_results.rds"))

message("\nSummary table (metamean τ vs Stan τ from best-fitting distribution):\n")
print(summary_tbl[, c("pathogen", "k", "I2_pct", "metamean_tau",
                       "stan_best_dist", "stan_tau_median", "stan_tau_lo", "stan_tau_hi")],
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


# =============================================================================
# LaTeX table
# =============================================================================

generate_meta_analysis_table <- function(
    summary_tbl,
    caption = paste0(
      "Classical meta-analysis of incubation period central estimates ",
      "(random-effects model, log-transformed mean, \\texttt{sm = MLN}). ",
      "$k$ = number of studies. ",
      "Pooled mean and 95\\% confidence interval (CI) are back-transformed to days. ",
      "$I^2$ measures the percentage of total variance due to between-study heterogeneity. ",
      "$\\tau_{\\text{meta}}$ is the between-study SD on the log scale (REML). ",
      "$\\tau_{\\text{Stan}}$ (median and 95\\% credible interval) is from the ",
      "best-fitting Bayesian model (``all data'' slot), also on the log scale."
    ),
    label = "tab:meta_analysis"
) {

  dist_labels <- c(
    lognormal = "Log-normal",
    gamma     = "Gamma",
    weibull   = "Weibull",
    burr      = "Burr XII",
    gengamma  = "Gen. Gamma"
  )

  pathogen_labels <- ddsynth:::.PATHOGEN_LABELS
  group_map       <- ddsynth:::.RESULT_KEY_TO_GROUP
  group_order     <- ddsynth:::.GROUP_ORDER

  # Format a number: fixed decimal places, dash for NA
  .fmt <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, digits = d, format = "f"))

  # ── Build body rows grouped by epidemiological category ───────────────────
  body_lines <- character(0L)
  first_group <- TRUE

  for (grp in group_order) {
    keys_in_grp <- names(group_map)[group_map == grp]
    rows_in_grp <- summary_tbl[summary_tbl$pathogen %in% keys_in_grp, , drop = FALSE]
    if (nrow(rows_in_grp) == 0L) next

    # Group header row
    if (!first_group) body_lines <- c(body_lines, "\\midrule")
    body_lines <- c(body_lines,
      paste0("\\rowcolor{gray!15}",
             "\\multicolumn{7}{l}{\\textit{", grp, "}} \\\\")
    )
    first_group <- FALSE

    # Data rows (order matches group_order within group)
    for (i in seq_len(nrow(rows_in_grp))) {
      r        <- rows_in_grp[i, ]
      plabel   <- pathogen_labels[r$pathogen]
      if (is.na(plabel)) plabel <- r$pathogen

      dist_lab <- dist_labels[r$stan_best_dist]
      if (is.na(dist_lab)) dist_lab <- r$stan_best_dist

      mean_ci  <- paste0(.fmt(r$pooled_mean), " [", .fmt(r$ci_lower), ", ", .fmt(r$ci_upper), "]")
      stan_ci  <- if (!is.na(r$stan_tau_median))
        paste0(.fmt(r$stan_tau_median, 3), " [", .fmt(r$stan_tau_lo, 3), ", ",
               .fmt(r$stan_tau_hi, 3), "]")
      else "--"

      body_lines <- c(body_lines, paste0(
        "\\quad ", plabel, " & ",
        r$k, " & ",
        mean_ci, " & ",
        .fmt(r$I2_pct, 1), "\\% & ",
        .fmt(r$metamean_tau, 3), " & ",
        stan_ci, " & ",
        dist_lab, " \\\\"
      ))
    }
  }

  # ── Column header ──────────────────────────────────────────────────────────
  col_header <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    "Pathogen & $k$ & ",
    "\\multicolumn{1}{c}{Pooled mean [95\\% CI] (days)} & ",
    "\\multicolumn{1}{c}{$I^2$} & ",
    "\\multicolumn{1}{c}{$\\tau_{\\text{meta}}$} & ",
    "\\multicolumn{1}{c}{$\\tau_{\\text{Stan}}$ [95\\% CrI]} & ",
    "\\multicolumn{1}{c}{Best model} \\\\"
  )

  # ── Assemble full table ────────────────────────────────────────────────────
  lines <- c(
    "% Requires \\usepackage{booktabs}, \\usepackage{longtable},",
    "% \\usepackage[table]{xcolor} in preamble.",
    "\\begin{longtable}{@{} l r l r r l l @{}}",
    paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
    "\\toprule",
    col_header,
    "\\midrule",
    "\\endfirsthead",
    "%",
    "\\multicolumn{7}{l}{\\small\\textit{continued from previous page}} \\\\",
    "\\toprule",
    col_header,
    "\\midrule",
    "\\endhead",
    "%",
    "\\multicolumn{7}{r}{\\small\\textit{continued on next page}} \\\\",
    "\\endfoot",
    "%",
    "\\bottomrule",
    "\\endlastfoot",
    "%",
    body_lines,
    "\\end{longtable}"
  )

  paste(lines, collapse = "\n")
}

tex <- generate_meta_analysis_table(summary_tbl)
tex_path <- file.path(OUTDIR, "meta_analysis_table.tex")
writeLines(tex, tex_path)
message("LaTeX table written to: ", tex_path)
cat(tex, "\n")


# =============================================================================
# Dataset composition table
# =============================================================================

# Classify a single dataset as "indiv" (D/E) or "sumstat" (A/B/C).
.composition_type <- function(d) {
  if (!is.null(d$freq_lower) && !is.null(d$freq_upper)) return("indiv")
  if (!is.null(d$freq_value) && !is.null(d$freq_count)) return("indiv")
  "sumstat"
}

# Aggregate composition stats from a list of datasets.
# Returns list(k_total, N_total, k_indiv, N_indiv, k_sumstat, N_sumstat).
.composition_stats <- function(datasets) {
  k_total   <- 0L; N_total   <- 0L
  k_indiv   <- 0L; N_indiv   <- 0L
  k_sumstat <- 0L; N_sumstat <- 0L

  for (d in datasets) {
    n  <- if (!is.null(d$n)) d$n else if (!is.null(d$freq_count)) sum(d$freq_count) else NA_integer_
    tp <- .composition_type(d)
    k_total <- k_total + 1L
    N_total <- N_total + if (!is.na(n)) n else 0L
    if (tp == "indiv") {
      k_indiv <- k_indiv + 1L; N_indiv <- N_indiv + if (!is.na(n)) n else 0L
    } else {
      k_sumstat <- k_sumstat + 1L; N_sumstat <- N_sumstat + if (!is.na(n)) n else 0L
    }
  }
  list(k_total=k_total, N_total=N_total, k_indiv=k_indiv,
       N_indiv=N_indiv, k_sumstat=k_sumstat, N_sumstat=N_sumstat)
}

# Format composition stats as a LaTeX data row string (6 cells, no pathogen label).
.fmt_comp_cells <- function(s) {
  paste(
    s$k_total, format(s$N_total, big.mark = ","),
    if (s$k_indiv > 0L) s$k_indiv else "--",
    if (s$k_indiv > 0L) format(s$N_indiv, big.mark = ",") else "--",
    if (s$k_sumstat > 0L) s$k_sumstat else "--",
    if (s$k_sumstat > 0L) format(s$N_sumstat, big.mark = ",") else "--",
    sep = " & "
  )
}

generate_dataset_composition_table <- function(
    full_registry,
    caption = paste0(
      "Composition of incubation period datasets included in the analysis. ",
      "For each pathogen (and subgroup where applicable), $k$ is the number of datasets ",
      "and $N$ is the total number of individual observations. ",
      "\\textit{Individual-level data} denotes datasets provided as frequency tables ",
      "(exact or interval-censored); \\textit{summary statistics only} denotes datasets ",
      "reported as mean $\\pm$ SD, median $\\pm$ IQR, or median $\\pm$ range."
    ),
    label = "tab:dataset_composition"
) {
  pathogen_labels <- ddsynth:::.PATHOGEN_LABELS
  subgroup_labels <- ddsynth:::.SUBGROUP_LABELS
  group_map       <- ddsynth:::.RESULT_KEY_TO_GROUP
  group_order     <- ddsynth:::.GROUP_ORDER

  # Column header (two-row: group spans on top, k/N labels below)
  hdr1 <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    " & \\multicolumn{2}{c}{\\textbf{Total}} &",
    " \\multicolumn{2}{c}{\\textbf{Individual-level}} &",
    " \\multicolumn{2}{c}{\\textbf{Summary statistics only}} \\\\"
  )
  hdr2 <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    "\\textbf{Pathogen} & $k$ & $N$ & $k$ & $N$ & $k$ & $N$ \\\\"
  )

  body_lines  <- character(0L)
  first_group <- TRUE

  # Index registry by pathogen name for quick lookup
  reg_index <- stats::setNames(
    lapply(full_registry, `[[`, "data"),
    sapply(full_registry, `[[`, "name")
  )

  for (grp in group_order) {
    keys_in_grp <- names(group_map)[group_map == grp]
    keys_present <- intersect(keys_in_grp, names(reg_index))
    if (length(keys_present) == 0L) next

    if (!first_group) body_lines <- c(body_lines, "\\midrule")
    body_lines <- c(body_lines,
      paste0("\\rowcolor{gray!15}",
             "\\multicolumn{7}{@{}l}{\\textbf{\\textit{", grp, "}}} \\\\")
    )
    first_group <- FALSE

    first_in_group <- TRUE
    for (nm in keys_in_grp) {
      if (!nm %in% names(reg_index)) next
      datasets <- reg_index[[nm]]
      plabel   <- pathogen_labels[[nm]]
      if (is.null(plabel) || is.na(plabel)) plabel <- nm

      # Overall totals for this pathogen
      stats_all <- .composition_stats(datasets)

      # Detect subgroups (NULL entries → NA, then omit)
      sgs <- unique(na.omit(sapply(datasets, function(d) {
        sg <- d$subgroup; if (is.null(sg)) NA_character_ else sg
      })))
      unclass_datasets <- Filter(function(d) is.null(d$subgroup), datasets)

      # Show subgroup rows only when ≥2 distinct subgroups OR ≥1 subgroup + unclassified.
      # Skip when all datasets belong to a single subgroup with no unclassified entries
      # (the parent row already captures everything).
      has_subgroups <- (length(sgs) >= 2L) ||
                       (length(sgs) == 1L && length(unclass_datasets) > 0L)

      if (!first_in_group)
        body_lines <- c(body_lines, "\\addlinespace[4pt]")

      body_lines <- c(body_lines, paste0(
        "\\textbf{", plabel, "} & ",
        .fmt_comp_cells(stats_all), " \\\\"
      ))
      first_in_group <- FALSE

      if (has_subgroups) {
        for (sg in sgs) {
          sg_datasets <- Filter(function(d) {
            dsg <- d$subgroup; !is.null(dsg) && identical(dsg, sg)
          }, datasets)
          if (length(sg_datasets) == 0L) next
          sg_stats <- .composition_stats(sg_datasets)
          sg_lab <- if (!is.na(subgroup_labels[sg])) subgroup_labels[[sg]] else sg
          body_lines <- c(body_lines, paste0(
            "\\quad \\textit{", sg_lab, "} & ",
            .fmt_comp_cells(sg_stats), " \\\\"
          ))
        }
        if (length(unclass_datasets) > 0L) {
          uc_stats <- .composition_stats(unclass_datasets)
          body_lines <- c(body_lines, paste0(
            "\\quad \\textit{Unclassified} & ",
            .fmt_comp_cells(uc_stats), " \\\\"
          ))
        }
      }
    }
  }

  # Assemble full table
  lines <- c(
    "% Requires \\usepackage{booktabs}, \\usepackage{longtable},",
    "% \\usepackage[table]{xcolor} in preamble.",
    "\\begin{longtable}{@{} l r r r r r r @{}}",
    paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
    "\\toprule",
    hdr1,
    hdr2,
    "\\midrule",
    "\\endfirsthead",
    "%",
    "\\multicolumn{7}{l}{\\small\\textit{continued from previous page}} \\\\",
    "\\toprule",
    hdr1,
    hdr2,
    "\\midrule",
    "\\endhead",
    "%",
    "\\multicolumn{7}{r}{\\small\\textit{continued on next page}} \\\\",
    "\\endfoot",
    "%",
    "\\bottomrule",
    "\\endlastfoot",
    "%",
    body_lines,
    "\\end{longtable}"
  )

  paste(lines, collapse = "\n")
}

# Full registry covering all pathogens (for composition table)
full_registry <- list(
  list(name = "EVD",      data = datasets_EVD),
  list(name = "MVD",      data = datasets_MVD),
  list(name = "Lassa",    data = datasets_Lassa),
  list(name = "CCHF",     data = datasets_CCHF_extended),
  list(name = "COVID_19", data = datasets_COVID_19),
  list(name = "SARS",     data = datasets_SARS),
  list(name = "MERS",     data = datasets_MERS),
  list(name = "Flu",      data = datasets_flu),
  list(name = "Dengue",   data = datasets_Dengue),
  list(name = "Zika",     data = datasets_Zika),
  list(name = "RVF",      data = datasets_RVF),
  list(name = "YFV",      data = datasets_YFV),
  list(name = "Nipah",    data = datasets_Nipah),
  list(name = "Mpox",     data = datasets_Mpox),
  list(name = "Measles",  data = datasets_Measles),
  list(name = "Smallpox", data = datasets_Smallpox),
  list(name = "Cholera",  data = datasets_Cholera),
  list(name = "Typhoid",  data = datasets_typhoid)
)

comp_tex      <- generate_dataset_composition_table(full_registry)
comp_tex_path <- file.path(OUTDIR, "dataset_composition_table.tex")
writeLines(comp_tex, comp_tex_path)
message("Dataset composition table written to: ", comp_tex_path)
cat(comp_tex, "\n")
