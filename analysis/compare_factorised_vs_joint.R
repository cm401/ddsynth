# =============================================================================
# compare_factorised_vs_joint.R
# -----------------------------------------------------------------------------
# Fits the joint-likelihood Stan model to the same real data used in
# main_analysis.R and produces a comparison table of posterior summaries
# under the two likelihoods.
#
# The likelihood differs only for summary_type 1 (median + range) and
# summary_type 2 (median + IQR), where the factorised model treats each
# reported order statistic independently while the joint model uses the
# true joint density of the three order statistics.  Pathogens with no
# type-1/2 datasets are included in the table with a note that the two
# models are equivalent.
#
# Strategy
# --------
# Rather than re-preparing data, the script reuses the stan_data objects
# stored in main_results.rds, ensuring an exact apples-to-apples comparison
# (same data, same priors, same MCMC settings).  The "filtered" analysis is
# used throughout, matching the primary results in the paper.
#
# Output
# ------
#   results/comparison_factorised_vs_joint.rds  — raw stanfit objects + summaries
#   results/comparison_factorised_vs_joint.csv  — formatted comparison table
#
# The RDS contains:
#   $joint_results   nested list: [[pathogen]][[dist]] with $fit, $stan_data
#   $comparison_tbl  data.frame with all summary columns for both models
# =============================================================================

devtools::load_all(here::here(), quiet = TRUE)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# ── 1. Settings ───────────────────────────────────────────────────────────────

CHAINS  <- 4
ITER    <- 12000
WARMUP  <- 2000
THIN    <- 1
SEED    <- 123
CONTROL <- list(adapt_delta = 0.999, max_treedepth = 12, stepsize = 0.01)

# Analysis label to compare (must exist in main_results.rds).
ANALYSIS_LABEL <- "filtered"

# Distributions to include.  Set to NULL to use every distribution present
# in the existing results for each pathogen.
DIST_NAMES <- c("lognormal", "gamma", "weibull", "burr", "gengamma")

# Credible interval width for formatted cells.
CI_PROBS <- c(0.025, 0.5, 0.975)

OUTPUT_DIR  <- here::here("results")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "comparison_factorised_vs_joint.rds")
CSV_FILE    <- file.path(OUTPUT_DIR, "comparison_factorised_vs_joint.csv")

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Wipe joint results and refit from scratch when TRUE.
FORCE_RERUN <- FALSE


# ── 2. Compile Stan models ────────────────────────────────────────────────────

stan_file_fact  <- system.file("stan",
                               "hierarchical_data_synthesis_summary_stats.stan",
                               package = "ddsynth")
stan_file_joint <- system.file("stan",
                               "hierarchical_data_synthesis_summary_stats_joint.stan",
                               package = "ddsynth")

message("Compiling factorised Stan model...")
stan_model_fact  <- rstan::stan_model(stan_file_fact)

message("Compiling joint Stan model...")
stan_model_joint <- rstan::stan_model(stan_file_joint)


# ── 3. Load main results ──────────────────────────────────────────────────────

main_results_file <- file.path(OUTPUT_DIR, "main_results.rds")
if (!file.exists(main_results_file))
  stop("main_results.rds not found in ", OUTPUT_DIR,
       " — run main_analysis.R first.")

message("Loading main results from: ", main_results_file)
main_results <- readRDS(main_results_file)


# ── 4. Load or initialise joint results ──────────────────────────────────────

if (!FORCE_RERUN && file.exists(OUTPUT_FILE)) {
  message("Resuming from existing comparison file: ", OUTPUT_FILE)
  stored <- readRDS(OUTPUT_FILE)
  joint_results <- stored$joint_results
} else {
  joint_results <- list()
}


# ── 5. Helper functions ───────────────────────────────────────────────────────

# Format a posterior sample vector as "median (lo, hi)".
.fmt <- function(x, digits = 2L) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_character_)
  q <- quantile(x, CI_PROBS)
  sprintf("%.*f (%.*f, %.*f)", digits, q[2L], digits, q[1L], digits, q[3L])
}

# Extract comparison row for one stanfit.
# Returns a named list of formatted strings and raw posterior medians.
.extract_summaries <- function(fit, dist_name, n_datasets) {
  sims <- tryCatch(rstan::extract(fit), error = function(e) NULL)
  if (is.null(sims)) return(NULL)

  tau_fmt <- if (n_datasets < 5L)
               "— (n<5)"
             else
               .fmt(sims$tau, digits = 2L)

  kappa_fmt <- if (dist_name %in% c("burr", "gengamma") && !is.null(sims$kappa))
                 .fmt(sims$kappa, digits = 2L)
               else
                 NA_character_

  list(
    pred_median_fmt = .fmt(sims$pred_median, digits = 1L),
    pred_q95_fmt    = .fmt(sims$pred_q95,    digits = 1L),
    mu0_fmt         = .fmt(sims$mu0,         digits = 2L),
    phi_fmt         = .fmt(sims$phi,         digits = 2L),
    tau_fmt         = tau_fmt,
    kappa_fmt       = kappa_fmt,
    # Raw posterior medians for computing deltas
    pred_median_med = median(sims$pred_median[is.finite(sims$pred_median)]),
    pred_q95_med    = median(sims$pred_q95[is.finite(sims$pred_q95)]),
    mu0_med         = median(sims$mu0),
    phi_med         = median(sims$phi),
    tau_med         = if (n_datasets >= 5L) median(sims$tau) else NA_real_
  )
}


# ── 6. Fit joint model for each pathogen × distribution ──────────────────────

pathogens <- names(main_results)

for (pathogen in pathogens) {

  if (is.null(joint_results[[pathogen]]))
    joint_results[[pathogen]] <- list()

  analysis_slot <- main_results[[pathogen]][[ANALYSIS_LABEL]]
  if (is.null(analysis_slot)) {
    message("\n  [SKIP] ", pathogen, ": no '", ANALYSIS_LABEL, "' analysis in main results.")
    next
  }

  dists_to_run <- if (is.null(DIST_NAMES)) names(analysis_slot) else DIST_NAMES

  for (dist_name in dists_to_run) {

    existing_fact <- analysis_slot[[dist_name]]

    # Skip distributions that were skipped or failed in the main analysis.
    if (is.null(existing_fact)) next
    if (isTRUE(existing_fact$skipped)) {
      message("\n  [SKIP] ", pathogen, " / ", dist_name,
              " — skipped in main analysis (", existing_fact$reason, ").")
      next
    }
    if (is.null(existing_fact$fit) || is.null(existing_fact$stan_data)) {
      message("\n  [SKIP] ", pathogen, " / ", dist_name, " — no fit or stan_data.")
      next
    }

    # Skip if already fitted in a previous run.
    if (!is.null(joint_results[[pathogen]][[dist_name]])) {
      message("\n  [SKIP] ", pathogen, " / ", dist_name,
              " — joint fit already in results.")
      next
    }

    stan_data <- existing_fact$stan_data
    n_type12  <- sum(stan_data$summary_type %in% c(1L, 2L))

    # When no type-1/2 datasets are present the two likelihoods are identical;
    # store a sentinel rather than wasting compute.
    if (n_type12 == 0L) {
      message("\n  [NOTE] ", pathogen, " / ", dist_name,
              " — no type-1/2 datasets; likelihoods identical; storing sentinel.")
      joint_results[[pathogen]][[dist_name]] <- list(
        identical_likelihood = TRUE,
        stan_data            = stan_data
      )
      next
    }

    message("\n  Fitting joint model: ", pathogen, " | ", dist_name,
            "  (n_datasets=", stan_data$n_datasets,
            ", n_type12=", n_type12, ")")

    fit_joint <- tryCatch(
      rstan::sampling(
        stan_model_joint,
        data    = stan_data,
        chains  = CHAINS,
        iter    = ITER,
        warmup  = WARMUP,
        thin    = THIN,
        control = CONTROL,
        seed    = SEED,
        refresh = 200
      ),
      error = function(e) {
        message("  [FAIL] Joint sampling failed: ", conditionMessage(e))
        NULL
      }
    )

    joint_results[[pathogen]][[dist_name]] <- list(
      fit                  = fit_joint,
      stan_data            = stan_data,
      identical_likelihood = FALSE
    )

    # Persist after each fit so a crash does not lose earlier work.
    saveRDS(list(joint_results = joint_results), OUTPUT_FILE)
  }
}

message("\n  All joint fits complete.")


# ── 7. Build comparison table ─────────────────────────────────────────────────

dist_labels <- c(
  lognormal = "Log-normal",
  gamma     = "Gamma",
  weibull   = "Weibull",
  burr      = "Burr XII",
  gengamma  = "Gen. gamma"
)

rows <- list()

for (pathogen in pathogens) {

  analysis_slot <- main_results[[pathogen]][[ANALYSIS_LABEL]]
  if (is.null(analysis_slot)) next

  dists_to_run <- if (is.null(DIST_NAMES)) names(analysis_slot) else DIST_NAMES

  for (dist_name in dists_to_run) {

    fact_res  <- analysis_slot[[dist_name]]
    joint_res <- joint_results[[pathogen]][[dist_name]]

    # Skip entirely absent or failed entries.
    if (is.null(fact_res))  next
    if (isTRUE(fact_res$skipped)) next
    if (is.null(fact_res$fit) || is.null(fact_res$stan_data)) next

    stan_data <- fact_res$stan_data
    n_type12  <- sum(stan_data$summary_type %in% c(1L, 2L))

    fact_summ  <- .extract_summaries(fact_res$fit,  dist_name, stan_data$n_datasets)
    if (is.null(fact_summ)) next

    # Joint summaries: either extracted from fit, or copied from factorised
    # when the likelihoods are identical (no type-1/2 data).
    if (!is.null(joint_res) && isTRUE(joint_res$identical_likelihood)) {
      joint_summ <- fact_summ  # identical by construction
      note <- "identical (no type-1/2 data)"
    } else if (!is.null(joint_res) && !is.null(joint_res$fit)) {
      joint_summ <- .extract_summaries(joint_res$fit, dist_name, stan_data$n_datasets)
      note <- ""
    } else {
      joint_summ <- NULL
      note <- "joint fit unavailable"
    }

    row <- data.frame(
      pathogen    = pathogen,
      dist        = dist_labels[[dist_name]],
      n_datasets  = stan_data$n_datasets,
      n_type12    = n_type12,
      note        = note,

      # Factorised model
      fact_pred_median = if (!is.null(fact_summ)) fact_summ$pred_median_fmt else NA,
      fact_pred_q95    = if (!is.null(fact_summ)) fact_summ$pred_q95_fmt    else NA,
      fact_mu0         = if (!is.null(fact_summ)) fact_summ$mu0_fmt         else NA,
      fact_phi         = if (!is.null(fact_summ)) fact_summ$phi_fmt         else NA,
      fact_tau         = if (!is.null(fact_summ)) fact_summ$tau_fmt         else NA,

      # Joint model
      joint_pred_median = if (!is.null(joint_summ)) joint_summ$pred_median_fmt else NA,
      joint_pred_q95    = if (!is.null(joint_summ)) joint_summ$pred_q95_fmt    else NA,
      joint_mu0         = if (!is.null(joint_summ)) joint_summ$mu0_fmt         else NA,
      joint_phi         = if (!is.null(joint_summ)) joint_summ$phi_fmt         else NA,
      joint_tau         = if (!is.null(joint_summ)) joint_summ$tau_fmt         else NA,

      # Absolute difference in posterior median point estimates (joint - factorised)
      delta_pred_median = if (!is.null(fact_summ) && !is.null(joint_summ))
                            round(joint_summ$pred_median_med - fact_summ$pred_median_med, 2L)
                          else NA_real_,
      delta_pred_q95    = if (!is.null(fact_summ) && !is.null(joint_summ))
                            round(joint_summ$pred_q95_med - fact_summ$pred_q95_med, 2L)
                          else NA_real_,

      stringsAsFactors = FALSE
    )

    rows[[length(rows) + 1L]] <- row
  }
}

comparison_tbl <- do.call(rbind, rows)
rownames(comparison_tbl) <- NULL


# ── 8. Save and print ─────────────────────────────────="──────────────────────

saveRDS(
  list(joint_results = joint_results, comparison_tbl = comparison_tbl),
  OUTPUT_FILE
)
message("\nRDS saved to: ", OUTPUT_FILE)

write.csv(comparison_tbl, CSV_FILE, row.names = FALSE)
message("CSV saved to: ", CSV_FILE)


# ── 9. Print formatted table ──────────────────────────────────────────────────

message("\n", strrep("=", 80))
message("COMPARISON: FACTORISED vs JOINT LIKELIHOOD  (analysis = '",
        ANALYSIS_LABEL, "')")
message(strrep("=", 80))
message(
  "\nColumns:\n",
  "  pathogen          : pathogen name\n",
  "  dist              : distribution\n",
  "  n_datasets        : total number of datasets\n",
  "  n_type12          : datasets with type-1/2 (order-statistic) summaries\n",
  "  note              : 'identical' when n_type12=0; blank otherwise\n",
  "  fact_*/joint_*    : posterior median (2.5%, 97.5%) for each model\n",
  "  delta_pred_median : joint minus factorised posterior median of pred_median\n",
  "  delta_pred_q95    : joint minus factorised posterior median of pred_q95\n"
)

# Print in two blocks: predicted quantities and parameters.
cols_pred   <- c("pathogen", "dist", "n_datasets", "n_type12",
                 "fact_pred_median", "joint_pred_median", "delta_pred_median",
                 "fact_pred_q95",    "joint_pred_q95",    "delta_pred_q95",
                 "note")
cols_params <- c("pathogen", "dist", "n_datasets", "n_type12",
                 "fact_mu0",  "joint_mu0",
                 "fact_phi",  "joint_phi",
                 "fact_tau",  "joint_tau",
                 "note")

message("\n--- Predicted quantities (days) ---\n")
print(comparison_tbl[, intersect(cols_pred, names(comparison_tbl))],
      row.names = FALSE, right = FALSE)

message("\n--- Parameter estimates ---\n")
print(comparison_tbl[, intersect(cols_params, names(comparison_tbl))],
      row.names = FALSE, right = FALSE)

message("\n", strrep("=", 80))
message("Done.")


# ── 10. LaTeX comparison table ────────────────────────────────────────────────
#
# Layout: two rows per (pathogen × distribution) pair — factorised on top,
# joint below.  This keeps the column count to 12, which fits comfortably in
# landscape A4 at \footnotesize with \tabcolsep = 4 pt.
#
# Columns:
#   Pathogen | Distribution | N | n_12 | Model | Median | p_95 |
#   mu_0 | phi | tau | Delta-Median | Delta-p_95
#
# Distribution, N, and n_12 use \multirow{2} to span both rows of each pair.
# Pathogen is shown in bold on the first factorised row of each group; the
# cell is left blank for subsequent rows (avoids variable-count multirow spans).
# A thin \cmidrule separates distribution pairs within a pathogen group;
# \midrule separates pathogen groups.
#
# Requires in the LaTeX preamble:
#   \usepackage{booktabs}
#   \usepackage{longtable}
#   \usepackage{multirow}
#   \usepackage{pdflscape}
#   \usepackage[table]{xcolor}   % for \rowcolor   % or lscape

.generate_comparison_latex_table <- function(
  tbl,
  pathogen_order = NULL,
  caption = paste0(
    "Comparison of posterior summaries under the factorised and joint ",
    "order-statistic likelihoods (\\textit{filtered} analysis). ",
    "Each pair of rows corresponds to one pathogen--distribution combination; ",
    "the upper row gives results under the factorised model and the lower row ",
    "under the joint model. ",
    "$N$: total number of datasets; ",
    "$n_{12}$: datasets contributing a type-1 (median\\,+\\,range) or ",
    "type-2 (median\\,+\\,IQR) summary, i.e.\\ those for which the two ",
    "likelihoods differ. ",
    "All posterior summaries are reported as median (2.5\\%ile, 97.5\\%ile). ",
    "$\\tau$ is shown as \\texttt{---} when $N < 5$ (prior-dominated). ",
    "$\\Delta$: joint minus factorised posterior median of the predictive ",
    "quantity (days). ",
    "$^{\\dagger}$\\,No type-1/2 datasets present; the two likelihoods are ",
    "identical and the joint row repeats the factorised estimates."
  ),
  label = "tab:comparison_factorised_joint"
) {

  # ── Helpers ─────────────────────────────────────────────────────────────────

  .esc <- function(x) {
    x <- gsub("_",  "\\_", x, fixed = TRUE)
    x <- gsub("&",  "\\&", x, fixed = TRUE)
    x <- gsub("%",  "\\%", x, fixed = TRUE)
    x
  }

  .cell <- function(x) {
    if (is.na(x) || x == "") return("---")
    x
  }

  .delta_fmt <- function(d) {
    if (is.na(d)) return("---")
    sprintf("%+.2f", d)
  }

  # ── Pathogen display labels ──────────────────────────────────────────────────

  plabels <- c(
    Nipah    = "Nipah",           MVD      = "Marburg (MVD)",
    EVD      = "Ebola (EVD)",     Lassa    = "Lassa fever",
    SARS     = "SARS",            MERS     = "MERS",
    Zika     = "Zika",            Measles  = "Measles",
    Mpox     = "Mpox",            Cholera  = "Cholera",
    RVF      = "Rift Valley fever", CCHF   = "CCHF",
    COVID_19 = "COVID-19",        Dengue   = "Dengue",
    YFV      = "Yellow fever",    Typhoid  = "Typhoid",
    Smallpox = "Smallpox",        Flu      = "Influenza"
  )

  if (is.null(pathogen_order))
    pathogen_order <- unique(tbl$pathogen)

  # ── Column spec ──────────────────────────────────────────────────────────────
  # Total nominal content width ≈ 194 mm; with \tabcolsep=4pt fits A4 landscape.

  col_spec <- paste0(
    "@{} ",
    "p{2.5cm} ",   # 1  Pathogen
    "p{1.6cm} ",   # 2  Distribution
    "r ",          # 3  N
    "r ",          # 4  n_12
    "l ",          # 5  Model
    "p{2.3cm} ",   # 6  Median (days)
    "p{2.3cm} ",   # 7  p_95 (days)
    "p{2.2cm} ",   # 8  mu_0
    "p{3.0cm} ",   # 9  phi
    "p{2.0cm} ",   # 10 tau
    "r ",          # 11 Delta Median
    "r ",          # 12 Delta p_95
    "@{}"
  )

  # ── Header rows ──────────────────────────────────────────────────────────────

  span_header <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    " & & & & ",
    "& \\multicolumn{5}{c}{Posterior summary --- median (95\\% CI)} ",
    "& \\multicolumn{2}{c}{$\\Delta$ (joint $-$ fact.)} \\\\"
  )
  span_cmidrule <- "\\cmidrule(lr){6-10} \\cmidrule(lr){11-12}"

  col_header <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    "Pathogen & Distribution & $N$ & $n_{12}$ & Model ",
    "& Median (days) & $p_{95}$ (days) & $\\mu_0$ & $\\phi$ & $\\tau$ ",
    "& Median & $p_{95}$ \\\\"
  )

  make_head <- function(is_first) {
    c(
      if (!is_first)
        "\\multicolumn{12}{l}{\\small\\textit{continued from previous page}} \\\\",
      "\\toprule",
      span_header,
      span_cmidrule,
      col_header,
      "\\midrule"
    )
  }

  # ── Table body ───────────────────────────────────────────────────────────────

  body <- character(0)
  first_pathogen <- TRUE

  for (pathogen in pathogen_order) {

    sub <- tbl[tbl$pathogen == pathogen, , drop = FALSE]
    if (nrow(sub) == 0L) next

    plabel <- if (!is.na(plabels[pathogen])) plabels[[pathogen]] else pathogen

    if (!first_pathogen)
      body <- c(body, "\\midrule")
    first_pathogen <- FALSE

    for (i in seq_len(nrow(sub))) {

      row <- sub[i, ]

      is_identical  <- grepl("identical",  row$note, fixed = FALSE)
      joint_unavail <- grepl("unavailable", row$note, fixed = FALSE)

      # ── Factorised row ──────────────────────────────────────────────────────

      # Pathogen label: bold on first dist row only; blank thereafter.
      path_cell <- if (i == 1L)
        paste0("\\textbf{", .esc(plabel), "}")
      else
        ""

      fact_cells <- paste(
        path_cell,
        paste0("\\multirow{2}{*}{", .esc(row$dist), "}"),
        paste0("\\multirow{2}{*}{", row$n_datasets, "}"),
        paste0("\\multirow{2}{*}{", row$n_type12,   "}"),
        "Factorised",
        .cell(row$fact_pred_median),
        .cell(row$fact_pred_q95),
        .cell(row$fact_mu0),
        .cell(row$fact_phi),
        .cell(row$fact_tau),
        "",   # Delta Median: blank on factorised row
        "",   # Delta p_95:   blank on factorised row
        sep = " & "
      )
      body <- c(body, paste0(fact_cells, " \\\\"))

      # ── Joint row ──────────────────────────────────────────────────────────

      model_lbl <- if (is_identical)  "Joint$^{\\dagger}$"
                   else if (joint_unavail) "Joint (failed)"
                   else               "Joint"

      j_median <- if (joint_unavail) "---" else .cell(row$joint_pred_median)
      j_q95    <- if (joint_unavail) "---" else .cell(row$joint_pred_q95)
      j_mu0    <- if (joint_unavail) "---" else .cell(row$joint_mu0)
      j_phi    <- if (joint_unavail) "---" else .cell(row$joint_phi)
      j_tau    <- if (joint_unavail) "---" else .cell(row$joint_tau)

      d_median <- if (joint_unavail || is_identical) "---"
                  else .delta_fmt(row$delta_pred_median)
      d_q95    <- if (joint_unavail || is_identical) "---"
                  else .delta_fmt(row$delta_pred_q95)

      joint_cells <- paste(
        "",          # Pathogen: blank (multirow from fact row covers it)
        "",          # Dist:     blank
        "",          # N:        blank
        "",          # n_12:     blank
        model_lbl,
        j_median, j_q95, j_mu0, j_phi, j_tau,
        d_median, d_q95,
        sep = " & "
      )
      body <- c(body, paste0(joint_cells, " \\\\"))

      # Thin rule between distribution pairs within the same pathogen group,
      # but not after the last distribution for this pathogen.
      if (i < nrow(sub))
        body <- c(body, "\\cmidrule(l{4pt}r{0pt}){2-12}")

    }
  }

  # ── Assemble full table ──────────────────────────────────────────────────────

  lines <- c(
    "% ============================================================",
    "% Comparison table: factorised vs joint order-statistic likelihood",
    "% Required LaTeX packages: booktabs, longtable, multirow, pdflscape, xcolor (table)",
    "% ============================================================",
    "{\\setlength{\\tabcolsep}{4pt}",
    "\\begin{landscape}",
    "\\footnotesize",
    paste0("\\begin{longtable}{", col_spec, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
    make_head(is_first = TRUE),
    "\\endfirsthead",
    make_head(is_first = FALSE),
    "\\endhead",
    "\\multicolumn{12}{r}{\\small\\textit{continued on next page}} \\\\",
    "\\endfoot",
    "\\bottomrule",
    paste0(
      "\\multicolumn{12}{l}{\\footnotesize",
      " $^{\\dagger}$\\,No type-1/2 datasets; joint and factorised",
      " likelihoods are identical.} \\\\"
    ),
    "\\endlastfoot",
    body,
    "\\end{longtable}",
    "\\end{landscape}",
    "}  % end \\setlength{\\tabcolsep}"
  )

  paste(lines, collapse = "\n")
}


# ── Generate, write, and echo the table ──────────────────────────────────────

latex_table <- .generate_comparison_latex_table(comparison_tbl)

latex_file <- file.path(OUTPUT_DIR, "comparison_factorised_vs_joint.tex")
writeLines(latex_table, latex_file)
message("\nLaTeX table written to: ", latex_file)

# Echo to console so it can be copied directly.
message("\n", strrep("-", 80))
message("LaTeX table (copy into document):\n")
cat(latex_table, "\n")
message(strrep("-", 80))
