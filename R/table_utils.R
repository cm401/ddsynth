# =============================================================================
# table_utils.R
# -----------------------------------------------------------------------------
# Functions for generating publication-ready LaTeX summary tables from the
# incubation period analysis results.
# =============================================================================


# ── Group-mapping constants ───────────────────────────────────────────────────

# Maps all_results keys (which use underscores) to the six epidemiological
# group names used in the map (map_utils.R .pathogen_group_map).
.RESULT_KEY_TO_GROUP <- c(
  EVD           = "Viral haemorrhagic fevers",
  MVD           = "Viral haemorrhagic fevers",
  Lassa         = "Viral haemorrhagic fevers",
  CCHF          = "Viral haemorrhagic fevers",
  COVID_19      = "Pandemic respiratory",
  SARS          = "Pandemic respiratory",
  MERS          = "Pandemic respiratory",
  Flu           = "Pandemic respiratory",
  Dengue        = "Arboviral / vector-borne",
  Zika          = "Arboviral / vector-borne",
  RVF           = "Arboviral / vector-borne",
  YFV           = "Arboviral / vector-borne",
  Nipah         = "Bat-reservoir zoonoses",
  Mpox          = "Human-to-human viral",
  Measles       = "Human-to-human viral",
  Smallpox      = "Human-to-human viral",
  Cholera       = "Environmental / zoonotic bacterial",
  Typhoid       = "Environmental / zoonotic bacterial"
)

.GROUP_ORDER <- c(
  "Viral haemorrhagic fevers",
  "Pandemic respiratory",
  "Arboviral / vector-borne",
  "Bat-reservoir zoonoses",
  "Human-to-human viral",
  "Environmental / zoonotic bacterial"
)

# Human-readable display labels for subgroup keys (the named keys used in
# subgroup_config in main_analysis.R).  Any key not listed here falls back
# to the raw key name.
.SUBGROUP_LABELS <- c(
  # Nipah
  Bangladesh    = "Bangladesh",
  # EVD
  WestAfrica    = "West Africa",
  nonWestAfrica = "non-West Africa",
  # SARS
  HongKong      = "Hong Kong",
  Canada        = "Canada",
  Taiwan        = "Taiwan",
  China         = "China",
  Singapore     = "Singapore",
  # MERS
  SaudiArabia   = "Saudi Arabia",
  SouthKorea    = "South Korea",
  # Cholera
  O1_El_Tor     = "O1 El Tor",
  O1_Classical  = "O1 Classical",
  # CCHF
  tick_bite     = "Tick-bite",
  nosocomial    = "Nosocomial",
  # COVID-19
  Wildtype      = "Wildtype",
  Alpha         = "Alpha",
  Delta         = "Delta",
  Omicron       = "Omicron",
  # Dengue
  inoculation   = "Inoculation",
  mosquito      = "Mosquito bite",
  # Typhoid
  CateredMeal   = "Catered meal",
  Experimental  = "Experimental",
  Other         = "Other",
  # Flu
  H1N1          = "H1N1",
  H3N2          = "H3N2",
  H5N1          = "H5N1",
  InfluenzaA    = "Influenza A",
  InfluenzaB    = "Influenza B"
)


# ── Internal helpers ──────────────────────────────────────────────────────────

# Escape characters that have special meaning in LaTeX.
.escape_latex <- function(x) {
  # Order matters: backslash must come first.
  x <- gsub("\\\\", "\\\\textbackslash{}", x, fixed = FALSE)
  x <- gsub("&",    "\\&",  x, fixed = TRUE)
  x <- gsub("%",    "\\%",  x, fixed = TRUE)
  x <- gsub("\\$",  "\\$",  x, fixed = FALSE)
  x <- gsub("#",    "\\#",  x, fixed = TRUE)
  x <- gsub("_",    "\\_",  x, fixed = TRUE)
  x
}

# Format a numeric vector as  "median (lo, hi)"  using the supplied quantile
# probabilities.  Returns "---" if x is NULL or all-NA.
.fmt_est <- function(x, probs, digits) {
  if (is.null(x)) return("---")
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("---")
  q <- stats::quantile(x, probs, na.rm = TRUE)
  sprintf("%.*f (%.*f, %.*f)", digits, q[2L], digits, q[1L], digits, q[3L])
}

# Extract the seven table cells from a single stanfit object.
.extract_fit_cells <- function(fit, dist_name, probs) {
  sims <- tryCatch(rstan::extract(fit), error = function(e) NULL)
  if (is.null(sims)) return(NULL)

  list(
    median = .fmt_est(sims$pred_median, probs, digits = 1L),
    q95    = .fmt_est(sims$pred_q95,    probs, digits = 1L),
    mu0    = .fmt_est(sims$mu0,         probs, digits = 2L),
    phi    = .fmt_est(sims$phi,         probs, digits = 2L),
    kappa  = if (dist_name %in% c("burr", "gengamma") && !is.null(sims$kappa))
               .fmt_est(sims$kappa, probs, digits = 2L)
             else
               "$-$",
    tau    = .fmt_est(sims$tau,         probs, digits = 2L)
  )
}


# Build the table_rows list for the pathogens present in `all_results`.
# `all_results` should already be filtered to the desired subset before calling.
.build_table_rows <- function(
  all_results,
  model_weights,
  show_subgroups,
  rhat_threshold,
  probs,
  p_labels
) {
  table_rows <- list()

  for (pathogen in names(all_results)) {
    plabel  <- p_labels[[pathogen]]
    if (is.null(plabel)) plabel <- pathogen

    weights <- model_weights[[pathogen]]
    if (is.null(weights) || length(weights) == 0L) next

    best_dist <- NULL
    best_res  <- NULL
    for (.cand in names(weights)) {
      .r <- all_results[[pathogen]][["filtered"]][[.cand]]
      if (is.null(.r) || isTRUE(.r$skipped) || is.null(.r$fit)) next
      if (!.fit_has_converged(.r$fit, rhat_threshold)) next
      best_dist <- .cand
      best_res  <- .r
      break
    }
    if (is.null(best_dist)) next

    cells <- .extract_fit_cells(best_res$fit, best_dist, probs)
    if (is.null(cells)) next

    table_rows[[length(table_rows) + 1L]] <- c(
      list(label = plabel, dist = .DIST_LABELS[[best_dist]], is_subgroup = FALSE),
      cells
    )

    if (show_subgroups) {
      sg_keys <- setdiff(names(all_results[[pathogen]]), c("all", "filtered"))
      for (sg in sg_keys) {
        sg_r <- all_results[[pathogen]][[sg]][[best_dist]]
        if (is.null(sg_r) || isTRUE(sg_r$skipped) || is.null(sg_r$fit)) next
        if (!.fit_has_converged(sg_r$fit, rhat_threshold)) next

        sg_cells <- .extract_fit_cells(sg_r$fit, best_dist, probs)
        if (is.null(sg_cells)) next

        sg_label <- if (!is.na(.SUBGROUP_LABELS[sg])) .SUBGROUP_LABELS[[sg]] else sg
        table_rows[[length(table_rows) + 1L]] <- c(
          list(label = sg_label, dist = "", is_subgroup = TRUE),
          sg_cells
        )
      }
    }
  }

  table_rows
}

# Render a list of table_rows into LaTeX lines.
# group_label: optional bold letter label emitted as a spanning header row.
# first_group: if TRUE, suppress the leading \midrule (header already present).
.render_group_rows <- function(table_rows, group_label = NULL, first_group = FALSE,
                               coloured = FALSE) {
  lines <- character(0)

  if (!first_group)
    lines <- c(lines, "\\midrule")

  if (!is.null(group_label)) {
    if (coloured) {
      lines <- c(lines,
        paste0("\\rowcolor{gray!15}",
               "\\multicolumn{8}{@{}l}{\\textbf{\\textit{",
               .escape_latex(group_label), "}}} \\\\"))
    } else {
      lines <- c(lines,
        paste0("\\multicolumn{8}{@{}l}{\\textbf{", group_label, "}} \\\\[2pt]"))
    }
  }

  for (i in seq_along(table_rows)) {
    r <- table_rows[[i]]

    if (i > 1L && !r$is_subgroup)
      lines <- c(lines, "\\addlinespace[4pt]")

    name_cell <- if (r$is_subgroup) {
      paste0("\\quad \\textit{", .escape_latex(r$label), "}")
    } else {
      paste0("\\textbf{", .escape_latex(r$label), "}")
    }

    row_str <- paste(
      name_cell,
      .escape_latex(r$dist),
      r$median,
      r$q95,
      r$mu0,
      r$phi,
      r$kappa,
      r$tau,
      sep = " & "
    )
    lines <- c(lines, paste0(row_str, " \\\\"))
  }

  lines
}


# ── Exported function ─────────────────────────────────────────────────────────

#' Generate a LaTeX summary table of incubation period estimates
#'
#' Produces a \pkg{booktabs}-style LaTeX table with one row per pathogen and
#' (optionally) indented italic sub-rows for each subgroup analysis.
#'
#' @details
#' For every row the table reports:
#' \itemize{
#'   \item the best-fitting distribution (selected with the same
#'         convergence-filtered LOO-weight logic as [plot_main_figure()]);
#'   \item the posterior predictive **median** with 95% credible interval;
#'   \item the posterior predictive **95th percentile** (\eqn{p_{95}}) with 95% CI;
#'   \item the population-level location parameter \eqn{\mu_0} with 95% CI;
#'   \item the shape parameter \eqn{\phi} with 95% CI;
#'   \item the second shape parameter \eqn{\kappa} with 95% CI (Burr XII and
#'         Generalised Gamma only; shown as \eqn{-} for other families).
#' }
#' All intervals are formatted as \samp{median (lower, upper)}.
#'
#' The table requires the \LaTeX{} packages \pkg{booktabs} and
#' \pkg{multirow} in the preamble.
#'
#' @param all_results Nested list produced by \code{analysis/main_analysis.R}.
#' @param model_weights Pre-computed output from
#'   [compute_pathogen_model_bayes_factors()].  Computed internally if
#'   \code{NULL} (slow — pre-compute and cache).
#' @param show_subgroups Logical.  Include subgroup sub-rows (default
#'   \code{TRUE}).
#' @param rhat_threshold Numeric.  Maximum acceptable Rhat for convergence;
#'   fits with max Rhat above this are skipped (default 1.05).
#' @param ci_prob Numeric.  Width of the credible interval, e.g. 0.95 gives
#'   a 2.5%–97.5% interval (default 0.95).
#' @param pathogen_labels Optional named character vector overriding the
#'   display names used in the table rows.
#' @param caption LaTeX \code{\\caption\{\}} string.
#' @param label LaTeX \code{\\label\{\}} string.
#'
#' @return A single character string containing the complete LaTeX table
#'   environment, ready to paste into an Overleaf document.
#'
#' @seealso [compute_pathogen_model_bayes_factors()], [plot_main_figure()]
#' @export
generate_results_table <- function(
  all_results,
  pathogens       = NULL,
  model_weights   = NULL,
  show_subgroups  = TRUE,
  rhat_threshold  = 1.05,
  ci_prob         = 0.95,
  pathogen_labels = NULL,
  caption = paste0(
    "Incubation period estimates by pathogen, grouped by epidemiological category. ",
    "Values are posterior medians with 95\\% credible intervals (CI). ",
    "$\\tau$ is the between-study heterogeneity SD. ",
    "Subgroup analyses are shown in italics beneath each pathogen."
  ),
  label = "tab:incubation_periods"
) {

  # ── Setup ───────────────────────────────────────────────────────────────────
  if (!is.null(pathogens))
    all_results <- all_results[intersect(pathogens, names(all_results))]

  probs    <- c((1 - ci_prob) / 2, 0.5, 1 - (1 - ci_prob) / 2)
  p_labels <- .PATHOGEN_LABELS
  if (!is.null(pathogen_labels))
    p_labels[names(pathogen_labels)] <- pathogen_labels

  if (is.null(model_weights)) {
    message(
      "No model_weights supplied — running ",
      "compute_pathogen_model_bayes_factors() now.\n",
      "Consider pre-computing: mw <- compute_pathogen_model_bayes_factors(all_results)"
    )
    model_weights <- compute_pathogen_model_bayes_factors(all_results)
  }

  # ── Build row data grouped by epidemiological category ───────────────────────
  body_lines <- character(0L)
  first_g    <- NA_integer_

  for (g in seq_along(.GROUP_ORDER)) {
    grp_name <- .GROUP_ORDER[g]
    grp_keys <- names(.RESULT_KEY_TO_GROUP)[.RESULT_KEY_TO_GROUP == grp_name]
    grp_keys <- intersect(grp_keys, names(all_results))
    if (length(grp_keys) == 0L) next

    sub  <- all_results[grp_keys]
    rows <- .build_table_rows(sub, model_weights, show_subgroups,
                              rhat_threshold, probs, p_labels)
    if (length(rows) == 0L) next

    is_first <- is.na(first_g)
    if (is_first) first_g <- g

    body_lines <- c(body_lines,
      .render_group_rows(rows, group_label = grp_name, first_group = is_first,
                         coloured = TRUE))
  }

  if (length(body_lines) == 0L)
    stop("No rows could be built — check that all_results contains valid fits.")

  # ── Shared header row (reused in first-head and continuation head) ───────────
  col_header <- paste0(
    "\\rowcolor[HTML]{F0E68C}",
    "Pathogen & Distribution & ",
    "\\multicolumn{1}{c}{Median (days)} & ",
    "\\multicolumn{1}{c}{$p_{95}$ (days)} & ",
    "\\multicolumn{1}{c}{$\\mu_0$} & ",
    "\\multicolumn{1}{c}{$\\phi$} & ",
    "\\multicolumn{1}{c}{$\\kappa$} & ",
    "\\multicolumn{1}{c}{$\\tau$} \\\\"
  )

  # ── Assemble LaTeX ──────────────────────────────────────────────────────────
  lines <- c(
    "% Requires \\usepackage{booktabs}, \\usepackage{longtable},",
    "% \\usepackage[table]{xcolor} and \\usepackage{pdflscape} in preamble.",
    "\\begin{landscape}",
    "\\small",
    paste0("\\begin{longtable}{@{} l l r r r r r r @{}}"),
    paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
    "\\toprule",
    col_header,
    "\\midrule",
    "\\endfirsthead",
    "%",
    "\\multicolumn{8}{l}{\\small\\textit{continued from previous page}} \\\\",
    "\\toprule",
    col_header,
    "\\midrule",
    "\\endhead",
    "%",
    "\\multicolumn{8}{r}{\\small\\textit{continued on next page}} \\\\",
    "\\endfoot",
    "%",
    "\\bottomrule",
    "\\endlastfoot",
    "%",
    body_lines,
    "\\end{longtable}",
    "\\end{landscape}"
  )

  paste(lines, collapse = "\n")
}


#' Generate a LaTeX summary table split into three panels by data availability
#'
#' Produces a single \pkg{booktabs}-style LaTeX table with the same columns as
#' [generate_results_table()], but with pathogens divided into three labelled
#' panels (A, B, C) separated by horizontal rules and bold panel-label header
#' rows.
#'
#' @param all_results Nested list produced by \code{analysis/main_analysis.R}.
#' @param group_a Character vector of pathogen keys for panel A (high data).
#' @param group_b Character vector of pathogen keys for panel B (moderate data).
#' @param group_c Character vector of pathogen keys for panel C (low data).
#' @param panel_labels Length-3 character vector of panel labels (default
#'   \code{c("A", "B", "C")}).
#' @param model_weights Pre-computed output from
#'   [compute_pathogen_model_bayes_factors()].  Computed once internally if
#'   \code{NULL}.
#' @param show_subgroups Logical.  Include subgroup sub-rows (default
#'   \code{TRUE}).
#' @param rhat_threshold Numeric.  Maximum acceptable Rhat (default 1.05).
#' @param ci_prob Numeric.  Width of the credible interval (default 0.95).
#' @param pathogen_labels Optional named character vector overriding display
#'   names.
#' @param caption LaTeX \code{\\caption\{\}} string.
#' @param label LaTeX \code{\\label\{\}} string.
#'
#' @return A single character string containing the complete LaTeX table
#'   environment.
#'
#' @seealso [generate_results_table()], [plot_main_figure_split()]
#' @export
generate_results_table_split <- function(
  all_results,
  group_a        = c("SARS", "COVID_19", "Measles", "Mpox"),
  group_b        = c("Nipah", "EVD", "MERS", "Cholera",
                     "CCHF", "Dengue", "YFV"),
  group_c        = c("MVD", "Lassa", "Zika", "RVF"),
  panel_labels   = c("A", "B", "C"),
  model_weights  = NULL,
  show_subgroups = TRUE,
  rhat_threshold = 1.05,
  ci_prob        = 0.95,
  pathogen_labels = NULL,
  caption = paste0(
    "Incubation period estimates by pathogen. ",
    "Values are posterior medians with 95\\% credible intervals (CI). ",
    "Pathogens are grouped by data availability: panel A (substantial data), ",
    "panel B (moderate data), and panel C (limited data). ",
    "Subgroup analyses are shown in italics beneath each pathogen."
  ),
  label = "tab:incubation_periods"
) {

  # ── Setup ───────────────────────────────────────────────────────────────────
  probs    <- c((1 - ci_prob) / 2, 0.5, 1 - (1 - ci_prob) / 2)
  p_labels <- .PATHOGEN_LABELS
  if (!is.null(pathogen_labels))
    p_labels[names(pathogen_labels)] <- pathogen_labels

  # Pre-compute model weights once — shared across all three groups
  if (is.null(model_weights)) {
    message(
      "No model_weights supplied — running ",
      "compute_pathogen_model_bayes_factors() now.\n",
      "Consider pre-computing: mw <- compute_pathogen_model_bayes_factors(all_results)"
    )
    model_weights <- compute_pathogen_model_bayes_factors(all_results)
  }

  # ── Build row data per group ─────────────────────────────────────────────────
  groups <- list(group_a, group_b, group_c)
  rows_per_group <- lapply(groups, function(grp) {
    sub <- all_results[intersect(grp, names(all_results))]
    .build_table_rows(sub, model_weights, show_subgroups, rhat_threshold,
                      probs, p_labels)
  })

  if (all(vapply(rows_per_group, length, integer(1L)) == 0L))
    stop("No rows could be built — check that all_results contains valid fits.")

  # ── Assemble LaTeX ──────────────────────────────────────────────────────────
  body_lines <- unlist(
    mapply(
      .render_group_rows,
      table_rows  = rows_per_group,
      group_label = as.list(panel_labels),
      first_group = list(TRUE, FALSE, FALSE),
      SIMPLIFY    = FALSE
    ),
    use.names = FALSE
  )

  lines <- c(
    "% Requires \\usepackage{booktabs} in preamble.",
    "\\begin{table}[ht]",
    "\\centering",
    "\\small",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{@{} l l r r r r r @{}}",
    "\\toprule",
    paste0(
      "Pathogen & Distribution & ",
      "\\multicolumn{1}{c}{Median (days)} & ",
      "\\multicolumn{1}{c}{$p_{95}$ (days)} & ",
      "\\multicolumn{1}{c}{$\\mu_0$} & ",
      "\\multicolumn{1}{c}{$\\phi$} & ",
      "\\multicolumn{1}{c}{$\\kappa$} \\\\"
    ),
    "\\midrule",
    body_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  paste(lines, collapse = "\n")
}


# ── Internal helpers (data table) ────────────────────────────────────────────

#' @noRd
.fmt_num <- function(x, digits = 1L) {
  if (is.null(x) || length(x) == 0L || !is.finite(x)) return("---")
  sprintf("%.*f", digits, x)
}

#' @noRd
.fmt_citation <- function(source) {
  if (is.null(source) || !nzchar(source)) return("---")
  m <- regmatches(source, regexpr("^[^,]+\\(\\d{4}\\)", source))
  if (length(m) == 1L && nzchar(m)) return(trimws(m))
  if (nchar(source) > 40L) paste0(substr(source, 1L, 37L), "...") else source
}

#' @noRd
.fmt_doi <- function(source) {
  if (is.null(source) || !nzchar(source)) return("---")
  m <- regmatches(source,
    regexec("doi:\\s*(?:doi\\.org/)?(10\\.[^,\\s]+)", source, perl = TRUE))
  if (length(m[[1L]]) >= 2L && nzchar(m[[1L]][2L])) return(m[[1L]][2L])
  "---"
}

#' @noRd
.detect_type <- function(d) {
  if (!is.null(d$freq_lower) && !is.null(d$freq_upper)) return("E")
  if (!is.null(d$freq_value) && !is.null(d$freq_count)) return("D")
  if (!is.null(d$mean)       && !is.null(d$sd))         return("C")
  if (!is.null(d$median)     && !is.null(d$Q1))         return("B")
  if (!is.null(d$median)     && !is.null(d$min))        return("A")
  "unknown"
}

#' @noRd
.derive_n <- function(d) {
  if (!is.null(d$n) && length(d$n) == 1L && is.finite(d$n)) return(as.integer(d$n))
  if (!is.null(d$freq_count)) return(as.integer(sum(d$freq_count)))
  NA_integer_
}

#' @noRd
.weighted_median <- function(values, counts) {
  stats::median(rep(values, times = as.integer(counts)))
}

#' @noRd
.fmt_summary_cell <- function(d, digits = 1L) {
  type <- .detect_type(d)
  if (type == "A") {
    sprintf("Median: %s (Range: %s--%s)",
      .fmt_num(d$median, digits), .fmt_num(d$min, digits), .fmt_num(d$max, digits))
  } else if (type == "B") {
    sprintf("Median: %s (IQR: %s--%s)",
      .fmt_num(d$median, digits), .fmt_num(d$Q1, digits), .fmt_num(d$Q3, digits))
  } else if (type == "C") {
    sprintf("Mean: %s (SD: %s)",
      .fmt_num(d$mean, digits), .fmt_num(d$sd, digits))
  } else if (type == "D") {
    wmed <- .weighted_median(d$freq_value, d$freq_count)
    sprintf("\\textit{Freq. table.} Median: %s (Range: %s--%s)",
      .fmt_num(wmed, digits),
      .fmt_num(min(d$freq_value), digits),
      .fmt_num(max(d$freq_value), digits))
  } else if (type == "E") {
    sprintf("\\textit{Freq. table (interval-censored).} Range: %s--%s",
      .fmt_num(min(d$freq_lower), digits),
      .fmt_num(max(d$freq_upper), digits))
  } else {
    "---"
  }
}

#' @noRd
.fmt_data_row <- function(d, digits = 1L) {
  ref      <- .escape_latex(.fmt_citation(d$source))
  country  <- .escape_latex(if (!is.null(d$country)  && nzchar(d$country))  d$country  else "---")
  subgroup <- .escape_latex(if (!is.null(d$subgroup) && nzchar(d$subgroup)) d$subgroup else "---")
  n_val    <- .derive_n(d)
  n_str    <- if (is.na(n_val)) "---" else as.character(n_val)
  summary  <- .fmt_summary_cell(d, digits)  # not escaped; all content is controlled
  doi      <- .escape_latex(.fmt_doi(d$source))
  paste(ref, country, subgroup, n_str, summary, doi, sep = " & ")
}

#' @noRd
.fmt_pathogen_header <- function(label) {
  c(paste0("\\rowcolor{gray!15}\\multicolumn{6}{|l|}{\\textbf{\\textit{",
           .escape_latex(label), "}}} \\\\"),
    "\\noalign{\\arrayrulecolor{gray}}",
    "\\hline",
    "\\noalign{\\arrayrulecolor{black}}")
}


# ── Exported function ─────────────────────────────────────────────────────────

#' Generate a multi-page landscape LaTeX data summary table
#'
#' Produces a landscape \pkg{longtable} listing all collected incubation period
#' observations across pathogens, with one bold section header row per
#' pathogen, a khaki-coloured column header, and a DOI column.
#'
#' @details
#' Each dataset entry produces one table row.  The Summary column reports:
#' \itemize{
#'   \item \strong{Type A} (median + range): \code{"Median: X (Range: lo--hi)"}
#'   \item \strong{Type B} (median + IQR): \code{"Median: X (IQR: Q1--Q3)"}
#'   \item \strong{Type C} (mean + SD): \code{"Mean: X (SD: s)"}
#'   \item \strong{Type D} (frequency table): prefixed \emph{Freq. table.},
#'         then weighted median and observed range
#'   \item \strong{Type E} (interval-censored frequency table): prefixed
#'         \emph{Freq. table (interval-censored).}, then overall range
#' }
#' For Type D entries where \code{n} is absent, \eqn{n} is computed as
#' \code{sum(freq_count)}.
#'
#' Required \LaTeX{} packages (add to preamble):
#' \preformatted{
#'   \usepackage{longtable}
#'   \usepackage{array}
#'   \usepackage[table]{xcolor}
#'   \usepackage{pdflscape}
#'   \newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}
#'   \newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}
#' }
#'
#' @param datasets A named list of pathogen datasets.  Each element must be a
#'   named list of entry lists (e.g.
#'   \code{list(Nipah = datasets_Nipah, MVD = datasets_MVD)}).
#'   Pathogens appear in the order of \code{names(datasets)}.
#' @param caption LaTeX \code{\\caption\{\}} string.
#' @param label   LaTeX \code{\\label\{\}} string.
#' @param digits  Integer.  Decimal places for numeric summaries (default 1).
#'
#' @return A single character string containing the complete landscape
#'   \code{longtable} environment, ready to paste into a manuscript.
#'
#' @seealso [generate_results_table()], [generate_results_table_split()]
#' @export
generate_data_table <- function(
  datasets,
  caption = "Collected incubation period datasets by pathogen.",
  label   = "tab:data_summary",
  digits  = 1L
) {
  stopifnot(is.list(datasets), length(datasets) > 0L, !is.null(names(datasets)))

  col_spec <- "|C{4cm}|C{2.0cm}|C{2.5cm}|R{0.7cm}|C{5.5cm}|C{6cm}|"

  header_cells <- paste(
    "\\textbf{Reference}",
    "\\textbf{Country}",
    "\\textbf{Subgroup}",
    "\\textbf{$n$}",
    "\\textbf{Summary}",
    "\\textbf{DOI}",
    sep = " & "
  )
  header_row <- paste0(
    "\\rowcolor[HTML]{F0E68C}", header_cells, " \\\\"
  )

  cont_head <- c(
    "\\multicolumn{6}{l}{continued from previous page} \\\\",
    "\\hline",
    header_row,
    "\\hline",
    "\\endhead"
  )

  body_lines <- character(0L)
  pathogen_names <- names(datasets)
  for (nm in pathogen_names) {
    entries    <- datasets[[nm]]
    n_entries  <- length(entries)
    body_lines <- c(body_lines, .fmt_pathogen_header(nm))
    for (j in seq_along(entries)) {
      # dotted line between rows within a section; solid line after the last
      # row (acts as top border of the next section header, or closes the table)
      sep <- if (j < n_entries) "\\hdashline[1pt/1pt]" else "\\hline"
      body_lines <- c(body_lines,
        paste0(.fmt_data_row(entries[[j]], digits), " \\\\"),
        sep)
    }
  }

  lines <- c(
    "% Requires in preamble:",
    "%   \\usepackage{longtable}, \\usepackage{array},",
    "%   \\usepackage[table]{xcolor}, \\usepackage{pdflscape},",
    "%   \\usepackage{arydshln}   % for \\hdashline",
    "% And define column types once (preamble or before first use):",
    "%   \\newcolumntype{C}[1]{>{\\centering\\arraybackslash}p{#1}}",
    "%   \\newcolumntype{R}[1]{>{\\raggedleft\\arraybackslash}p{#1}}",
    "\\begin{landscape}",
    paste0("\\begin{longtable}{", col_spec, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "} \\\\"),
    "\\hline",
    header_row,
    "\\hline",
    "\\endfirsthead",
    "%",
    cont_head,
    "%",
    "\\multicolumn{6}{r}{continued on next page} \\\\",
    "\\endfoot",
    "\\hline",
    "\\endlastfoot",
    "%",
    body_lines,
    "\\end{longtable}",
    "\\end{landscape}"
  )

  paste(lines, collapse = "\n")
}
