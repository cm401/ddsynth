# =============================================================================
# table_utils.R
# -----------------------------------------------------------------------------
# Functions for generating publication-ready LaTeX summary tables from the
# incubation period analysis results.
# =============================================================================


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
               "$-$"
  )
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
  model_weights   = NULL,
  show_subgroups  = TRUE,
  rhat_threshold  = 1.05,
  ci_prob         = 0.95,
  pathogen_labels = NULL,
  caption = paste0(
    "Incubation period estimates by pathogen. ",
    "Values are posterior medians with 95\\% credible intervals (CI). ",
    "Subgroup analyses are shown in italics beneath each pathogen."
  ),
  label = "tab:incubation_periods"
) {

  # ── Setup ───────────────────────────────────────────────────────────────────
  probs  <- c((1 - ci_prob) / 2, 0.5, 1 - (1 - ci_prob) / 2)

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

  # ── Build row data ───────────────────────────────────────────────────────────
  # Each element of `table_rows` is a list with:
  #   $label       display name (character)
  #   $dist        distribution label (character; "" for subgroup rows)
  #   $median      formatted string
  #   $q95         formatted string
  #   $mu0         formatted string  (population-level location parameter)
  #   $phi         formatted string
  #   $kappa       formatted string  ("$-$" when not applicable)
  #   $is_subgroup logical

  table_rows <- list()

  for (pathogen in names(all_results)) {
    plabel  <- p_labels[[pathogen]]
    if (is.null(plabel)) plabel <- pathogen

    weights <- model_weights[[pathogen]]
    if (is.null(weights) || length(weights) == 0L) next

    # Walk weight-ranked distributions to find a converged fit
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

    # ── Subgroup sub-rows ────────────────────────────────────────────────────
    if (show_subgroups) {
      sg_keys <- setdiff(names(all_results[[pathogen]]), c("all", "filtered"))
      for (sg in sg_keys) {
        sg_r <- all_results[[pathogen]][[sg]][[best_dist]]
        if (is.null(sg_r) || isTRUE(sg_r$skipped) || is.null(sg_r$fit)) next
        if (!.fit_has_converged(sg_r$fit, rhat_threshold)) next

        sg_cells <- .extract_fit_cells(sg_r$fit, best_dist, probs)
        if (is.null(sg_cells)) next

        table_rows[[length(table_rows) + 1L]] <- c(
          list(label = sg, dist = "", is_subgroup = TRUE),
          sg_cells
        )
      }
    }
  }

  if (length(table_rows) == 0L)
    stop("No rows could be built — check that all_results contains valid fits.")

  # ── Assemble LaTeX ──────────────────────────────────────────────────────────
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
    "\\midrule"
  )

  for (i in seq_along(table_rows)) {
    r <- table_rows[[i]]

    # Vertical space before every new pathogen block (not before subgroup rows)
    if (i > 1L && !r$is_subgroup)
      lines <- c(lines, "\\addlinespace[4pt]")

    # Pathogen / subgroup cell
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
      sep = " & "
    )
    lines <- c(lines, paste0(row_str, " \\\\"))
  }

  lines <- c(lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  paste(lines, collapse = "\n")
}
