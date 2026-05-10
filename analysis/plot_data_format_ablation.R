# =============================================================================
# plot_data_format_ablation.R
# -----------------------------------------------------------------------------
# Visualises the data-format ablation results produced by
# data_format_ablation.R.
#
# Three figures, each answering a distinct question:
#
#   Figure 1 — "Do the arms agree?"
#     Three-arm forest plot.  Pathogens on the y-axis, pred_median on a
#     log-scaled x-axis.  Arms A, B, C dodged vertically per pathogen, with
#     horizontal 95% CrI bars.  A vertical reference line marks the federated
#     (arm C) point estimate.  Faceted by distribution family.
#
#   Figure 2 — "How much does federation help?"
#     Side-by-side metrics strip using the same pathogen ordering.
#     Left panel  : interval ratio  (A_width / C_width, B_width / C_width);
#                   reference line at 1 (equal precision).
#     Middle panel: Jensen-Shannon divergence  (JS_CA, JS_CB; bits).
#     Right panel : overlap coefficient  (OVL_CA, OVL_CB).
#     All three panels share the y-axis (pathogens) and are assembled with
#     patchwork.  Restricted to one distribution (PRIMARY_DIST) for clarity;
#     a supplementary version facets over all distributions.
#
#   Figure 3 — "What does the gain look like in practice?"
#     Posterior predictive density overlays for the TOP_N_DENSITY pathogens
#     with the largest JS_CA divergence (i.e. where individual-level data
#     pulls the federated estimate furthest from the summary-stat arm).
#     Three shaded ribbons (A, B, C) per pathogen panel, computed from
#     compute_predictive_cdf().
#
# Output
# ------
#   results/figures/fig_ablation_forest.pdf  (+ .png)
#   results/figures/fig_ablation_gain.pdf    (+ .png)
#   results/figures/fig_ablation_densities.pdf (+ .png)
# =============================================================================

devtools::load_all(here::here(), quiet = TRUE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── 1. Settings ───────────────────────────────────────────────────────────────

# Distribution shown in Figure 2 (gain strip) and Figure 3 (densities).
PRIMARY_DIST <- "Log-normal"

# Number of pathogens shown in Figure 3, chosen by largest JS_CA divergence.
TOP_N_DENSITY <- 6L

# Save dimensions (inches) and resolution.
FIG_WIDTH_WIDE  <- 14
FIG_WIDTH_HALF  <- 7
FIG_HEIGHT_ROW  <- 0.55   # height per pathogen row in Figures 1 & 2
FIG_DPI         <- 300

# Arm colours (Wong colorblind-safe palette).
ARM_COLOURS <- c(
  "A" = "#0072B2",   # blue
  "B" = "#D55E00",   # vermilion
  "C" = "#009E73"    # green
)
ARM_LABELS <- c(
  "A" = "Individual-level only",
  "B" = "Summary-statistics only",
  "C" = "Federated"
)
ARM_SHAPES <- c("A" = 19, "B" = 17, "C" = 18)  # circle, triangle, diamond

OUTPUT_DIR <- here::here("results", "figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 2. Load results ───────────────────────────────────────────────────────────

rds_file <- here::here("results", "data_format_ablation.rds")
if (!file.exists(rds_file))
  stop("Results file not found: ", rds_file,
       "\nRun analysis/data_format_ablation.R first.")

stored         <- readRDS(rds_file)
comparison_tbl <- stored$comparison_tbl
ablation_fits  <- stored$ablation_fits

if (nrow(comparison_tbl) == 0L)
  stop("comparison_tbl is empty — check that data_format_ablation.R completed.")


# ── 3. Shared helpers ─────────────────────────────────────────────────────────

# Parse "9.1 (8.2, 10.3)" -> named vector c(med, lo, hi).
.parse_cri <- function(s) {
  if (is.na(s) || grepl("—|NA", s))
    return(c(med = NA_real_, lo = NA_real_, hi = NA_real_))
  nums <- suppressWarnings(
    as.numeric(unlist(regmatches(s, gregexpr("[0-9]+\\.?[0-9]*", s))))
  )
  if (length(nums) >= 3L) c(med = nums[1L], lo = nums[2L], hi = nums[3L])
  else                     c(med = NA_real_, lo = NA_real_,  hi = NA_real_)
}

# Consistent theme for all figures.
theme_ablation <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text         = element_text(face = "bold"),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      axis.title.y       = element_blank()
    )
}

# Order pathogens by federated pred_median (log-normal, descending).
.pathogen_order <- function(tbl, dist = PRIMARY_DIST) {
  tbl |>
    filter(.data$dist == !!dist) |>
    rowwise() |>
    mutate(C_med = .parse_cri(C_pred_median)[["med"]]) |>
    ungroup() |>
    arrange(C_med) |>
    pull(pathogen) |>
    unique()
}

pathogen_order <- .pathogen_order(comparison_tbl)


# ── 4. Figure 1: Three-arm forest plot ───────────────────────────────────────
#
# For each arm, parse the formatted CrI string into med/lo/hi and reshape to
# long format.  Dodge three arms vertically within each pathogen × distribution
# cell.  Log-scale x-axis with a vertical reference line at arm C's median.

.parse_arm_col <- function(tbl, col, arm_label) {
  tbl |>
    select(pathogen, dist, cri_str = all_of(col)) |>
    rowwise() |>
    mutate(parsed = list(.parse_cri(cri_str))) |>
    ungroup() |>
    mutate(
      med = vapply(parsed, `[[`, numeric(1L), "med"),
      lo  = vapply(parsed, `[[`, numeric(1L), "lo"),
      hi  = vapply(parsed, `[[`, numeric(1L), "hi"),
      arm = arm_label
    ) |>
    select(pathogen, dist, arm, med, lo, hi)
}

forest_long <- bind_rows(
  .parse_arm_col(comparison_tbl, "A_pred_median", "A"),
  .parse_arm_col(comparison_tbl, "B_pred_median", "B"),
  .parse_arm_col(comparison_tbl, "C_pred_median", "C")
) |>
  filter(!is.na(med)) |>
  mutate(
    pathogen = factor(pathogen, levels = pathogen_order),
    arm      = factor(arm, levels = c("A", "B", "C"))
  )

# Federated reference line (arm C) per facet cell.
ref_lines <- forest_long |>
  filter(arm == "C") |>
  select(pathogen, dist, ref_med = med)

dodge_height <- 0.55   # vertical spread per pathogen row

fig1 <- ggplot(forest_long,
               aes(x = med, y = pathogen,
                   colour = arm, shape = arm)) +
  # Reference line at arm C's point estimate
  geom_vline(data  = ref_lines,
             aes(xintercept = ref_med),
             colour = ARM_COLOURS[["C"]], linewidth = 0.25,
             linetype = "dashed") +
  # CrI bars
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height   = 0,
                 linewidth = 0.5,
                 position = position_dodge(width = dodge_height)) +
  # Point estimates
  geom_point(size     = 1.8,
             position = position_dodge(width = dodge_height)) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 20, 50),
    labels = c("1", "2", "5", "10", "20", "50")
  ) +
  scale_colour_manual(values = ARM_COLOURS, labels = ARM_LABELS) +
  scale_shape_manual(values  = ARM_SHAPES,  labels = ARM_LABELS) +
  facet_wrap(~ dist, nrow = 1L, scales = "free_x") +
  labs(
    x        = "Posterior predictive median (days, log scale)",
    caption  = "Horizontal bars: 95% credible intervals.  Dashed line: federated (C) point estimate."
  ) +
  theme_ablation() +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

n_pathogens <- length(pathogen_order)
fig1_height <- max(4, n_pathogens * FIG_HEIGHT_ROW + 1.5)

ggsave(file.path(OUTPUT_DIR, "fig_ablation_forest.pdf"),
       fig1, width = FIG_WIDTH_WIDE, height = fig1_height, device = "pdf")
ggsave(file.path(OUTPUT_DIR, "fig_ablation_forest.png"),
       fig1, width = FIG_WIDTH_WIDE, height = fig1_height, dpi = FIG_DPI)
message("Figure 1 saved.")


# ── 5. Figure 2: Federated-gain metrics strip ─────────────────────────────────
#
# Restricted to PRIMARY_DIST.  Three ggplot objects assembled side-by-side with
# patchwork.  Each panel shares the y-axis; only the left panel shows pathogen
# labels.

gain_df <- comparison_tbl |>
  filter(dist == PRIMARY_DIST) |>
  mutate(pathogen = factor(pathogen, levels = pathogen_order)) |>
  select(pathogen,
         ratio_CA, ratio_CB,
         JS_CA,    JS_CB,
         OVL_CA,   OVL_CB) |>
  pivot_longer(
    cols      = -pathogen,
    names_to  = c("metric", "comparison"),
    names_sep = "_",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  mutate(
    comparison = factor(comparison,
                        levels = c("CA", "CB"),
                        labels = c("vs A (individual-level)", "vs B (summary-stats)")),
    metric = factor(metric,
                    levels = c("ratio", "JS", "OVL"),
                    labels = c("Interval ratio\n(arm / federated)",
                               "Jensen-Shannon\ndivergence (bits)",
                               "Overlap\ncoefficient"))
  )

# Shared y scale
y_scale <- scale_y_discrete(drop = FALSE)

# Helper: one panel
.gain_panel <- function(df, metric_label, x_lab, ref_val = NULL,
                        show_y = TRUE) {
  d <- filter(df, metric == metric_label)
  p <- ggplot(d, aes(x = value, y = pathogen,
                     colour = comparison, shape = comparison)) +
    geom_point(size = 2.2, alpha = 0.85) +
    y_scale +
    labs(x = x_lab) +
    scale_colour_manual(
      values = c("vs A (individual-level)" = ARM_COLOURS[["A"]],
                 "vs B (summary-stats)"    = ARM_COLOURS[["B"]])
    ) +
    scale_shape_manual(
      values = c("vs A (individual-level)" = 19,
                 "vs B (summary-stats)"    = 17)
    ) +
    theme_ablation() +
    theme(legend.position = if (show_y) "none" else "bottom")

  if (!is.null(ref_val))
    p <- p + geom_vline(xintercept = ref_val,
                        linetype = "dashed", colour = "grey50",
                        linewidth = 0.4)
  if (!show_y)
    p <- p + theme(axis.text.y = element_blank())
  p
}

p2a <- .gain_panel(gain_df,
                   "Interval ratio\n(arm / federated)",
                   "Interval ratio  (> 1 = federated is tighter)",
                   ref_val  = 1,
                   show_y   = TRUE)

p2b <- .gain_panel(gain_df,
                   "Jensen-Shannon\ndivergence (bits)",
                   "JS divergence (bits)",
                   ref_val  = 0,
                   show_y   = FALSE)

p2c <- .gain_panel(gain_df,
                   "Overlap\ncoefficient",
                   "Overlap coefficient  (1 = identical)",
                   ref_val  = 1,
                   show_y   = FALSE) +
  theme(legend.position = "bottom")

fig2 <- (p2a | p2b | p2c) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "fig_ablation_gain.pdf"),
       fig2, width = FIG_WIDTH_WIDE * 0.7, height = fig1_height,
       device = "pdf")
ggsave(file.path(OUTPUT_DIR, "fig_ablation_gain.png"),
       fig2, width = FIG_WIDTH_WIDE * 0.7, height = fig1_height,
       dpi = FIG_DPI)
message("Figure 2 saved.")


# ── 6. Figure 3: Posterior predictive density overlays ───────────────────────
#
# Selects the TOP_N_DENSITY pathogens with the largest JS_CA (arm A vs
# federated) for PRIMARY_DIST.  Calls compute_predictive_cdf() for each of the
# three arms and overlays the posterior-median density with a ±95% CrI ribbon.
# Falls back gracefully if a fit is absent.

dist_code_map <- c(
  "Log-normal" = "lognormal", "Gamma" = "gamma", "Weibull" = "weibull",
  "Burr XII"   = "burr",      "Gen. gamma" = "gg"
)

# Identify the top pathogens to plot.
top_pathogens <- comparison_tbl |>
  filter(dist == PRIMARY_DIST, !is.na(JS_CA)) |>
  arrange(desc(JS_CA)) |>
  slice_head(n = TOP_N_DENSITY) |>
  pull(pathogen)

if (length(top_pathogens) == 0L) {
  message("No JS_CA values available — skipping Figure 3.")
} else {

  cdf_dist_name <- dist_code_map[[PRIMARY_DIST]]

  # Determine a common x upper bound for each pathogen using the maximum
  # pred_q95 upper CrI bound across all three arms.
  .x_upper <- function(pathogen_name) {
    row <- filter(comparison_tbl,
                  pathogen == pathogen_name & dist == PRIMARY_DIST)
    vals <- sapply(c("A_pred_q95", "B_pred_q95", "C_pred_q95"), function(col) {
      p <- tryCatch(.parse_cri(row[[col]])[["hi"]], error = function(e) NA_real_)
      p
    })
    v <- max(vals, na.rm = TRUE)
    if (!is.finite(v)) 50 else v * 1.25
  }

  # Build density data for one arm of one pathogen.
  .get_density_df <- function(fit, arm_label, pathogen_name) {
    if (is.null(fit)) return(NULL)
    x_upper <- .x_upper(pathogen_name)
    x_seq   <- seq(0.001, x_upper, length.out = 400L)
    cdf_obj <- tryCatch(
      compute_predictive_cdf(fit, cdf_dist_name, x_seq = x_seq, n_draws = 300L),
      error = function(e) NULL
    )
    if (is.null(cdf_obj)) return(NULL)

    dx <- diff(x_seq)[1L]
    .to_pdf <- function(F_vals) {
      f <- c(diff(pmax(pmin(F_vals, 1 - 1e-9), 1e-9)), 0) / dx
      f[f < 0] <- 0
      f
    }

    # Convert CDF summary columns to PDF by forward-differencing.
    data.frame(
      x        = x_seq,
      pdf_med  = .to_pdf(cdf_obj$summary$median),
      pdf_lo   = .to_pdf(cdf_obj$summary$low),     # CrI on density
      pdf_hi   = .to_pdf(cdf_obj$summary$high),
      arm      = arm_label,
      pathogen = pathogen_name
    )
  }

  density_list <- list()

  for (pg in top_pathogens) {
    arm_names   <- c("individual_only", "summary_only", "federated")
    arm_labels  <- c("A", "B", "C")
    primary_key <- names(dist_code_map[dist_code_map == cdf_dist_name &
                                         names(dist_code_map) == PRIMARY_DIST])
    # Map PRIMARY_DIST label back to DIST_CODES key
    dist_key_map <- c(
      "Log-normal" = "lognormal", "Gamma" = "gamma", "Weibull" = "weibull",
      "Burr XII"   = "burr",      "Gen. gamma" = "gengamma"
    )
    dk <- dist_key_map[[PRIMARY_DIST]]

    for (i in seq_along(arm_names)) {
      res <- ablation_fits[[pg]][[arm_names[i]]][[dk]]
      if (is.null(res) || isTRUE(res$skipped) || is.null(res$fit)) next
      df <- .get_density_df(res$fit, arm_labels[i], pg)
      if (!is.null(df)) density_list[[length(density_list) + 1L]] <- df
    }
  }

  if (length(density_list) == 0L) {
    message("No density data could be computed — skipping Figure 3.")
  } else {

    density_df <- bind_rows(density_list) |>
      mutate(
        arm      = factor(arm, levels = c("A", "B", "C")),
        pathogen = factor(pathogen, levels = top_pathogens)
      )

    # Also annotate each panel with the JS_CA and JS_CB values.
    js_labels <- comparison_tbl |>
      filter(pathogen %in% top_pathogens, dist == PRIMARY_DIST) |>
      mutate(
        label = sprintf(
          "JS[CA] == %.3f~~~~~JS[CB] == %.3f",
          round(JS_CA, 3L), round(JS_CB, 3L)
        ),
        pathogen = factor(pathogen, levels = top_pathogens)
      ) |>
      select(pathogen, label)

    # x position for label: 70th percentile of x range per pathogen
    label_x <- density_df |>
      group_by(pathogen) |>
      summarise(x_pos = quantile(x, 0.72), .groups = "drop")
    y_pos <- density_df |>
      group_by(pathogen) |>
      summarise(y_pos = max(pdf_hi, na.rm = TRUE) * 0.97, .groups = "drop")
    js_labels <- left_join(js_labels, label_x, by = "pathogen") |>
      left_join(y_pos, by = "pathogen")

    # Ribbon alpha by arm: federated slightly more prominent.
    ribbon_alpha <- c("A" = 0.20, "B" = 0.20, "C" = 0.25)
    line_lwd     <- c("A" = 0.6,  "B" = 0.6,  "C" = 0.9)

    fig3 <- ggplot(density_df,
                   aes(x = x, group = arm, fill = arm, colour = arm)) +
      geom_ribbon(aes(ymin = pdf_lo, ymax = pdf_hi),
                  alpha = 0.18, colour = NA) +
      geom_line(aes(y = pdf_med, linewidth = arm)) +
      geom_text(data  = js_labels,
                aes(x = x_pos, y = y_pos, label = label),
                inherit.aes = FALSE,
                size  = 2.6, hjust = 0, parse = TRUE, colour = "grey30") +
      scale_fill_manual(values = ARM_COLOURS,  labels = ARM_LABELS) +
      scale_colour_manual(values = ARM_COLOURS, labels = ARM_LABELS) +
      scale_linewidth_manual(values = line_lwd,  labels = ARM_LABELS) +
      facet_wrap(~ pathogen, scales = "free", ncol = 2L) +
      labs(
        x       = "Incubation period (days)",
        y       = "Density",
        caption = paste0(
          "Top ", TOP_N_DENSITY, " pathogens by JS divergence (arm A vs federated).  ",
          PRIMARY_DIST, " distribution.  ",
          "Shaded bands: 95% posterior credible intervals."
        )
      ) +
      theme_ablation() +
      theme(
        axis.title.y    = element_text(),
        legend.position = "bottom"
      ) +
      guides(
        fill      = guide_legend(override.aes = list(alpha = 0.4)),
        linewidth = guide_legend(override.aes = list(linewidth = 1))
      )

    n_rows_fig3 <- ceiling(length(top_pathogens) / 2L)
    fig3_height <- max(4, n_rows_fig3 * 3 + 1.5)

    ggsave(file.path(OUTPUT_DIR, "fig_ablation_densities.pdf"),
           fig3, width = FIG_WIDTH_HALF * 1.5, height = fig3_height,
           device = "pdf")
    ggsave(file.path(OUTPUT_DIR, "fig_ablation_densities.png"),
           fig3, width = FIG_WIDTH_HALF * 1.5, height = fig3_height,
           dpi = FIG_DPI)
    message("Figure 3 saved.")
  }
}


# ── 7. Supplementary: gain strip faceted over all distributions ───────────────
#
# Same as Figure 2 but faceted by distribution so the reader can check whether
# the federated-gain pattern is consistent across distributional assumptions.

gain_all_dists <- comparison_tbl |>
  mutate(pathogen = factor(pathogen, levels = pathogen_order)) |>
  select(pathogen, dist,
         ratio_CA, ratio_CB,
         JS_CA,    JS_CB,
         OVL_CA,   OVL_CB) |>
  pivot_longer(
    cols      = c(ratio_CA, ratio_CB, JS_CA, JS_CB, OVL_CA, OVL_CB),
    names_to  = c("metric", "comparison"),
    names_sep = "_",
    values_to = "value"
  ) |>
  filter(!is.na(value)) |>
  mutate(
    comparison = factor(comparison,
                        levels = c("CA", "CB"),
                        labels = c("vs A (individual-level)",
                                   "vs B (summary-stats)")),
    metric = factor(metric,
                    levels = c("ratio", "JS", "OVL"),
                    labels = c("Interval ratio",
                               "JS divergence (bits)",
                               "Overlap coefficient"))
  )

ref_vals <- c("Interval ratio" = 1, "JS divergence (bits)" = 0,
              "Overlap coefficient" = 1)

fig_supp <- ggplot(gain_all_dists,
                   aes(x = value, y = pathogen,
                       colour = comparison, shape = comparison)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_vline(
    data = data.frame(
      metric  = names(ref_vals),
      ref_val = unname(ref_vals)
    ),
    aes(xintercept = ref_val),
    colour   = "grey50",
    linetype = "dashed",
    linewidth = 0.35,
  ) +
  scale_colour_manual(
    values = c("vs A (individual-level)" = ARM_COLOURS[["A"]],
               "vs B (summary-stats)"    = ARM_COLOURS[["B"]])
  ) +
  scale_shape_manual(
    values = c("vs A (individual-level)" = 19,
               "vs B (summary-stats)"    = 17)
  ) +
  facet_grid(dist ~ metric, scales = "free_x") +
  labs(x = NULL) +
  theme_ablation() +
  theme(
    legend.position = "bottom",
    strip.text.y    = element_text(angle = 0, hjust = 0)
  )

ggsave(file.path(OUTPUT_DIR, "fig_ablation_gain_all_dists.pdf"),
       fig_supp,
       width  = FIG_WIDTH_WIDE * 0.7,
       height = max(5, n_pathogens * FIG_HEIGHT_ROW * 5 + 1.5),
       device = "pdf")
ggsave(file.path(OUTPUT_DIR, "fig_ablation_gain_all_dists.png"),
       fig_supp,
       width  = FIG_WIDTH_WIDE * 0.7,
       height = max(5, n_pathogens * FIG_HEIGHT_ROW * 5 + 1.5),
       dpi = FIG_DPI)
message("Supplementary gain figure (all distributions) saved.")


# ── 8. Figure 4 (combined): Multi-panel federated gain summary ───────────────
#
# A single three-column figure suitable for the main paper body.
#
#   Panel A (wide) — Dumbbell forest plot (PRIMARY_DIST only)
#     A grey segment connects the arm A and arm B point estimates; its length
#     encodes how much the two single-format arms disagree.  Three
#     point + CrI layers (A, B, C) are drawn on top.  Arm C is rendered more
#     prominently (larger point, heavier error bar) because it is the focus.
#     A log-scaled x-axis accommodates the range across pathogens.
#
#   Panel B (narrow) — Interval ratio  (arm_width / C_width)
#     One dot per arm × pathogen pair (CA in blue, CB in orange).
#     Reference line at 1: points to the right mean the single-format arm is
#     wider (less precise) than the federated model.
#
#   Panel C (narrow) — Jensen-Shannon divergence  (bits)
#     JS_CA and JS_CB on the same scale.  Larger values signal a bigger
#     distributional shift between that arm and the federated result.
#     Reference line at 0 (identical distributions).
#
#   Panels A/B/C share the pathogen y-axis via patchwork alignment.
#   Panel labels A/B/C are added automatically by plot_annotation().
#   The combined legend is collected at the bottom.
#
#   Output: results/figures/fig_ablation_combined.{pdf,png}

# ── 8a. Best-fitting distribution per pathogen (from main analysis) ──────────
#
# Load model weights pre-computed from main_results.rds via
# compute_pathogen_model_bayes_factors().  The best distribution is the
# first entry in each pathogen's weight vector (highest pseudo-Bayes factor
# weight).  This is the same criterion used in plot_main_figure() and
# make_supplementary_figures.R, ensuring consistency across the paper.
#
# Falls back to PRIMARY_DIST for pathogens absent from model_weights or
# whose best distribution was not fitted in the ablation.

# Map internal Stan dist names -> display labels used in comparison_tbl$dist
MAIN_DIST_TO_DISPLAY <- c(
  lognormal = "Log-normal",
  gamma     = "Gamma",
  weibull   = "Weibull",
  burr      = "Burr XII",
  gengamma  = "Gen. gamma"
)

model_weights_file <- here::here("results", "model_weights.rds")
if (!file.exists(model_weights_file))
  stop("model_weights.rds not found at: ", model_weights_file,
       "\nPre-compute via: mw <- compute_pathogen_model_bayes_factors(all_results); ",
       "saveRDS(mw, 'results/model_weights.rds')")

main_model_weights <- readRDS(model_weights_file)

# For each ablation pathogen, extract the highest-weight distribution from
# the main analysis.
ablation_pathogens <- unique(comparison_tbl$pathogen)
best_internal <- vapply(ablation_pathogens, function(p) {
  w <- main_model_weights[[p]]
  if (is.null(w) || length(w) == 0L) return(NA_character_)
  names(w)[1L]   # best = highest pseudo-BF weight
}, character(1L))

best_dist_tbl <- data.frame(
  pathogen  = ablation_pathogens,
  best_dist = MAIN_DIST_TO_DISPLAY[best_internal],
  stringsAsFactors = FALSE
)

# Fall back to PRIMARY_DIST where the main-analysis best dist is unavailable
# in the ablation results (e.g. Gen. gamma was not fitted in arm B/C).
available_dists <- unique(comparison_tbl$dist)
needs_fallback  <- is.na(best_dist_tbl$best_dist) |
                   !(best_dist_tbl$best_dist %in% available_dists)
if (any(needs_fallback)) {
  message("Note: falling back to PRIMARY_DIST for: ",
          paste(best_dist_tbl$pathogen[needs_fallback], collapse = ", "))
  best_dist_tbl$best_dist[needs_fallback] <- PRIMARY_DIST
}

# ── 8b. Wide data using best distribution per pathogen ───────────────────────

wide_best <- comparison_tbl |>
  inner_join(best_dist_tbl, by = c("pathogen", "dist" = "best_dist")) |>
  mutate(pathogen = factor(pathogen, levels = pathogen_order)) |>
  rowwise() |>
  mutate(
    # Median (P50) posterior predictive estimates + 95 % CrI bounds
    A_med    = .parse_cri(A_pred_median)[["med"]],
    A_lo     = .parse_cri(A_pred_median)[["lo"]],
    A_hi     = .parse_cri(A_pred_median)[["hi"]],
    B_med    = .parse_cri(B_pred_median)[["med"]],
    B_lo     = .parse_cri(B_pred_median)[["lo"]],
    B_hi     = .parse_cri(B_pred_median)[["hi"]],
    C_med    = .parse_cri(C_pred_median)[["med"]],
    C_lo     = .parse_cri(C_pred_median)[["lo"]],
    C_hi     = .parse_cri(C_pred_median)[["hi"]],
    # 95th percentile (P95) posterior predictive estimates + 95 % CrI bounds
    A_q95    = .parse_cri(A_pred_q95)[["med"]],
    A_q95_lo = .parse_cri(A_pred_q95)[["lo"]],
    A_q95_hi = .parse_cri(A_pred_q95)[["hi"]],
    B_q95    = .parse_cri(B_pred_q95)[["med"]],
    B_q95_lo = .parse_cri(B_pred_q95)[["lo"]],
    B_q95_hi = .parse_cri(B_pred_q95)[["hi"]],
    C_q95    = .parse_cri(C_pred_q95)[["med"]],
    C_q95_lo = .parse_cri(C_pred_q95)[["lo"]],
    C_q95_hi = .parse_cri(C_pred_q95)[["hi"]]
  ) |>
  ungroup()

# y-axis labels: "Pathogen (Distribution)".
dist_label_lookup <- setNames(
  paste0(wide_best$pathogen, "\n(", wide_best$dist, ")"),
  as.character(wide_best$pathogen)
)

# ── 8c. Forest data — P50 and P95 facets ─────────────────────────────────────

forest_p50 <- wide_best |>
  select(pathogen, A_med, A_lo, A_hi, B_med, B_lo, B_hi, C_med, C_lo, C_hi) |>
  pivot_longer(-pathogen, names_to = c("arm", ".value"),
               names_pattern = "^(.)_(.*)") |>
  filter(!is.na(med)) |>
  mutate(arm = factor(arm, levels = c("A", "B", "C")), metric = "Median (P50)")

forest_p95 <- wide_best |>
  select(pathogen,
         A_med = A_q95,    A_lo = A_q95_lo, A_hi = A_q95_hi,
         B_med = B_q95,    B_lo = B_q95_lo, B_hi = B_q95_hi,
         C_med = C_q95,    C_lo = C_q95_lo, C_hi = C_q95_hi) |>
  pivot_longer(-pathogen, names_to = c("arm", ".value"),
               names_pattern = "^(.)_(.*)") |>
  filter(!is.na(med)) |>
  mutate(arm = factor(arm, levels = c("A", "B", "C")), metric = "95th percentile (P95)")

# Dumbbell backbone segments for both facets (A-to-B range, un-dodged).
dumbbell_segs <- bind_rows(
  wide_best |> filter(!is.na(A_med), !is.na(B_med)) |>
    transmute(pathogen, x = A_med, xend = B_med, metric = "Median (P50)"),
  wide_best |> filter(!is.na(A_q95), !is.na(B_q95)) |>
    transmute(pathogen, x = A_q95, xend = B_q95, metric = "95th percentile (P95)")
) |>
  mutate(metric = factor(metric, levels = c("Median (P50)", "95th percentile (P95)")))

# Dodging width and size scales (used by panels A and C).
COMB_DODGE  <- 0.5
ARM_SIZES_C <- c("A" = 1.8, "B" = 1.8, "C" = 2.8)
ARM_EBW_C   <- c("A" = 0.35, "B" = 0.35, "C" = 0.75)

# ── 8d. Extract mu0 CrI widths from Stan fits ─────────────────────────────────
#
# mu0 is the population-level location parameter (distribution's log scale).
# Its 95% CrI width captures pure estimation uncertainty, independent of τ.
# Ratio arm / C > 1 means arm C has tighter knowledge of the population mean.

.extract_mu0_cri <- function(fit_slot) {
  empty <- c(med = NA_real_, lo = NA_real_, hi = NA_real_, width = NA_real_)
  if (is.null(fit_slot) || isTRUE(fit_slot$skipped) || is.null(fit_slot$fit))
    return(empty)
  draws <- tryCatch(
    rstan::extract(fit_slot$fit, pars = "mu0")$mu0,
    error = function(e) NULL
  )
  if (is.null(draws) || length(draws) == 0L) return(empty)
  q <- quantile(draws, c(0.025, 0.5, 0.975), na.rm = TRUE)
  c(med = q[2L], lo = q[1L], hi = q[3L], width = q[3L] - q[1L])
}

DISPLAY_TO_INTERNAL <- setNames(names(MAIN_DIST_TO_DISPLAY), MAIN_DIST_TO_DISPLAY)
ARM_FIT_KEYS        <- c("A" = "individual_only", "B" = "summary_only",
                         "C" = "federated")

mu0_tbl <- do.call(rbind, lapply(best_dist_tbl$pathogen, function(p) {
  bd_int <- DISPLAY_TO_INTERNAL[[ best_dist_tbl$best_dist[best_dist_tbl$pathogen == p] ]]
  if (is.na(bd_int)) return(NULL)
  do.call(rbind, lapply(c("A", "B", "C"), function(arm) {
    slot <- ablation_fits[[p]][[ ARM_FIT_KEYS[[arm]] ]][[ bd_int ]]
    mt   <- .extract_mu0_cri(slot)
    data.frame(pathogen = p, arm = arm,
               mu0_med = mt["med"], mu0_lo = mt["lo"],
               mu0_hi  = mt["hi"], mu0_width = mt["width"],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
})) |>
  filter(!is.na(mu0_width)) |>
  mutate(arm      = factor(arm, levels = c("A", "B", "C")),
         pathogen = factor(pathogen, levels = levels(wide_best$pathogen)))

# mu0 CrI width ratio: arm / C  (> 1 → arm C is more precise about μ).
mu0_ratio <- mu0_tbl |>
  select(pathogen, arm, mu0_width) |>
  pivot_wider(names_from = arm, values_from = mu0_width,
              names_prefix = "w_") |>
  filter(!is.na(w_C)) |>
  pivot_longer(cols = c(w_A, w_B), names_to = "arm_key", values_to = "ratio") |>
  filter(!is.na(ratio)) |>
  mutate(
    ratio      = ratio / w_C,
    comparison = factor(arm_key,
                        levels = c("w_A", "w_B"),
                        labels = c("Individual-level only",
                                   "Summary-statistics only"))
  )

# ── 8e. Parse τ from comparison_tbl ──────────────────────────────────────────
#
# τ is the between-study heterogeneity SD.  Only available when n_datasets ≥ 5
# in that arm (stored as "— (n<5)" otherwise); missing rows are dropped.

tau_long <- wide_best |>
  rowwise() |>
  mutate(
    A_tau_med = .parse_cri(A_tau)[["med"]],
    A_tau_lo  = .parse_cri(A_tau)[["lo"]],
    A_tau_hi  = .parse_cri(A_tau)[["hi"]],
    B_tau_med = .parse_cri(B_tau)[["med"]],
    B_tau_lo  = .parse_cri(B_tau)[["lo"]],
    B_tau_hi  = .parse_cri(B_tau)[["hi"]],
    C_tau_med = .parse_cri(C_tau)[["med"]],
    C_tau_lo  = .parse_cri(C_tau)[["lo"]],
    C_tau_hi  = .parse_cri(C_tau)[["hi"]]
  ) |>
  ungroup() |>
  select(pathogen, matches("^[ABC]_tau_(med|lo|hi)$")) |>
  pivot_longer(-pathogen,
               names_to      = c("arm", ".value"),
               names_pattern = "^(.)_tau_(.*)$") |>
  filter(!is.na(med)) |>
  mutate(arm = factor(arm, levels = c("A", "B", "C")))

# ── 8f. Unified arm labels for legend merging ─────────────────────────────────
#
# Convert the "A"/"B"/"C" factor to full label text in all forest/tau/mu0 data.
# All panels then use the same colour/shape scale keyed by label text, so
# patchwork's guides = "collect" produces a single merged legend.

ARM_COLOURS_FULL <- setNames(ARM_COLOURS, ARM_LABELS[names(ARM_COLOURS)])
ARM_SHAPES_FULL  <- setNames(ARM_SHAPES,  ARM_LABELS[names(ARM_SHAPES)])
ARM_SIZES_C_FULL <- setNames(ARM_SIZES_C, ARM_LABELS[names(ARM_SIZES_C)])
ARM_EBW_C_FULL   <- setNames(ARM_EBW_C,   ARM_LABELS[names(ARM_EBW_C)])

COMP_COLOURS <- ARM_COLOURS_FULL[c("Individual-level only",
                                    "Summary-statistics only")]
COMP_SHAPES  <- ARM_SHAPES_FULL[ c("Individual-level only",
                                    "Summary-statistics only")]

.relabel_arm <- function(arm_fac) {
  factor(ARM_LABELS[as.character(arm_fac)], levels = unname(ARM_LABELS))
}

forest_both <- bind_rows(forest_p50, forest_p95) |>
  mutate(metric = factor(metric,
                          levels = c("Median (P50)", "95th percentile (P95)")),
         arm    = .relabel_arm(arm))

tau_long <- tau_long |> mutate(arm = .relabel_arm(arm))
mu0_tbl  <- mu0_tbl  |> mutate(arm = .relabel_arm(arm))

# ── 8g. Build all four panels ─────────────────────────────────────────────────

# Shared y-axis theme for the narrow strip panels (labels suppressed).
y_shared <- list(
  scale_y_discrete(drop = FALSE),
  theme_ablation(),
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
)

# Panel A: faceted dumbbell forest (P50 | P95).
pA <- ggplot(forest_both,
             aes(x = med, y = pathogen, colour = arm, shape = arm)) +
  geom_segment(
    data        = dumbbell_segs,
    aes(x = x, xend = xend, y = pathogen, yend = pathogen),
    colour      = "grey78", linewidth = 0.9, lineend = "round",
    inherit.aes = FALSE
  ) +
  geom_errorbarh(
    aes(xmin = lo, xmax = hi, linewidth = arm),
    height   = 0,
    position = position_dodge(width = COMB_DODGE)
  ) +
  geom_point(aes(size = arm), position = position_dodge(width = COMB_DODGE)) +
  facet_wrap(~ metric, nrow = 1L, scales = "free_x") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100),
                labels = c("1", "2", "5", "10", "20", "50", "100")) +
  scale_y_discrete(labels = dist_label_lookup) +
  scale_colour_manual(values = ARM_COLOURS_FULL) +
  scale_shape_manual( values = ARM_SHAPES_FULL) +
  scale_size_manual(  values = ARM_SIZES_C_FULL) +
  scale_linewidth_manual(values = ARM_EBW_C_FULL) +
  labs(x        = "Days (log scale)",
       subtitle = "Best-fitting distribution per pathogen (main analysis model weights)") +
  theme_ablation() +
  guides(colour    = guide_legend(override.aes = list(size = 2.5)),
         size      = "none",
         linewidth = "none",
         shape     = "none")

# Panel B: μ₀ CrI width ratio — pure information gain.
# Values > 1 indicate arm C has tighter posterior for the population mean.
pB <- ggplot(mu0_ratio,
             aes(x = ratio, y = pathogen,
                 colour = comparison, shape = comparison)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = COMP_COLOURS) +
  scale_shape_manual( values = COMP_SHAPES) +
  labs(x = expression(mu[0]~"CrI ratio (arm / C)"),
       subtitle = "Information gain") +
  y_shared

# Panel C: τ per arm — between-study heterogeneity.
# Arm C may detect larger τ when cross-data-type contrast reveals more
# between-study variation; its τ posterior is also better estimated
# (narrower CrI) due to more studies.
pC <- ggplot(tau_long,
             aes(x = med, y = pathogen, colour = arm, shape = arm)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height    = 0, linewidth = 0.35,
                 position  = position_dodge(width = COMB_DODGE)) +
  geom_point(size = 2.0, position = position_dodge(width = COMB_DODGE)) +
  scale_colour_manual(values = ARM_COLOURS_FULL) +
  scale_shape_manual( values = ARM_SHAPES_FULL) +
  labs(x        = expression(tau~"(heterogeneity SD)"),
       subtitle = "Between-study heterogeneity") +
  y_shared

# Panel D: predictive CrI ratio — the combined (confounded) signal.
pred_ratio <- wide_best |>
  select(pathogen, ratio_CA, ratio_CB) |>
  pivot_longer(cols      = c(ratio_CA, ratio_CB),
               names_to  = "comparison",
               values_to = "ratio") |>
  filter(!is.na(ratio)) |>
  mutate(comparison = factor(comparison,
                              levels = c("ratio_CA", "ratio_CB"),
                              labels = c("Individual-level only",
                                         "Summary-statistics only")))

pD <- ggplot(pred_ratio,
             aes(x = ratio, y = pathogen,
                 colour = comparison, shape = comparison)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = COMP_COLOURS) +
  scale_shape_manual( values = COMP_SHAPES) +
  labs(x        = "Predictive CrI ratio (arm / C)",
       subtitle = "Combined effect") +
  y_shared

# ── 8h. Assemble and save ─────────────────────────────────────────────────────

fig4 <- (pA | pB | pC | pD) +
  plot_layout(widths = c(2, 1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag        = element_text(face = "bold", size = 10))

fig4_height <- max(4, n_pathogens * FIG_HEIGHT_ROW + 1.8)

ggsave(file.path(OUTPUT_DIR, "fig_ablation_combined.pdf"),
       fig4, width = FIG_WIDTH_WIDE * 1.35, height = fig4_height, device = "pdf")
ggsave(file.path(OUTPUT_DIR, "fig_ablation_combined.png"),
       fig4, width = FIG_WIDTH_WIDE * 1.35, height = fig4_height, dpi = FIG_DPI)
message("Figure 4 (combined) saved.")


message("\nAll figures written to: ", OUTPUT_DIR)
