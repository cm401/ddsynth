# analysis/make_supplementary_figures.R
# =============================================================================
# Generates per-pathogen supplementary figures for the incubation period paper.
#
# Each figure has up to five panels (not all panels appear for every pathogen):
#
#   A: Posterior predictive CDF comparison across all fitted distributions
#      (using the filtered dataset)
#   B: All-data vs filtered-data CDF for the best-fit distribution
#      (only shown when filtering removed at least one dataset)
#   C: Overall vs subgroup CDFs for the best-fit distribution
#      (only shown when subgroup analyses exist)
#   D: Forest-plot of raw summary statistics (median ± range/IQR, mean ± SD)
#      (only shown when at least one type-A/B/C dataset exists)
#   F: Frequency table data — bar charts for exact tables and horizontal-
#      segment plots for interval-censored observations
#      (only shown when at least one type-D/E dataset exists)
#
# Usage:
#   Rscript analysis/make_supplementary_figures.R
#   -- or --
#   source("analysis/make_supplementary_figures.R")
#
# Outputs are written to results/supplementary_figures/<pathogen>.pdf
# =============================================================================

devtools::load_all(here::here(), quiet = TRUE)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── Grab unexported package objects needed below ──────────────────────────────

.DIST_COLORS       <- ddsynth:::.DIST_COLORS
.DIST_LABELS       <- ddsynth:::.DIST_LABELS
.DIST_CDF_NAME     <- ddsynth:::.DIST_CDF_NAME
.SUBGROUP_PALETTE  <- ddsynth:::.SUBGROUP_PALETTE
.SUBGROUP_LABELS   <- ddsynth:::.SUBGROUP_LABELS
.PATHOGEN_LABELS   <- ddsynth:::.PATHOGEN_LABELS
.fit_has_converged <- ddsynth:::.fit_has_converged

# ── Colour maps for data panels ───────────────────────────────────────────────

.SUMMARY_TYPE_COLOURS <- c(
  "Median + Range" = "#2166AC",
  "Median + IQR"   = "#1A9641",
  "Mean + SD"      = "#E08214"
)

# Shape: circle (19) for median-based; square (15) for mean-based
.SUMMARY_SHAPES <- c(
  "Median + Range" = 19,
  "Median + IQR"   = 19,
  "Mean + SD"      = 15
)

# Linetype: solid for range, dashed for IQR, dotted for SD
.SUMMARY_LINETYPES <- c(
  "Median + Range" = "solid",
  "Median + IQR"   = "dashed",
  "Mean + SD"      = "dotted"
)


# =============================================================================
# Utility helpers
# =============================================================================

# Return the dataset list from the first non-skipped fit in an analysis slot.
.get_datasets <- function(analysis_slot) {
  for (r in analysis_slot) {
    if (!is.null(r) && !isTRUE(r$skipped) && !is.null(r$datasets))
      return(r$datasets)
  }
  NULL
}

# Detect the summary-stat type of a single dataset entry.
.detect_type <- function(d) {
  if (!is.null(d$freq_lower) && !is.null(d$freq_upper)) return("E")
  if (!is.null(d$freq_value) && !is.null(d$freq_count)) return("D")
  if (!is.null(d$mean)   && !is.null(d$sd))  return("C")
  if (!is.null(d$median) && !is.null(d$Q1))  return("B")
  if (!is.null(d$median) && !is.null(d$min)) return("A")
  "unknown"
}

# Short citation label: "Surname (year)" with truncation fallback.
.short_cite <- function(source) {
  if (is.null(source) || !nzchar(source)) return("Unknown")
  m <- regmatches(source, regexpr("^[^,]+\\(\\d{4}\\)", source))
  if (length(m) == 1L && nzchar(m)) return(trimws(m))
  if (nchar(source) > 40L) paste0(substr(source, 1L, 37L), "...") else source
}

# Derive x_max for CDF plots: nearest 5-day boundary above the 99th percentile.
.auto_xmax <- function(fit, cdf_dname, fallback = 60L) {
  cdf_coarse <- tryCatch(
    compute_predictive_cdf(fit, cdf_dname,
                           x_seq   = seq(0, 200, length.out = 100L),
                           n_draws = 40L),
    error = function(e) NULL
  )
  if (is.null(cdf_coarse)) return(fallback)
  idx <- which(cdf_coarse$summary$median >= 0.99)[1L]
  if (is.na(idx)) return(fallback)
  x99 <- cdf_coarse$summary$x[idx]
  max(10L, ceiling(x99 / 5) * 5L + 5L)
}

# Shared ggplot2 theme for CDF panels.
.cdf_theme <- function(base_size = 9) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(colour = "grey92"),
      legend.position   = "right",
      legend.text       = ggplot2::element_text(size = base_size * 0.80),
      legend.key.width  = ggplot2::unit(0.9, "cm"),
      axis.title        = ggplot2::element_text(size = base_size * 0.90),
      axis.text         = ggplot2::element_text(size = base_size * 0.80),
      plot.title        = ggplot2::element_text(
        size = base_size, face = "bold",
        margin = ggplot2::margin(b = 3)
      )
    )
}

# Shared CDF y-axis scale (0–1 with P50 and P95 marked).
.cdf_y_scale <- function() {
  ggplot2::scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1),
    labels = c("0", ".25", ".5", ".75", ".95", "1"),
    limits = c(0, 1),
    expand = c(0, 0)
  )
}

# Build a filled-polygon data.frame for an empirical step CDF.
# Returns x/y columns suitable for geom_polygon (ymin = 0, ymax = step).
.ecdf_polygon <- function(ecdf_x, ecdf_y, x_start = 0, x_end = NULL) {
  n <- length(ecdf_x)
  if (n == 0L) return(data.frame(x = numeric(0), y = numeric(0)))
  if (is.null(x_end)) x_end <- ecdf_x[n] * 1.05

  # Build upper staircase boundary
  upper_x <- c(x_start, ecdf_x[1L])
  upper_y <- c(0, 0)
  for (i in seq_len(n)) {
    next_x  <- if (i < n) ecdf_x[i + 1L] else x_end
    upper_x <- c(upper_x, ecdf_x[i], next_x)
    upper_y <- c(upper_y, ecdf_y[i], ecdf_y[i])
  }
  # Close the polygon along the x-axis
  data.frame(x = c(upper_x, x_end, x_start),
             y = c(upper_y, 0,     0))
}

# Build a filled-polygon data.frame for the uncertainty band between two
# step-function ECDFs (e.g. lower-bound and upper-bound ECDFs for interval-
# censored data).  The polygon traces the top staircase forward, then the
# bottom staircase in reverse, forming a closed band.
#
# top_x/top_y : staircase for the top boundary   (higher CDF values)
# bot_x/bot_y : staircase for the bottom boundary (lower  CDF values)
.ecdf_band_polygon <- function(top_x, top_y, bot_x, bot_y,
                                x_start = 0, x_end = NULL) {
  if (is.null(x_end)) x_end <- max(c(top_x, bot_x)) * 1.05

  make_stair <- function(x, y) {
    n  <- length(x)
    px <- c(x_start, x[1L])
    py <- c(0, 0)
    for (i in seq_len(n)) {
      nx <- if (i < n) x[i + 1L] else x_end
      px <- c(px, x[i], nx)
      py <- c(py, y[i], y[i])
    }
    list(x = px, y = py)
  }

  top <- make_stair(top_x, top_y)
  bot <- make_stair(bot_x, bot_y)

  # Forward along top boundary, backward along bottom → closed polygon
  data.frame(x = c(top$x, rev(bot$x)),
             y = c(top$y, rev(bot$y)))
}


# =============================================================================
# Panel A: Distribution comparison (all converged fits, filtered data)
# =============================================================================

.make_panel_A <- function(pathogen_results, x_max, best_dist = NULL,
                          n_draws = 200L, base_size = 9) {

  dist_results <- pathogen_results[["filtered"]]
  if (is.null(dist_results)) return(NULL)

  x_seq   <- seq(0, x_max, length.out = 300L)
  cdf_dfs <- list()

  for (dist_name in names(dist_results)) {
    r <- dist_results[[dist_name]]
    if (is.null(r) || isTRUE(r$skipped) || is.null(r$fit)) next
    if (!.fit_has_converged(r$fit))                          next
    cdf_dname <- .DIST_CDF_NAME[[dist_name]]
    cdf <- tryCatch(
      compute_predictive_cdf(r$fit, cdf_dname, x_seq = x_seq, n_draws = n_draws),
      error = function(e) NULL
    )
    if (!is.null(cdf))
      cdf_dfs[[dist_name]] <- dplyr::mutate(cdf$summary, dist = dist_name)
  }

  if (length(cdf_dfs) == 0) return(NULL)

  all_df <- dplyr::bind_rows(cdf_dfs)
  dist_order  <- c("lognormal", "gamma", "weibull", "burr", "gengamma")
  present     <- intersect(dist_order, unique(all_df$dist))
  col_vals    <- setNames(.DIST_COLORS[present], .DIST_LABELS[present])
  all_df$dist <- factor(all_df$dist, levels = present,
                        labels = .DIST_LABELS[present])

  p <- ggplot2::ggplot(all_df, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low, ymax = high, fill = dist),
      alpha = 0.10, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median, colour = dist),
      linewidth = 0.75
    ) +
    ggplot2::scale_colour_manual(values = col_vals, name = "Distribution") +
    ggplot2::scale_fill_manual(  values = col_vals, name = "Distribution") +
    .cdf_y_scale() +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0.01, 0)) +
    ggplot2::labs(
      title = "A. Distribution comparison (filtered data)",
      x     = "Days",
      y     = "Cumulative probability"
    ) +
    .cdf_theme(base_size)

  # Annotate the best-fitting distribution name in the lower-right corner
  if (!is.null(best_dist) && best_dist %in% names(.DIST_LABELS)) {
    p <- p + ggplot2::annotate(
      "text",
      x        = x_max * 0.97,
      y        = 0.02,
      label    = .DIST_LABELS[[best_dist]],
      colour   = .DIST_COLORS[[best_dist]],
      size     = 2.5,
      hjust    = 1,
      vjust    = 0,
      fontface = "italic"
    )
  }

  p
}


# =============================================================================
# Panel B: All data vs filtered data (best-fit distribution only)
# =============================================================================

.make_panel_B <- function(pathogen_results, best_dist, x_max,
                          n_draws = 200L, base_size = 9) {

  res_all  <- pathogen_results[["all"]][[best_dist]]
  res_filt <- pathogen_results[["filtered"]][[best_dist]]

  if (is.null(res_all)  || isTRUE(res_all$skipped)  || is.null(res_all$fit))  return(NULL)
  if (is.null(res_filt) || isTRUE(res_filt$skipped) || is.null(res_filt$fit)) return(NULL)

  # Skip if the dataset lists are identical (no filtering occurred)
  if (setequal(names(res_all$datasets), names(res_filt$datasets))) return(NULL)

  x_seq     <- seq(0, x_max, length.out = 300L)
  cdf_dname <- .DIST_CDF_NAME[[best_dist]]

  cdf_all  <- tryCatch(
    compute_predictive_cdf(res_all$fit,  cdf_dname, x_seq = x_seq, n_draws = n_draws),
    error = function(e) NULL
  )
  cdf_filt <- tryCatch(
    compute_predictive_cdf(res_filt$fit, cdf_dname, x_seq = x_seq, n_draws = n_draws),
    error = function(e) NULL
  )
  if (is.null(cdf_all) || is.null(cdf_filt)) return(NULL)

  n_removed <- length(res_all$datasets) - length(res_filt$datasets)
  subtitle  <- sprintf("%d dataset%s removed by filtering",
                       n_removed, if (n_removed == 1L) "" else "s")

  all_df <- dplyr::bind_rows(
    dplyr::mutate(cdf_all$summary,  analysis = "All data"),
    dplyr::mutate(cdf_filt$summary, analysis = "Filtered")
  )
  col_vals <- c("All data" = "grey50", "Filtered" = .DIST_COLORS[[best_dist]])

  ggplot2::ggplot(all_df, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = low, ymax = high, fill = analysis),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median, colour = analysis),
      linewidth = 0.75
    ) +
    ggplot2::scale_colour_manual(values = col_vals, name = NULL) +
    ggplot2::scale_fill_manual(  values = col_vals, name = NULL) +
    .cdf_y_scale() +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0.01, 0)) +
    ggplot2::labs(
      title    = "B. All data vs. filtered data",
      subtitle = subtitle,
      x        = "Days",
      y        = "Cumulative probability"
    ) +
    .cdf_theme(base_size)
}


# =============================================================================
# Panel C: Subgroup comparison (best-fit distribution)
# =============================================================================

.make_panel_C <- function(pathogen_results, best_dist, x_max,
                          n_draws = 200L, base_size = 9) {

  sg_keys <- setdiff(names(pathogen_results), c("all", "filtered"))
  if (length(sg_keys) == 0L) return(NULL)

  res_filt <- pathogen_results[["filtered"]][[best_dist]]
  if (is.null(res_filt) || isTRUE(res_filt$skipped) || is.null(res_filt$fit)) return(NULL)

  x_seq     <- seq(0, x_max, length.out = 300L)
  cdf_dname <- .DIST_CDF_NAME[[best_dist]]

  cdf_ov <- tryCatch(
    compute_predictive_cdf(res_filt$fit, cdf_dname, x_seq = x_seq, n_draws = n_draws),
    error = function(e) NULL
  )
  if (is.null(cdf_ov)) return(NULL)

  cdf_list <- list(dplyr::mutate(cdf_ov$summary, analysis = "Overall"))

  for (sg in sg_keys) {
    r <- pathogen_results[[sg]][[best_dist]]
    if (is.null(r) || isTRUE(r$skipped) || is.null(r$fit)) next
    sg_label <- if (sg %in% names(.SUBGROUP_LABELS)) .SUBGROUP_LABELS[[sg]] else sg
    cdf <- tryCatch(
      compute_predictive_cdf(r$fit, cdf_dname, x_seq = x_seq,
                             n_draws = ceiling(n_draws / 2L)),
      error = function(e) NULL
    )
    if (!is.null(cdf))
      cdf_list[[sg]] <- dplyr::mutate(cdf$summary, analysis = sg_label)
  }

  if (length(cdf_list) <= 1L) return(NULL)  # no subgroups could be computed

  all_df    <- dplyr::bind_rows(cdf_list)
  sg_labels <- setdiff(unique(all_df$analysis), "Overall")
  n_sg      <- length(sg_labels)

  col_vals <- c(
    "Overall" = .DIST_COLORS[[best_dist]],
    setNames(.SUBGROUP_PALETTE[seq_len(n_sg)], sg_labels)
  )
  lwd_vals <- c("Overall" = 0.9, setNames(rep(0.55, n_sg), sg_labels))

  ggplot2::ggplot(all_df, ggplot2::aes(x = x)) +
    # 95% credible ribbon for the overall estimate only
    ggplot2::geom_ribbon(
      data = dplyr::filter(all_df, analysis == "Overall"),
      ggplot2::aes(ymin = low, ymax = high),
      fill   = .DIST_COLORS[[best_dist]],
      alpha  = 0.15,
      colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median, colour = analysis, linewidth = analysis)
    ) +
    ggplot2::scale_colour_manual(
      values = col_vals, name = NULL,
      breaks = sg_labels   # "Overall" identified by ribbon; keep legend compact
    ) +
    ggplot2::scale_linewidth_manual(values = lwd_vals, guide = "none") +
    .cdf_y_scale() +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0.01, 0)) +
    ggplot2::labs(
      title = "C. Subgroup comparison",
      x     = "Days",
      y     = "Cumulative probability"
    ) +
    .cdf_theme(base_size)
}


# =============================================================================
# Panel D: Summary statistics (forest-plot style)
# =============================================================================

.make_panel_D <- function(datasets_all, datasets_filtered, base_size = 9) {

  filtered_names <- names(datasets_filtered)
  rows <- list()

  for (nm in names(datasets_all)) {
    d    <- datasets_all[[nm]]
    type <- .detect_type(d)
    if (!type %in% c("A", "B", "C")) next

    included <- nm %in% filtered_names
    n_val    <- if (!is.null(d$n)) as.integer(d$n) else NA_integer_
    cite     <- .short_cite(d$source)
    subgrp   <- if (!is.null(d$subgroup) && nzchar(d$subgroup)) d$subgroup else NA_character_
    country  <- if (!is.null(d$country)  && nzchar(d$country))  d$country  else NA_character_

    lbl <- cite
    if (!is.na(n_val))   lbl <- paste0(lbl, " (n=", n_val, ")")
    if (!is.na(subgrp))  lbl <- paste0(lbl, " [", subgrp, "]")
    if (!is.na(country)) lbl <- paste0(lbl, " — ", country)

    row <- switch(type,
      A = data.frame(
        label    = lbl,
        central  = d$median,
        lower    = d$min,
        upper    = d$max,
        type     = "Median + Range",
        included = included,
        stringsAsFactors = FALSE
      ),
      B = data.frame(
        label    = lbl,
        central  = d$median,
        lower    = d$Q1,
        upper    = d$Q3,
        type     = "Median + IQR",
        included = included,
        stringsAsFactors = FALSE
      ),
      C = data.frame(
        label    = lbl,
        central  = d$mean,
        lower    = d$mean - d$sd,
        upper    = d$mean + d$sd,
        type     = "Mean + SD",
        included = included,
        stringsAsFactors = FALSE
      )
    )
    rows[[nm]] <- row
  }

  if (length(rows) == 0L) return(NULL)

  df <- dplyr::bind_rows(rows)
  # Reverse order so first dataset appears at top of the plot
  df$label    <- factor(df$label, levels = rev(unique(df$label)))
  df$included <- factor(ifelse(df$included, "Included", "Excluded"),
                        levels = c("Included", "Excluded"))

  ggplot2::ggplot(df, ggplot2::aes(y = label, colour = type, alpha = included)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = lower, xend = upper, yend = label, linetype = type),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = central, shape = type),
      size = 2.5
    ) +
    ggplot2::scale_colour_manual(
      values = .SUMMARY_TYPE_COLOURS,
      name   = "Summary type"
    ) +
    ggplot2::scale_shape_manual(
      values = .SUMMARY_SHAPES,
      name   = "Summary type"
    ) +
    ggplot2::scale_linetype_manual(
      values = .SUMMARY_LINETYPES,
      name   = "Summary type"
    ) +
    ggplot2::scale_alpha_manual(
      values = c("Included" = 1.0, "Excluded" = 0.25),
      name   = NULL,
      labels = c("Included" = "In filtered analysis",
                 "Excluded" = "Excluded by filter")
    ) +
    ggplot2::labs(
      title = "D. Raw summary statistics by dataset",
      x     = "Days",
      y     = NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = base_size * 0.75),
      legend.position    = "right",
      plot.title         = ggplot2::element_text(
        size = base_size, face = "bold",
        margin = ggplot2::margin(b = 3)
      )
    )
}


# =============================================================================
# Panel F: Frequency table data — empirical CDF + fitted CDF overlay
# =============================================================================
#
# Each freq-table dataset (type D = exact values; type E = interval-censored)
# gets its own sub-panel showing:
#   • Green filled polygon  — empirical step CDF
#   • Black step line       — empirical step CDF outline
#   • Coloured ribbon + line — posterior predictive CDF (best-fit distribution)
#
# For interval-censored data (type E) the empirical CDF is computed at the
# upper bound of each exposure interval (conservative / "at most by" CDF).
# =============================================================================

.make_panel_F <- function(datasets_all, datasets_filtered,
                           pathogen_results = NULL, best_dist = NULL,
                           x_max = 60L, n_draws = 100L,
                           base_size = 9) {

  filtered_names <- names(datasets_filtered)

  # Collect all frequency-table entries
  freq_entries <- list()
  for (nm in names(datasets_all)) {
    d    <- datasets_all[[nm]]
    type <- .detect_type(d)
    if (type %in% c("D", "E"))
      freq_entries[[nm]] <- list(d = d, type = type,
                                 included = nm %in% filtered_names)
  }
  if (length(freq_entries) == 0L) return(NULL)

  # ── Find the maximum data value across all freq datasets ────────────────────
  all_maxvals <- vapply(freq_entries, function(fe) {
    d <- fe$d
    if (fe$type == "D") max(d$freq_value) else max(d$freq_upper)
  }, numeric(1L))
  x_plot_max <- max(x_max, max(all_maxvals) * 1.05)

  # ── Compute fitted CDF once (filtered results, best-fit distribution) ───────
  fitted_cdf <- NULL
  dist_col   <- NULL
  dist_label <- NULL

  if (!is.null(pathogen_results) && !is.null(best_dist)) {
    res_filt <- pathogen_results[["filtered"]][[best_dist]]
    if (!is.null(res_filt) && !isTRUE(res_filt$skipped) && !is.null(res_filt$fit)) {
      cdf_dname <- .DIST_CDF_NAME[[best_dist]]
      x_seq     <- seq(0, x_plot_max, length.out = 300L)
      fitted_cdf <- tryCatch(
        compute_predictive_cdf(res_filt$fit, cdf_dname,
                               x_seq = x_seq, n_draws = n_draws),
        error = function(e) NULL
      )
      dist_col   <- .DIST_COLORS[[best_dist]]
      dist_label <- .DIST_LABELS[[best_dist]]
    }
  }

  # ── Build one CDF panel per freq-table dataset ──────────────────────────────
  panel_plots <- list()
  first_panel <- TRUE

  for (nm in names(freq_entries)) {
    fe   <- freq_entries[[nm]]
    d    <- fe$d
    type <- fe$type
    cite <- .short_cite(d$source)
    n_obs <- sum(d$freq_count)

    # Compute empirical CDF values
    # For interval-censored data (type E) we also compute the lower-bound ECDF
    # (at freq_lower) which gives the upper bound on the true CDF; the band
    # between the two ECDFs represents uncertainty about the true CDF.
    ecdf_lower_x <- NULL
    ecdf_lower_y <- NULL
    band_df      <- NULL

    if (type == "D") {
      sorted_idx  <- order(d$freq_value)
      ecdf_x      <- d$freq_value[sorted_idx]
      ecdf_counts <- d$freq_count[sorted_idx]
      x_caption   <- NULL
    } else {
      # Upper-bound ECDF (at freq_upper) = lower bound on true CDF
      sorted_upper <- order(d$freq_upper)
      ecdf_x       <- d$freq_upper[sorted_upper]
      ecdf_counts  <- d$freq_count[sorted_upper]

      # Lower-bound ECDF (at freq_lower) = upper bound on true CDF
      sorted_lower  <- order(d$freq_lower)
      ecdf_lower_x  <- d$freq_lower[sorted_lower]
      ecdf_lower_y  <- cumsum(d$freq_count[sorted_lower]) / sum(d$freq_count)

      x_caption <- paste0(
        "Interval-censored: black step = lower bound on true CDF (at interval upper end); ",
        "dashed step = upper bound on true CDF (at interval lower end); ",
        "band = uncertainty region"
      )
    }
    ecdf_y <- cumsum(ecdf_counts) / sum(ecdf_counts)

    # Polygon (filled area under step function)
    poly_df <- .ecdf_polygon(ecdf_x, ecdf_y,
                             x_start = 0, x_end = x_plot_max)

    # Uncertainty band polygon between the two ECDFs (type E only)
    if (!is.null(ecdf_lower_x)) {
      band_df <- .ecdf_band_polygon(
        top_x   = ecdf_lower_x, top_y = ecdf_lower_y,
        bot_x   = ecdf_x,       bot_y = ecdf_y,
        x_start = 0, x_end = x_plot_max
      )
    }

    # Step-line data (includes origin so the line starts at 0)
    step_df <- data.frame(x = c(0, ecdf_x), y = c(0, ecdf_y))

    panel_title <- paste0(
      if (first_panel) "F. " else "",
      cite, " (n=", n_obs, ")"
    )
    first_panel <- FALSE

    p <- ggplot2::ggplot() +
      # Fitted CDF ribbon (drawn first so ECDF sits on top)
      { if (!is.null(fitted_cdf))
          ggplot2::geom_ribbon(
            data = fitted_cdf$summary,
            ggplot2::aes(x = x, ymin = low, ymax = high),
            fill = dist_col, alpha = 0.15, colour = NA
          )
      } +
      # Uncertainty band between lower- and upper-bound ECDFs (type E only)
      { if (!is.null(band_df))
          ggplot2::geom_polygon(
            data = band_df,
            ggplot2::aes(x = x, y = y),
            fill = "#CCEBC5", alpha = 0.80, colour = NA
          )
      } +
      # Filled empirical step CDF polygon (upper-bound ECDF)
      ggplot2::geom_polygon(
        data = poly_df,
        ggplot2::aes(x = x, y = y),
        fill = "#B2DF8A", alpha = 0.55, colour = NA
      ) +
      # Empirical step CDF outline (upper-bound / lower bound on true CDF)
      ggplot2::geom_step(
        data = step_df,
        ggplot2::aes(x = x, y = y),
        colour = "black", linewidth = 0.75, direction = "hv"
      ) +
      # Lower-bound ECDF step line (upper bound on true CDF, type E only)
      { if (!is.null(ecdf_lower_x))
          ggplot2::geom_step(
            data = data.frame(x = c(0, ecdf_lower_x), y = c(0, ecdf_lower_y)),
            ggplot2::aes(x = x, y = y),
            colour = "#33A02C", linewidth = 0.60, linetype = "dashed", direction = "hv"
          )
      } +
      # Fitted CDF median line (drawn last so it is visible over the polygon)
      { if (!is.null(fitted_cdf))
          ggplot2::geom_line(
            data = fitted_cdf$summary,
            ggplot2::aes(x = x, y = median,
                         colour = dist_label, linetype = dist_label),
            linewidth = 0.8
          )
      } +
      # Colour/linetype scale for the fitted-CDF legend entry
      { if (!is.null(dist_label))
          list(
            ggplot2::scale_colour_manual(
              values   = setNames(dist_col, dist_label),
              name     = "Fitted CDF"
            ),
            ggplot2::scale_linetype_manual(
              values   = setNames("solid", dist_label),
              name     = "Fitted CDF"
            )
          )
      } +
      .cdf_y_scale() +
      ggplot2::scale_x_continuous(
        limits = c(0, x_plot_max), expand = c(0.01, 0)
      ) +
      ggplot2::labs(
        title   = panel_title,
        x       = "Days",
        y       = "Cumulative probability",
        caption = x_caption
      ) +
      .cdf_theme(base_size) +
      ggplot2::theme(
        legend.position = if (!is.null(fitted_cdf)) "right" else "none",
        legend.text     = ggplot2::element_text(size = base_size * 0.75),
        plot.caption    = ggplot2::element_text(size = base_size * 0.70,
                                                hjust = 0, colour = "grey50")
      )

    panel_plots[[nm]] <- p
  }

  if (length(panel_plots) == 0L) return(NULL)

  n_plots <- length(panel_plots)
  n_cols  <- if (n_plots >= 3L) 2L else min(n_plots, 2L)

  patchwork::wrap_plots(panel_plots, ncol = n_cols)
}


# =============================================================================
# Assembly: build the complete supplementary figure for one pathogen
# =============================================================================

#' Build the supplementary figure for a single pathogen
#'
#' @param pathogen      Character string: the pathogen key (e.g. "COVID_19").
#' @param pathogen_results  Named list: all_results[[pathogen]].
#' @param model_weights Named list from [compute_pathogen_model_bayes_factors()].
#' @param n_draws       Number of posterior draws for CDF computation.
#' @param base_size     Base font size passed to ggplot2 themes.
#'
#' @return A patchwork ggplot object, or NULL if no valid fits exist.
build_supp_figure <- function(pathogen,
                              pathogen_results,
                              model_weights,
                              n_draws   = 300L,
                              base_size = 9) {

  label <- .PATHOGEN_LABELS[[pathogen]]
  if (is.null(label)) label <- pathogen

  # ── Determine best-fit distribution (same logic as plot_main_figure) ────────
  weights   <- model_weights[[pathogen]]
  best_dist <- NULL

  if (!is.null(weights) && length(weights) > 0L) {
    for (.cand in names(weights)) {
      .r <- pathogen_results[["filtered"]][[.cand]]
      if (is.null(.r) || isTRUE(.r$skipped) || is.null(.r$fit)) next
      if (!.fit_has_converged(.r$fit)) next
      best_dist <- .cand
      break
    }
  }

  if (is.null(best_dist)) {
    message("  [SKIP] ", pathogen, ": no converged fit found.")
    return(NULL)
  }

  # ── Derive x_max from the best-fit distribution ─────────────────────────────
  x_max <- tryCatch(
    .auto_xmax(pathogen_results[["filtered"]][[best_dist]]$fit,
               .DIST_CDF_NAME[[best_dist]]),
    error = function(e) 60L
  )

  # ── Retrieve full and filtered dataset lists ─────────────────────────────────
  datasets_all  <- .get_datasets(pathogen_results[["all"]])
  datasets_filt <- .get_datasets(pathogen_results[["filtered"]])
  if (is.null(datasets_all))  datasets_all  <- list()
  if (is.null(datasets_filt)) datasets_filt <- list()

  # ── Build individual panels ──────────────────────────────────────────────────
  message("  Building panels...")

  panel_A <- tryCatch(
    .make_panel_A(pathogen_results, x_max, best_dist = best_dist,
                  n_draws = n_draws, base_size = base_size),
    error = function(e) { message("    [WARN] Panel A: ", conditionMessage(e)); NULL }
  )

  panel_B <- tryCatch(
    .make_panel_B(pathogen_results, best_dist, x_max, n_draws = n_draws, base_size = base_size),
    error = function(e) { message("    [WARN] Panel B: ", conditionMessage(e)); NULL }
  )

  panel_C <- tryCatch(
    .make_panel_C(pathogen_results, best_dist, x_max, n_draws = n_draws, base_size = base_size),
    error = function(e) { message("    [WARN] Panel C: ", conditionMessage(e)); NULL }
  )

  panel_D <- tryCatch(
    .make_panel_D(datasets_all, datasets_filt, base_size = base_size),
    error = function(e) { message("    [WARN] Panel D: ", conditionMessage(e)); NULL }
  )

  panel_F <- tryCatch(
    .make_panel_F(datasets_all, datasets_filt,
                  pathogen_results = pathogen_results,
                  best_dist        = best_dist,
                  x_max            = x_max,
                  n_draws          = ceiling(n_draws / 3L),
                  base_size        = base_size),
    error = function(e) { message("    [WARN] Panel F: ", conditionMessage(e)); NULL }
  )

  # ── Assemble with patchwork ──────────────────────────────────────────────────
  # Layout: [A | B] on the top row (B omitted if NULL → A takes full width)
  #          [C]    full-width (omitted if NULL)
  #          [D]    full-width (omitted if NULL)
  #          [F]    full-width (omitted if NULL)
  #
  # Heights are set relative to the number of datasets for D/F.

  n_summary_ds <- sum(vapply(datasets_all,
                             function(d) .detect_type(d) %in% c("A","B","C"),
                             logical(1L)))
  n_freq_ds    <- sum(vapply(datasets_all,
                             function(d) .detect_type(d) %in% c("D","E"),
                             logical(1L)))

  h_D <- if (!is.null(panel_D)) max(4, min(18, n_summary_ds * 0.70)) else 0
  h_F <- if (!is.null(panel_F)) max(4, min(16, n_freq_ds   * 3.0))   else 0

  # Row 1: Panel A (and optionally B)
  if (!is.null(panel_A) && !is.null(panel_B)) {
    row1 <- panel_A + panel_B + patchwork::plot_layout(ncol = 2L)
  } else if (!is.null(panel_A)) {
    row1 <- panel_A
  } else {
    row1 <- NULL
  }

  # Collect non-NULL rows in order
  rows   <- Filter(Negate(is.null), list(row1, panel_C, panel_D, panel_F))
  heights <- c(
    if (!is.null(row1))    6   else NULL,
    if (!is.null(panel_C)) 5.5 else NULL,
    if (!is.null(panel_D)) h_D else NULL,
    if (!is.null(panel_F)) h_F else NULL
  )

  if (length(rows) == 0L) return(NULL)

  fig <- patchwork::wrap_plots(rows, ncol = 1L,
                               heights = heights) +
    patchwork::plot_annotation(
      title = label,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face  = "bold",
          size  = base_size + 4L,
          hjust = 0.5,
          margin = ggplot2::margin(b = 6)
        )
      )
    )

  fig
}


# =============================================================================
# Main loop — generate and save one PDF per pathogen
# =============================================================================

message("Loading main results (this may take a moment for large files)...")
all_results <- readRDS(here::here("results", "main_results.rds"))

mw_path <- here::here("results", "model_weights.rds")
if (file.exists(mw_path)) {
  message("Loading cached model weights...")
  model_weights <- readRDS(mw_path)
} else {
  message("Computing model weights (LOO-CV)...")
  model_weights <- compute_pathogen_model_bayes_factors(all_results)
  saveRDS(model_weights, mw_path)
  message("  Saved model weights to: ", mw_path)
}

output_dir <- here::here("results", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", output_dir)

for (pathogen in names(all_results)) {
  message("\n", strrep("-", 60))
  message("Pathogen: ", pathogen)

  fig <- tryCatch(
    build_supp_figure(
      pathogen         = pathogen,
      pathogen_results = all_results[[pathogen]],
      model_weights    = model_weights,
      n_draws          = 300L,
      base_size        = 9
    ),
    error = function(e) {
      message("  [ERROR] ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fig)) {
    message("  No figure produced — skipping.")
    next
  }

  # Determine appropriate figure height from the assembled patchwork
  # (each panel's height was set in cm; total ≈ sum of heights + title margin)
  pdf_path <- file.path(output_dir, paste0(pathogen, "_supplementary.pdf"))

  # Collect panel info for this pathogen to set PDF dimensions
  datasets_all  <- tryCatch(
    .get_datasets(all_results[[pathogen]][["all"]]),
    error = function(e) list()
  )
  if (is.null(datasets_all)) datasets_all <- list()

  # Width: wider when A and B are shown side-by-side (filtering removed data)
  datasets_filt_names <- names(tryCatch(
    .get_datasets(all_results[[pathogen]][["filtered"]]),
    error = function(e) list()
  ))
  datasets_all_names <- names(datasets_all)
  fig_width_cm <- if (!setequal(datasets_all_names, datasets_filt_names)) 22 else 18

  has_C <- length(setdiff(names(all_results[[pathogen]]),
                          c("all", "filtered"))) > 0L
  has_D <- any(vapply(datasets_all,
                      function(d) .detect_type(d) %in% c("A","B","C"), logical(1L)))
  has_F <- any(vapply(datasets_all,
                      function(d) .detect_type(d) %in% c("D","E"), logical(1L)))

  n_summary_ds <- sum(vapply(datasets_all,
                             function(d) .detect_type(d) %in% c("A","B","C"), logical(1L)))
  n_freq_ds    <- sum(vapply(datasets_all,
                             function(d) .detect_type(d) %in% c("D","E"), logical(1L)))

  fig_height_cm <- 1.5 +   # title margin
    6.0 +                   # row1 (A ± B)
    (if (has_C) 5.5  else 0) +
    (if (has_D) max(4, min(18, n_summary_ds * 0.70)) else 0) +
    (if (has_F) max(4, min(16, n_freq_ds    * 3.0))  else 0)

  fig_height_cm <- max(fig_height_cm, 10)

  ggplot2::ggsave(
    filename = pdf_path,
    plot     = fig,
    width    = fig_width_cm,
    height   = fig_height_cm,
    units    = "cm",
    device   = "pdf"
  )
  message("  Saved: ", pdf_path)
}

message("\n", strrep("=", 60))
message("Done. Figures written to: ", output_dir)
message(strrep("=", 60))
