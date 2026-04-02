
#' @export
create_results_summary <- function(results, scenarios)
{
  summary_results <- results %>% filter(converged) %>%
    group_by(scenario_name, dist_type, n_datasets) %>%
    summarise( n_converged = n(), # Coverage (should be close to 0.95)
               coverage_mu0 = mean(coverage_mu0, na.rm = TRUE),
               coverage_tau = mean(coverage_tau, na.rm = TRUE),
               coverage_phi = mean(coverage_phi, na.rm = TRUE), # Median bias (should be close to 0)
               median_bias_mu0 = median(bias_mu0, na.rm = TRUE),
               median_bias_tau = median(bias_tau, na.rm = TRUE),
               median_bias_phi = median(bias_phi, na.rm = TRUE), # Mean absolute bias
               mae_mu0 = mean(abs(bias_mu0), na.rm = TRUE),
               mae_tau = mean(abs(bias_tau), na.rm = TRUE),
               mae_phi = mean(abs(bias_phi), na.rm = TRUE), # IQD
               mean_iqd = mean(iqd, na.rm = TRUE),
               median_iqd = median(iqd, na.rm = TRUE),
               sd_iqd = sd(iqd, na.rm = TRUE),.groups = "drop" )

  summary_results <- summary_results %>%
    left_join(dplyr::select(scenarios, scenario_name, scenario_group,
                            summary_type, n_obs_mean, n_obs_sd,
                            n_obs_min, n_obs_max, vary_n),
              by = "scenario_name")

  return(summary_results)
}

#' @export
create_coverage_plot <- function(summary_res)
{
  coverage_long <- summary_res %>%
    select(dist_type, summary_type, n_datasets, n_obs, coverage_mu0, coverage_tau, coverage_phi) %>%
    pivot_longer(cols = starts_with("coverage_"), names_to = "parameter", values_to = "coverage") %>%
    mutate(parameter = str_remove(parameter, "coverage_"),
           summary_type_label = factor(summary_type, levels = 1:5, labels = c("Median+Range", "Median+IQR", "Mean+SD","Freq Table", "Mixed"))) %>%
    mutate(
      n_obs_bucket = case_when(
        n_obs == 5  ~ "5",
        n_obs == 10 ~ "10",
        n_obs == 20 ~ "20",
        n_obs > 25  ~ "25+",
      ) %>% factor(levels = c("5", "10", "20", "25+")),
      n_datasets_bucket = case_when(
        n_datasets < 10 ~ "<10",
        n_datasets < 20 ~ "<20",
        n_datasets < 30 ~ "<30",
        n_datasets >= 30 ~ "30+",
      ) %>% factor(levels = c("<10", "<20", "<30","30+"))
    )


  plt <- ggplot(coverage_long, aes(x = parameter, y = coverage, color = n_datasets_bucket, shape = n_obs_bucket)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)) +
    scale_shape_manual(values = c("5" = 4, "10" = 3, "20" = 8, "25+" = 5)) +
    scale_color_aaas() +
    scale_x_discrete(labels = c("mu0" = expression(mu[0]), "phi" = expression(phi), "tau" = expression(tau))) +
    facet_grid(summary_type_label ~ dist_type, labeller = label_value) +
    labs(title = "Coverage of 95% Credible Intervals", subtitle = "Red line indicates nominal 95% coverage",
         x = "Parameter", y = "Empirical Coverage", color = "N datasets", shape = "N obs") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))

  return(plt)
}


#' @export
create_bias_plot <- function(summary_res)
{
  bias_long <- summary_res %>%
    select(dist_type, summary_type, n_datasets, n_obs, median_bias_mu0, median_bias_tau, median_bias_phi) %>%
    pivot_longer(cols = starts_with("median_bias_"), names_to = "parameter", values_to = "bias") %>%
    mutate(parameter = str_remove(parameter, "median_bias_"),
           summary_type_label = factor(summary_type, levels = 1:5, labels = c("Median+Range", "Median+IQR", "Mean+SD","Freq Table", "Mixed"))) %>%
    mutate(
      n_obs_bucket = case_when(
        n_obs == 5  ~ "5",
        n_obs == 10 ~ "10",
        n_obs == 20 ~ "20",
        n_obs > 25  ~ "25+",
      ) %>% factor(levels = c("5", "10", "20", "25+")),
      n_datasets_bucket = case_when(
        n_datasets < 10 ~ "<10",
        n_datasets < 20 ~ "<20",
        n_datasets < 30 ~ "<30",
        n_datasets >= 30 ~ "30+",
      ) %>% factor(levels = c("<10", "<20", "<30","30+"))
    )

  plt <- ggplot(bias_long, aes(x = parameter, y = bias, color = n_datasets_bucket, shape = n_obs_bucket)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)) +
    scale_shape_manual(values = c("5" = 4, "10" = 3, "20" = 8, "25+" = 5)) +
    scale_color_aaas() +
    scale_x_discrete(labels = c("mu0" = expression(mu[0]), "phi" = expression(phi), "tau" = expression(tau))) +
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scale='free_y') +
    labs(title = "Median Bias of Parameter Estimates", subtitle = "Red line indicates zero bias",
         x = "Parameter", y = "Median Bias", color = "N datasets", shape = "N obs") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))

  return(plt)
}

#' @export
create_mae_plot <- function(summary_res)
{
  bias_long <- summary_res %>%
    select(dist_type, summary_type, n_datasets, n_obs, mae_mu0, mae_tau, mae_phi) %>%
    pivot_longer(cols = starts_with("mae_"), names_to = "parameter", values_to = "mae") %>%
    mutate(parameter = str_remove(parameter, "mae_"),
           summary_type_label = factor(summary_type, levels = 1:5, labels = c("Median+Range", "Median+IQR", "Mean+SD","Freq Table", "Mixed"))) %>%
    mutate(
      n_obs_bucket = case_when(
        n_obs == 5  ~ "5",
        n_obs == 10 ~ "10",
        n_obs == 20 ~ "20",
        n_obs > 25  ~ "25+",
      ) %>% factor(levels = c("5", "10", "20", "25+")),
      n_datasets_bucket = case_when(
        n_datasets < 10 ~ "<10",
        n_datasets < 20 ~ "<20",
        n_datasets < 30 ~ "<30",
        n_datasets >= 30 ~ "30+",
      ) %>% factor(levels = c("<10", "<20", "<30","30+"))
    )

  plt <- ggplot(bias_long, aes(x = parameter, y = mae, color = n_datasets_bucket, shape = n_obs_bucket)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)) +
    scale_shape_manual(values = c("5" = 4, "10" = 3, "20" = 8, "25+" = 5)) +
    scale_color_aaas() +
    scale_x_discrete(labels = c("mu0" = expression(mu[0]), "phi" = expression(phi), "tau" = expression(tau))) +
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scale='free_y') +
    labs(title = "Mean Absolute Error of Parameter Estimates", subtitle = "Red line indicates Mean Absolute Error",
         x = "Parameter", y = "Mean Absolute Error", color = "N datasets", shape = "N obs") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))

  return(plt)
}

#' @export
create_iqd_plot <- function(summary_res)
{
  summary_res <- summary_res %>% mutate(summary_type_label = factor(summary_type, levels = 1:5, labels = c("Median+Range", "Median+IQR", "Mean+SD","Freq Table", "Mixed"))) %>%
    mutate(
      n_obs_bucket = case_when(
        n_obs == 5  ~ "5",
        n_obs == 10 ~ "10",
        n_obs == 20 ~ "20",
        n_obs > 25  ~ "25+",
      ) %>% factor(levels = c("5", "10", "20", "25+")),
      n_datasets_bucket = case_when(
        n_datasets < 10 ~ "<10",
        n_datasets < 20 ~ "<20",
        n_datasets < 30 ~ "<30",
        n_datasets >= 30 ~ "30+",
      ) %>% factor(levels = c("<10", "<20", "<30","30+"))
    )

  plt <- ggplot(summary_res, aes(x = n_obs_bucket, y = mean_iqd, color = n_datasets_bucket, shape = n_obs_bucket)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_shape_manual(values = c("5" = 4, "10" = 3, "20" = 8, "25+" = 5)) +
    scale_color_aaas() +
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scales = "free_y") +
    labs(title = "Integrated Quadratic Distance (IQD)",
         subtitle = "Lower values indicate better predictive performance",
         x = "Number of observations", y = "Mean IQD", color = "N datasets", shape = "N obs") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))

  return(plt)
}

#' @export
create_convergence_plot <- function(res_out, scenarios)
{
  # Join only summary_type from scenarios (n_obs_mean and the other n_obs_*
  # columns are already present in res_out from both the old and new runners).
  # Drop summary_type from res_out first to avoid .x/.y suffixes in case the
  # new runner already includes it.
  res_tmp <- res_out %>%
    dplyr::select(-dplyr::any_of("summary_type")) %>%
    left_join(dplyr::select(scenarios, scenario_name, summary_type),
              by = "scenario_name") %>%
    mutate(summary_type_label = factor(summary_type, levels = 1:5,
                                       labels = c("Median+Range", "Median+IQR",
                                                  "Mean+SD", "Freq Table", "Mixed")),
           n_obs = n_obs_mean) %>%
    mutate(
      n_obs_bucket = case_when(
        n_obs == 5  ~ "5",
        n_obs == 10 ~ "10",
        n_obs == 20 ~ "20",
        n_obs > 25  ~ "25+",
      ) %>% factor(levels = c("5", "10", "20", "25+")),
      n_datasets_bucket = case_when(
        n_datasets < 10 ~ "<10",
        n_datasets < 20 ~ "<20",
        n_datasets < 30 ~ "<30",
        n_datasets >= 30 ~ "30+",
      ) %>% factor(levels = c("<10", "<20", "<30","30+"))
    )


  convergence_summary <- res_tmp %>% group_by(scenario_idx, dist_type, summary_type_label, n_datasets_bucket, n_obs_bucket) %>%
    summarise( convergence_rate = mean(converged),.groups = "drop" ) %>%
    filter(!is.na(summary_type_label))

  plt <- ggplot(convergence_summary, aes(x = n_obs_bucket, y = convergence_rate, color = n_datasets_bucket, shape = n_obs_bucket)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_point(alpha = 0.7, size = 1.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)) +
    scale_shape_manual(values = c("5" = 4, "10" = 3, "20" = 8, "25+" = 5)) +
    scale_color_aaas() +
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scales = "free_y") +
    labs(title = "Convergence ratio of Bayesian data-synthesis models",
         x = "Number of observations", y = "Convergence ratio", color = "N datasets", shape = "N obs") +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1))

  return(plt)
}


#' Plot true CDFs for simulation study input scenarios
#'
#' Generates a faceted CDF plot for all distribution families used in the
#' simulation study, coloured by family (Lancet palette) and differentiated
#' by linetype for Burr XII and generalised gamma sub-scenarios. CDFs are
#' evaluated at the population mean (\code{loc = mu0}); between-study
#' heterogeneity (\code{tau = 0.4}) is not shown.
#'
#' @return A \code{ggplot} object.
#' @export
generate_cdf_plot_simulation_study_input <- function()
{
  # ── Parameter table ──────────────────────────────────────────────────────────
  scenarios <- tibble::tribble(
    ~label,                              ~dist_type,  ~mu0,    ~phi, ~kappa,
    "Lognormal",                         "lognormal", 2.00000, 0.35,   1.0,
    "Gamma",                             "gamma",     2.00000, 3.50,   1.0,
    "Weibull",                           "weibull",   2.00000, 2.00,   1.0,
    "Burr XII (\u03c6=2, \u03ba=2)",    "burr12",    1.94591, 2.00,   2.0,
    "Burr XII (\u03c6=3, \u03ba=2)",    "burr12",    1.94591, 3.00,   2.0,
    "Burr XII (\u03c6=2, \u03ba=5)",    "burr12",    1.94591, 2.00,   5.0,
    "Burr XII (\u03c6=3, \u03ba=5)",    "burr12",    1.94591, 3.00,   5.0,
    "Burr XII (\u03c6=2.5, \u03ba=3)",  "burr12",    1.94591, 2.50,   3.0,
    "GG (Q=0.1)",                        "gengamma",  1.94591, 0.50,   0.1,
    "GG (Q=0.5)",                        "gengamma",  1.94591, 0.50,   0.5,
    "GG (Q=1.0)",                        "gengamma",  1.94591, 0.50,   1.0,
    "GG (Q=2.0)",                        "gengamma",  1.94591, 0.50,   2.0,
    "GG (Q=3.0)",                        "gengamma",  1.94591, 0.50,   3.0
  ) |>
    mutate(dist_label = factor(
      dist_type,
      levels = c("lognormal","gamma","weibull","burr12","gengamma"),
      labels = c("Lognormal","Gamma","Weibull","Burr XII","Generalised Gamma")
    ))
  
  # ── CDF function ─────────────────────────────────────────────────────────────
  pcdf <- function(x, dist_type, mu0, phi, kappa) {
    switch(dist_type,
           lognormal = plnorm(x, meanlog = mu0, sdlog = phi),
           gamma     = pgamma(x, shape = phi, rate = phi / exp(mu0)),
           weibull   = pweibull(x, shape = phi, scale = exp(mu0)),
           burr12    = 1 - (1 + (x / exp(mu0))^phi)^(-kappa),
           gengamma  = {
             gs <- 1 / kappa^2
             u  <- exp((log(x) - mu0) * kappa / phi) / kappa^2
             pgamma(u, shape = gs, rate = 1)
           }
    )
  }
  
  # ── Build CDF data ────────────────────────────────────────────────────────────
  x_grid <- seq(0.1, 30, length.out = 600)
  
  cdf_data <- scenarios |>
    rowwise() |>
    mutate(cdf = list(tibble(x = x_grid,
                             y = pcdf(x_grid, dist_type, mu0, phi, kappa)))) |>
    ungroup() |>
    unnest(cdf)
  
  # ── Scales ───────────────────────────────────────────────────────────────────
  # Linetypes: single-line families → solid; multi-line families → varied
  lt_values <- c(
    "Lognormal"                         = "solid",
    "Gamma"                             = "solid",
    "Weibull"                           = "solid",
    "Burr XII (\u03c6=2, \u03ba=2)"    = "solid",
    "Burr XII (\u03c6=3, \u03ba=2)"    = "dashed",
    "Burr XII (\u03c6=2, \u03ba=5)"    = "dotted",
    "Burr XII (\u03c6=3, \u03ba=5)"    = "dotdash",
    "Burr XII (\u03c6=2.5, \u03ba=3)"  = "longdash",
    "GG (Q=0.1)"                        = "solid",
    "GG (Q=0.5)"                        = "dashed",
    "GG (Q=1.0)"                        = "dotted",
    "GG (Q=2.0)"                        = "dotdash",
    "GG (Q=3.0)"                        = "longdash"
  )
  
  # Only show Burr XII and GG sub-scenarios in the linetype legend
  lt_breaks <- c(
    "Burr XII (\u03c6=2, \u03ba=2)",  "Burr XII (\u03c6=3, \u03ba=2)",
    "Burr XII (\u03c6=2, \u03ba=5)",  "Burr XII (\u03c6=3, \u03ba=5)",
    "Burr XII (\u03c6=2.5, \u03ba=3)",
    "GG (Q=0.1)", "GG (Q=0.5)", "GG (Q=1.0)", "GG (Q=2.0)", "GG (Q=3.0)"
  )
  
  panel_labels <- c(
    lognormal = "Lognormal", gamma = "Gamma", weibull = "Weibull",
    burr12 = "Burr XII",    gengamma = "Generalised Gamma"
  )
  
  # ── Plot ─────────────────────────────────────────────────────────────────────
  ggplot(cdf_data, aes(x = x, y = y,
                       colour   = dist_label,
                       linetype = label)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(
      ~factor(dist_type, levels = c("lognormal","gamma","weibull","burr12","gengamma")),
      ncol     = 2,
      scales   = "free_x",
      labeller = as_labeller(panel_labels)
    ) +
    ggsci::scale_colour_lancet(name = "Distribution") +
    scale_linetype_manual(
      values = lt_values,
      breaks = lt_breaks,
      name   = "Sub-scenario"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, 0.25)) +
    labs(
      x       = "Delay (days)",
      y       = "Cumulative probability",
      title   = "True CDFs at population mean (\u03bc\u2080)",
      caption = "loc\u00a0=\u00a0\u03bc\u2080; between-study heterogeneity (\u03c4\u00a0=\u00a00.4) not shown"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position      = c(0.75, 0.17),   # centre of empty bottom-right cell
      legend.justification = c(0.5, 0.5),
      legend.box           = "vertical",
      legend.box.spacing   = unit(2, "pt"),
      legend.text          = element_text(size = 8),
      legend.title         = element_text(size = 9, face = "bold"),
      legend.background    = element_rect(fill = "white", colour = "grey80"),
      legend.key.width     = unit(1.2, "cm"),
      strip.background     = element_rect(fill = "grey92"),
      panel.grid.minor     = element_blank()
    ) +
    guides(
      colour   = guide_legend(
        order        = 1,
        nrow         = 1,          # single horizontal row of 5
        override.aes = list(linetype = "solid", linewidth = 1)
      ),
      linetype = guide_legend(
        order        = 2,
        ncol         = 2,          # 2 cols × 5 rows beneath the colour row
        override.aes = list(colour = "grey30")
      )
    )
}


#' Matched-moments distribution comparison plot
#'
#' Solves for the parameters of each distribution family such that all curves
#' share the same mean and coefficient of variation (CV) as the reference
#' Gamma scenario (mu0\,=\,2, phi\,=\,3.5). The resulting PDFs are plotted on
#' a faceted panel so that differences in distributional \emph{shape} — tail
#' heaviness, skewness, modality — are visible even when the first two moments
#' are identical.
#'
#' @param burr_kappas Numeric vector of kappa (k) values for Burr XII curves.
#'   Default: \code{c(2, 3, 5)}.
#' @param gg_kappas Numeric vector of Q values for the generalised gamma curves.
#'   Default: \code{c(0.1, 0.5, 1.0, 2.0, 3.0)}.
#' @return A \code{ggplot} object.
#' @export
generate_matched_moments_plot <- function(burr_kappas = 2,
                                          gg_kappas   = c(0.5, 1.5)) {

  # ── Target moments (Gamma: mu0=2, phi=3.5) ────────────────────────────────
  target_mean <- exp(2.0)
  target_cv   <- 1 / sqrt(3.5)   # ≈ 0.5345

  # ── CV functions (exp(mu0) cancels in the ratio) ──────────────────────────
  cv_weibull <- function(phi) {
    sqrt(gamma(1 + 2/phi) / gamma(1 + 1/phi)^2 - 1)
  }
  cv_burr12 <- function(phi, kappa) {
    m1 <- kappa * exp(lbeta(kappa - 1/phi, 1 + 1/phi))
    m2 <- kappa * exp(lbeta(kappa - 2/phi, 1 + 2/phi))
    sqrt(m2 / m1^2 - 1)
  }
  cv_gengamma <- function(phi, kappa) {
    gs      <- 1 / kappa^2
    log_cv2 <- lgamma(gs + 2*phi/kappa) + lgamma(gs) - 2*lgamma(gs + phi/kappa)
    sqrt(exp(log_cv2) - 1)
  }

  # ── Solve for (mu0, phi) per distribution ────────────────────────────────
  # Lognormal — closed form
  phi_ln  <- sqrt(log(1 + target_cv^2))
  mu0_ln  <- log(target_mean) - phi_ln^2 / 2

  # Gamma — already at reference
  phi_gam <- 1 / target_cv^2   # = 3.5
  mu0_gam <- log(target_mean)  # = 2.0

  # Weibull — numerical
  phi_w <- uniroot(function(phi) cv_weibull(phi) - target_cv, c(0.1, 20))$root
  mu0_w <- log(target_mean / gamma(1 + 1/phi_w))

  # Burr XII — one curve per kappa
  burr_list <- lapply(burr_kappas, function(kappa) {
    lb  <- 2 / kappa + 1e-6   # requires kappa*phi > 2 for variance
    phi <- uniroot(function(phi) cv_burr12(phi, kappa) - target_cv,
                   c(lb, 20))$root
    m1  <- kappa * exp(lbeta(kappa - 1/phi, 1 + 1/phi))
    mu0 <- log(target_mean) - log(m1)
    list(label = paste0("Burr XII (\u03ba=", kappa, ")"),
         dist_type = "burr12", mu0 = mu0, phi = phi, kappa = kappa)
  })

  # GG — one curve per Q (kappa)
  gg_list <- lapply(gg_kappas, function(kappa) {
    gs  <- 1 / kappa^2
    phi <- uniroot(function(phi) cv_gengamma(phi, kappa) - target_cv,
                   c(1e-4, 20))$root
    mu0 <- log(target_mean) -
      (2*phi/kappa) * log(kappa) -
      lgamma(gs + phi/kappa) + lgamma(gs)
    list(label = paste0("GG (Q=", kappa, ")"),
         dist_type = "gengamma", mu0 = mu0, phi = phi, kappa = kappa)
  })

  # ── Parameter table ───────────────────────────────────────────────────────
  list_to_row <- function(p)
    tibble::tibble(label = p$label, dist_type = p$dist_type,
                   mu0 = p$mu0, phi = p$phi, kappa = p$kappa)

  params <- dplyr::bind_rows(
    tibble::tibble(
      label     = c("Lognormal", "Gamma", "Weibull"),
      dist_type = c("lognormal", "gamma", "weibull"),
      mu0       = c(mu0_ln, mu0_gam, mu0_w),
      phi       = c(phi_ln, phi_gam, phi_w),
      kappa     = 1.0
    ),
    dplyr::bind_rows(lapply(burr_list, list_to_row)),
    dplyr::bind_rows(lapply(gg_list,   list_to_row))
  ) |>
    dplyr::mutate(
      dist_label = factor(
        dist_type,
        levels = c("lognormal","gamma","weibull","burr12","gengamma"),
        labels = c("Lognormal","Gamma","Weibull","Burr XII","Generalised Gamma")
      )
    )

  # ── PDF function ──────────────────────────────────────────────────────────
  dpdf <- function(x, dist_type, mu0, phi, kappa) {
    switch(dist_type,
      lognormal = dlnorm(x, meanlog = mu0, sdlog = phi),
      gamma     = dgamma(x, shape = phi, rate = phi / exp(mu0)),
      weibull   = dweibull(x, shape = phi, scale = exp(mu0)),
      burr12    = {
        lambda <- exp(mu0)
        phi * kappa / lambda * (x/lambda)^(phi-1) *
          (1 + (x/lambda)^phi)^(-(kappa+1))
      },
      gengamma  = {
        gs <- 1 / kappa^2
        w  <- (log(x) - mu0) / phi
        exp(log(kappa) - log(phi) - log(x) +
              gs*log(gs) + gs*kappa*w - gs*exp(kappa*w) - lgamma(gs))
      }
    )
  }

  # ── Build PDF data ─────────────────────────────────────────────────────────
  x_grid <- seq(0.01, 30, length.out = 600)

  pdf_data <- params |>
    dplyr::rowwise() |>
    dplyr::mutate(pdf = list(tibble::tibble(
      x       = x_grid,
      density = dpdf(x_grid, dist_type, mu0, phi, kappa)
    ))) |>
    dplyr::ungroup() |>
    tidyr::unnest(pdf)

  # ── Linetype mapping ───────────────────────────────────────────────────────
  lt_pool    <- c("solid","dashed","dotted","dotdash","longdash","twodash")
  burr_lbls  <- sapply(burr_list, `[[`, "label")
  gg_lbls    <- sapply(gg_list,   `[[`, "label")

  lt_values <- c(
    stats::setNames(rep("solid", 3), c("Lognormal","Gamma","Weibull")),
    stats::setNames(lt_pool[seq_along(burr_kappas)], burr_lbls),
    stats::setNames(lt_pool[seq_along(gg_kappas)],   gg_lbls)
  )

  # ── Plot ───────────────────────────────────────────────────────────────────
  ggplot2::ggplot(
    pdf_data,
    ggplot2::aes(x = x, y = density, colour = dist_label, linetype = label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(
      xintercept = target_mean, linetype = "dotted",
      colour = "grey50", linewidth = 0.4
    ) +
    ggsci::scale_colour_lancet(name = "Distribution") +
    ggplot2::scale_linetype_manual(
      values = lt_values,
      breaks = c(burr_lbls, gg_lbls),
      name   = "Sub-scenario"
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
    ggplot2::labs(
      x        = "Delay (days)",
      y        = "Density",
      title    = "Matched-moments distribution comparison",
      subtitle = sprintf(
        "All curves: mean\u00a0=\u00a0%.2f days, CV\u00a0=\u00a0%.3f (SD\u00a0\u2248\u00a0%.2f days)",
        target_mean, target_cv, target_mean * target_cv
      ),
      caption  = paste0(
        "Dotted line marks the common mean (",
        round(target_mean, 2), " days). ",
        "Parameters solved by moment matching."
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 9),
      legend.title     = ggplot2::element_text(size = 10, face = "bold"),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour   = ggplot2::guide_legend(
        order        = 1,
        override.aes = list(linetype = "solid", linewidth = 1)
      ),
      linetype = ggplot2::guide_legend(
        order        = 2,
        override.aes = list(colour = "grey30")
      )
    )
}
