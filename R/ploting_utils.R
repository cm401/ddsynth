
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
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scales = 'free_y') +
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
    facet_grid(summary_type_label ~ dist_type, labeller = label_value, scales = 'free_y') +
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


# =============================================================================
# Main analysis figure
# =============================================================================

# ── Internal constants ────────────────────────────────────────────────────────

# One colour per distribution (matches the vignette palette)
.DIST_COLORS <- c(
  lognormal = "#2166AC",   # blue
  gamma     = "#1A9641",   # green
  weibull   = "#D73027",   # red
  burr      = "#762A83",   # purple
  gengamma  = "#E08214"    # orange
)

# Human-readable distribution labels for in-panel annotation
.DIST_LABELS <- c(
  lognormal = "Log-normal",
  gamma     = "Gamma",
  weibull   = "Weibull",
  burr      = "Burr\u00a0XII",
  gengamma  = "Gen.\u00a0Gamma"
)

# compute_predictive_cdf() uses "gg" for generalised gamma, not "gengamma"
.DIST_CDF_NAME <- c(
  lognormal = "lognormal",
  gamma     = "gamma",
  weibull   = "weibull",
  burr      = "burr",
  gengamma  = "gg"
)

# Qualitative palette for subgroup lines (up to 8; ColorBrewer Set1 + extras)
.SUBGROUP_PALETTE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
  "#984EA3", "#A65628", "#F781BF", "#666666"
)

# Default pathogen display labels
.PATHOGEN_LABELS <- c(
  Nipah         = "Nipah",
  MVD           = "Marburg (MVD)",
  EVD           = "Ebola (EVD)",
  Lassa         = "Lassa fever",
  SARS          = "SARS",
  MERS          = "MERS",
  Zika          = "Zika",
  Measles       = "Measles",
  Mpox          = "Mpox",
  Cholera       = "Cholera",
  RVF           = "Rift Valley fever",
  CCHF          = "CCHF",
  COVID_19      = "COVID-19",
  Dengue        = "Dengue",
  YFV           = "Yellow fever",
  Typhoid       = "Typhoid",
  Smallpox      = "Smallpox",
  Flu           = "Influenza"
)

# ── Internal panel helpers ────────────────────────────────────────────────────

# Placeholder panel shown for pathogens whose results are not yet available
.pending_panel <- function(label, base_size = 9) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.55, label = label,
                      hjust = 0.5, vjust = 0.5, size = 3.2,
                      fontface = "bold", colour = "gray40") +
    ggplot2::annotate("text", x = 0.5, y = 0.42, label = "No results yet",
                      hjust = 0.5, vjust = 0.5, size = 2.4,
                      colour = "gray60", fontface = "italic") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "gray97", colour = "gray80", linewidth = 0.4
      )
    )
}

# Panel listing pathogens for which no systematic review data are available
.no_data_panel <- function(pathogen_names, base_size = 9) {
  n      <- length(pathogen_names)
  n_col1 <- ceiling(n / 2L)
  col1   <- pathogen_names[seq_len(n_col1)]
  col2   <- if (n > n_col1) pathogen_names[(n_col1 + 1L):n] else character(0L)

  y_start <- 0.86
  y_step  <- min(0.11, 0.80 / max(n_col1, 1L))

  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.97,
      label = "No systematic review data available",
      hjust = 0.5, vjust = 1, fontface = "bold",
      size = 2.8, colour = "gray20"
    ) +
    ggplot2::annotate(
      "segment", x = 0.05, xend = 0.95, y = 0.90, yend = 0.90,
      colour = "gray70", linewidth = 0.3
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(
        fill = "gray97", colour = "gray80", linewidth = 0.4
      )
    )

  for (i in seq_along(col1)) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.04, y = y_start - (i - 1L) * y_step,
      label = paste0("\u2022 ", col1[i]),
      hjust = 0, vjust = 1, size = 2.3, colour = "gray30"
    )
  }
  for (i in seq_along(col2)) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.52, y = y_start - (i - 1L) * y_step,
      label = paste0("\u2022 ", col2[i]),
      hjust = 0, vjust = 1, size = 2.3, colour = "gray30"
    )
  }
  p
}

# Internal helper: TRUE when all monitored parameters have Rhat <= threshold.
# Returns FALSE on any error (treat uncertain convergence as non-converged).
.fit_has_converged <- function(fit, rhat_threshold = 1.05) {
  tryCatch({
    s     <- rstan::summary(fit)$summary
    rhats <- s[, "Rhat"]
    rhats <- rhats[!is.na(rhats)]
    if (length(rhats) == 0L) return(FALSE)
    max(rhats) <= rhat_threshold
  }, error = function(e) FALSE)
}


# Build the CDF panel for one pathogen.
#
# Shows the posterior predictive CDF (ribbon + median line) for the overall
# "filtered" analysis, with dashed/dotted segments marking the P50 and P95
# uncertainty ranges.  Subgroup median CDFs are overlaid as thinner lines in a
# qualitative colour palette; the overall analysis is identified by the thicker
# line and its ribbon — it does not appear in the legend to keep the legend
# compact.  The chosen distribution is annotated in the bottom-right corner.
.build_pathogen_panel <- function(
  pathogen_label,
  result_filtered,    # all_results[[p]][["filtered"]][[best_dist]]
  subgroup_results,   # named list: sg → all_results[[p]][[sg]][[best_dist]]
  best_dist,
  best_weight   = NULL,
  x_max         = NULL,   # NULL → auto-derive from 99th pct; numeric → use as-is
  n_draws       = 500L,
  max_subgroups = 8L,
  base_size     = 9,
  show_title    = TRUE,   # set FALSE to suppress the panel title
  point_highlights = NULL # numeric vector of x-days to mark with red stars
) {
  dist_col   <- .DIST_COLORS[[best_dist]]
  dist_label <- .DIST_LABELS[[best_dist]]
  cdf_dname  <- .DIST_CDF_NAME[[best_dist]]

  # ── Dynamic x_max: nearest 5-day step at or above the 99th percentile ─────
  if (is.null(x_max)) {
    cdf_coarse <- tryCatch(
      compute_predictive_cdf(
        result_filtered$fit, cdf_dname,
        x_seq   = seq(0, 150, length.out = 100L),
        n_draws = 50L
      ),
      error = function(e) NULL
    )
    x99 <- if (!is.null(cdf_coarse)) {
      idx <- which(cdf_coarse$summary$median >= 0.99)[1L]
      if (!is.na(idx)) cdf_coarse$summary$x[idx] else 60
    } else 60
    # Round up to the nearest 5-day boundary, then add an extra 5-day buffer
    x_max <- max(10, ceiling(x99 / 5) * 5 + 5L)
  }

  x_seq <- seq(0, x_max, length.out = 300L)

  # ── Overall (filtered) CDF ─────────────────────────────────────────────────
  cdf_ov <- compute_predictive_cdf(
    result_filtered$fit, cdf_dname,
    x_seq = x_seq, n_draws = n_draws
  )
  q_df <- extract_quantiles(
    cdf_ov$summary, probs = c(0.5, 0.95), cdf_mat = cdf_ov$cdf_mat
  )
  overall_df <- dplyr::mutate(cdf_ov$summary, analysis_ = "Overall")

  # ── Subgroup CDFs ──────────────────────────────────────────────────────────
  sg_show <- head(names(subgroup_results), max_subgroups)
  sg_dfs  <- list()
  for (sg in sg_show) {
    r <- subgroup_results[[sg]]
    if (is.null(r) || is.null(r$fit)) next
    sg_cdf <- tryCatch(
      compute_predictive_cdf(
        r$fit, cdf_dname,
        x_seq   = x_seq,
        n_draws = ceiling(n_draws / 2L)
      ),
      error = function(e) NULL
    )
    if (!is.null(sg_cdf))
      sg_dfs[[sg]] <- dplyr::mutate(sg_cdf$summary, analysis_ = sg)
  }

  n_sg   <- length(sg_dfs)
  all_df <- dplyr::bind_rows(c(list(overall_df), sg_dfs))
  all_df$analysis_ <- factor(all_df$analysis_,
                              levels = c("Overall", names(sg_dfs)))

  # ── Colour / size scales ───────────────────────────────────────────────────
  sg_cols   <- if (n_sg > 0L)
    setNames(.SUBGROUP_PALETTE[seq_len(n_sg)], names(sg_dfs))
  else
    character(0L)

  color_map <- c("Overall" = dist_col, sg_cols)
  lwd_map   <- c("Overall" = 0.9,  setNames(rep(0.55, n_sg), names(sg_dfs)))
  alpha_map <- c("Overall" = 1.0,  setNames(rep(0.80, n_sg), names(sg_dfs)))

  # ── Plot ───────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(all_df, ggplot2::aes(x = x)) +

    # 95 % credible ribbon — overall only
    ggplot2::geom_ribbon(
      data  = dplyr::filter(all_df, analysis_ == "Overall"),
      ggplot2::aes(ymin = low, ymax = high),
      fill  = dist_col, alpha = 0.15, colour = NA
    ) +

    # Median CDF lines for all analyses
    ggplot2::geom_line(
      ggplot2::aes(y         = median,
                   colour    = analysis_,
                   linewidth = analysis_,
                   alpha     = analysis_)
    ) +

    # P50 and P95 uncertainty segments (overall only, dashed)
    ggplot2::geom_segment(
      data = q_df,
      ggplot2::aes(x = x_low, xend = x_high,
                   y = quantile, yend = quantile),
      colour    = dist_col,
      linetype  = "dashed",
      linewidth = 0.35,
      alpha     = 0.70
    ) +

    # Highlighted observations (e.g. external case reports not used in inference)
    # Rendered as upward-pointing red arrows terminating at the CDF curve.
    {
      if (!is.null(point_highlights) && length(point_highlights) > 0L) {
        hi_y  <- stats::approx(overall_df$x, overall_df$median,
                               xout = point_highlights, rule = 2)$y
        hi_df <- data.frame(x = point_highlights, y = hi_y)
        ggplot2::geom_segment(
          data        = hi_df,
          ggplot2::aes(x = x, xend = x, y = y - 0.20, yend = y - 0.01),
          colour      = "red",
          linewidth   = 0.5,
          arrow       = ggplot2::arrow(
            length = ggplot2::unit(0.15, "cm"),
            type   = "closed"
          ),
          inherit.aes = FALSE
        )
      } else {
        NULL
      }
    } +

    # Scales — legend shows only subgroups; "Overall" is identified by the ribbon
    ggplot2::scale_colour_manual(
      values = color_map,
      breaks = names(sg_cols)
    ) +
    ggplot2::scale_linewidth_manual(values = lwd_map,   guide = "none") +
    ggplot2::scale_alpha_manual(     values = alpha_map, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1),
      labels = c("0", ".25", ".5", ".75", ".95", "1"),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max), expand = c(0.01, 0)
    ) +
    ggplot2::labs(
      title  = if (isTRUE(show_title)) pathogen_label else NULL,
      x      = "Days",
      y      = "Cumulative probability",
      colour = NULL
    ) +

    # Distribution label + parameters — bottom-right, colour-coded, italic.
    # Built as a plotmath expression (parse = TRUE) so Greek letters (phi, kappa)
    # render correctly in every PDF device, not just cairo_pdf.
    {
      sims <- tryCatch(rstan::extract(result_filtered$fit), error = function(e) NULL)

      # Sanitise the distribution label: replace non-breaking spaces (\u00a0)
      # with ordinary spaces so the plotmath parser accepts them inside quotes.
      dl_safe <- gsub("\u00a0", " ", dist_label)

      if (!is.null(sims)) {
        # Line 1: distribution name (quoted plain text in plotmath)
        l1 <- sprintf("'%s'", dl_safe)

        # Line 2: mean (optional; quoted plain text)
        l2 <- if (!is.null(sims$pred_mean))
          sprintf("'Mean %.1f d'", stats::median(sims$pred_mean))
        else NULL

        # Line 3: shape parameters using plotmath Greek symbols.
        #   phi   → φ  (all distributions)
        #   kappa → κ  (Burr XII and GG only)
        has_phi   <- !is.null(sims$phi)
        has_kappa <- best_dist %in% c("burr", "gengamma") && !is.null(sims$kappa)
        l3 <- if (has_phi && has_kappa)
          sprintf("phi==%.2f~','~kappa==%.2f",
                  stats::median(sims$phi), stats::median(sims$kappa))
        else if (has_phi)
          sprintf("phi==%.2f", stats::median(sims$phi))
        else NULL

        # Stack non-NULL lines using nested atop(), then wrap in italic().
        # Reduce builds: atop(atop(l1, l2), l3) for 3 lines, atop(l1, l2) for 2, etc.
        active <- Filter(Negate(is.null), list(l1, l2, l3))
        stacked <- Reduce(function(acc, x) sprintf("atop(%s, %s)", acc, x), active)
        annot_label <- sprintf("italic(%s)", stacked)
      } else {
        annot_label <- sprintf("italic('%s')", dl_safe)
      }

      ggplot2::annotate(
        "text",
        x = x_max * 0.97, y = 0.02,
        label  = annot_label,
        colour = dist_col,
        size   = 2.1,
        hjust  = 1, vjust = 0,
        lineheight = 1.1,
        parse  = TRUE
      )
    }

  # Model weight annotation — only when one distribution clearly dominates
  if (!is.null(best_weight) && best_weight >= 0.75) {
    p <- p + ggplot2::annotate(
      "text",
      x = x_max * 0.03, y = 0.97,
      label  = sprintf("w = %.2f", best_weight),
      colour = "gray50", size = 2.0, hjust = 0, vjust = 1
    )
  }

  p +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = base_size,
        margin = ggplot2::margin(b = 2)
      ),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(colour = "gray92"),
      # Subgroup legend: bottom-right, anchored at its bottom-right corner so
      # it sits just above the distribution / parameter annotation text.
      legend.position      = if (n_sg > 0L) c(0.97, 0.26) else "none",
      legend.justification = c(1, 0),
      legend.text        = ggplot2::element_text(size = base_size * 0.72),
      legend.key.size    = ggplot2::unit(0.28, "cm"),
      legend.key.width   = ggplot2::unit(0.45, "cm"),
      legend.background  = ggplot2::element_rect(
        fill = "white", colour = "gray80", linewidth = 0.3
      ),
      legend.margin      = ggplot2::margin(2, 4, 2, 4),
      axis.title         = ggplot2::element_text(size = base_size * 0.85),
      axis.text          = ggplot2::element_text(size = base_size * 0.75),
      plot.margin        = ggplot2::margin(3, 5, 3, 5)
    )
}


# ── Exported functions ────────────────────────────────────────────────────────

#' Compute LOO-based model weights for all pathogens
#'
#' @description
#' For each pathogen, computes Leave-One-Out cross-validation (LOO-CV) for
#' every valid distribution fit in the chosen `analysis` slot, then converts
#' the ELPD differences into pseudo-Bayes-factor model weights.
#'
#' **Why not bridge sampling?**  `bridgesampling::bridge_sampler()` requires
#' the compiled C++ Stan model to be present in memory.  When a `stanfit` is
#' saved to RDS and reloaded the internal model object is not serialised and
#' the call fails.  LOO-CV operates entirely on the stored MCMC draws
#' (`log_lik` parameter) so it works correctly on reloaded fits.
#'
#' **Method**: `loo::loo()` is called on the non-placeholder `log_lik` columns
#' of each fit (the Stan model stores zero-filled placeholder columns for
#' unused summary-statistic slots; these are stripped before passing to LOO).
#' The resulting ELPD estimates are converted to weights via softmax, giving
#' relative model probabilities under equal model priors — analogous to Bayes
#' factors computed from marginal likelihoods.
#'
#' @param all_results Nested list produced by `analysis/main_analysis.R`.
#' @param analysis Character scalar.  Which analysis slot to use
#'   (default `"filtered"`).
#'
#' @return A named list, one entry per pathogen.  Each entry is a named
#'   numeric vector of model weights (summing to 1) sorted in decreasing
#'   order, keyed by distribution name.  Pathogens with fewer than one valid
#'   fit return `NULL`.
#'
#' @seealso [plot_main_figure()]
#' @export
compute_pathogen_model_bayes_factors <- function(all_results, analysis = "filtered") {
  pathogens     <- names(all_results)
  result        <- vector("list", length(pathogens))
  names(result) <- pathogens

  for (pathogen in pathogens) {
    dist_results <- all_results[[pathogen]][[analysis]]
    if (is.null(dist_results)) next

    elpds <- c()

    for (dist_name in names(dist_results)) {
      r <- dist_results[[dist_name]]
      # Skip: NULL slot, skipped-GG sentinel, or failed fit
      if (is.null(r) || isTRUE(r$skipped) || is.null(r$fit)) next

      elpd <- tryCatch({
        # Extract log_lik matrix (draws × log_lik slots)
        ll_mat <- rstan::extract(r$fit, "log_lik")$log_lik

        # The Stan model stores zero-filled placeholder columns for unused
        # summary-statistic slots (e.g. slot 3 for mean+SD datasets, slots
        # 2 and 3 for frequency-table datasets).  Remove them before LOO so
        # they are not treated as observations with log-likelihood = 0.
        valid_cols <- which(colMeans(abs(ll_mat)) > 1e-10)
        if (length(valid_cols) == 0L) stop("no non-zero log_lik columns")

        lo <- loo::loo(ll_mat[, valid_cols, drop = FALSE])
        lo$estimates["elpd_loo", "Estimate"]
      }, error = function(e) {
        message("  [WARN] loo() failed for ",
                pathogen, "/", dist_name, ": ", conditionMessage(e))
        NA_real_
      })

      if (!is.na(elpd)) elpds[dist_name] <- elpd
    }

    if (length(elpds) == 0L) next
    if (length(elpds) == 1L) {
      result[[pathogen]] <- setNames(1.0, names(elpds))
      next
    }

    # Softmax of ELPD differences → pseudo-BF model weights
    elpd_c             <- elpds - max(elpds)
    w                  <- exp(elpd_c) / sum(exp(elpd_c))
    result[[pathogen]] <- sort(w, decreasing = TRUE)
  }

  result
}


#' Build the main incubation-period analysis figure
#'
#' @description
#' Produces a multi-panel figure with one tile per pathogen.  Each tile shows
#' the posterior predictive CDF (95 % credible ribbon + median line) for the
#' best-fitting distribution selected by LOO-based model weights (pseudo-Bayes
#' factors; see [compute_pathogen_model_bayes_factors()]).  Dashed segments at the P50 and P95 marks span the 95 % posterior
#' uncertainty of those quantiles.
#'
#' When subgroup analyses are present (e.g. country or variant strata), their
#' median CDFs are overlaid as thin coloured lines without ribbons; the overall
#' "filtered" analysis remains the primary visual element.  The chosen
#' distribution is labelled in the bottom-right corner together with the
#' posterior median of the mean delay and shape parameter(s).  A model-weight
#' annotation (\eqn{w}) is shown in the top-left when one model clearly
#' dominates (\eqn{w \ge 0.75}).
#'
#' @param all_results Nested list produced by `analysis/main_analysis.R`.
#' @param pathogen_labels Optional named character vector overriding display
#'   labels (e.g. `c(COVID_19 = "COVID-19")`).
#' @param x_max Upper limit of the x-axis (days).  `NULL` (default) derives
#'   each panel's limit automatically as the nearest 5-day step at or above the
#'   99th percentile of the posterior predictive CDF (minimum 10 days).
#'   Alternatively, supply a single numeric applied to all panels, or a named
#'   numeric vector with per-pathogen overrides plus an optional `"default"`
#'   entry.
#' @param n_draws Integer.  Posterior draws for [compute_predictive_cdf()]
#'   for the overall analysis.  Subgroup CDFs use `ceiling(n_draws / 2)`.
#' @param model_weights Pre-computed output from [compute_pathogen_model_bayes_factors()].
#'   If `NULL`, weights are computed internally — this is slow (1–2 hours for
#'   a full set of pathogens).  Pre-compute and cache this object.
#' @param show_subgroups Logical.  Overlay subgroup median CDFs (default TRUE).
#' @param max_subgroups Integer.  Maximum subgroup lines per panel (default 8).
#' @param ncol Integer.  Columns in the assembled figure (default 4).
#' @param base_size Numeric.  Base font size for [ggplot2::theme_bw()]
#'   (default 9).
#'
#' @return A \pkg{patchwork} ggplot object.
#'
#' @seealso [compute_pathogen_model_bayes_factors()]
#' @export
plot_main_figure <- function(
  all_results,
  pathogens            = NULL,
  pathogen_labels      = NULL,
  x_max                = NULL,
  n_draws              = 500L,
  model_weights        = NULL,
  show_subgroups       = TRUE,
  max_subgroups        = 8L,
  ncol                 = 4L,
  base_size            = 9,
  show_title           = TRUE,
  pathogen_highlights  = NULL  # named list: pathogen key → numeric vector of x-days
) {

  # ── Pathogen subset filter ─────────────────────────────────────────────────
  if (!is.null(pathogens))
    all_results <- all_results[intersect(pathogens, names(all_results))]

  # ── Display labels ─────────────────────────────────────────────────────────
  labels <- .PATHOGEN_LABELS
  if (!is.null(pathogen_labels))
    labels[names(pathogen_labels)] <- pathogen_labels

  # ── Model weights ──────────────────────────────────────────────────────────
  if (is.null(model_weights)) {
    message(
      "No model_weights supplied — running compute_pathogen_model_bayes_factors() now.\n",
      "Consider pre-computing and caching: mw <- compute_pathogen_model_bayes_factors(all_results)\n",
      "This may take 1-2 hours for a full result set."
    )
    model_weights <- compute_pathogen_model_bayes_factors(all_results)
  }

  # ── Build per-pathogen panels ──────────────────────────────────────────────
  panels <- list()

  for (pathogen in names(all_results)) {
    label   <- labels[[pathogen]]
    if (is.null(label)) label <- pathogen

    weights <- model_weights[[pathogen]]
    if (is.null(weights) || length(weights) == 0L) {
      panels[[pathogen]] <- .pending_panel(label, base_size)
      next
    }

    # Walk through weight-ranked distributions until we find one that
    # (a) has a valid fit, and (b) has converged (max Rhat ≤ 1.05).
    # This guards against cases like Nipah where the top-weighted model
    # has poor mixing and produces implausible credible bands.
    best_dist    <- NULL
    best_wt      <- NULL
    filtered_res <- NULL

    for (.cand in names(weights)) {
      .r <- all_results[[pathogen]][["filtered"]][[.cand]]
      if (is.null(.r) || isTRUE(.r$skipped) || is.null(.r$fit)) next
      if (!.fit_has_converged(.r$fit)) {
        message("  [CONV] ", pathogen, "/", .cand,
                " — max Rhat > 1.05; skipping for figure.")
        next
      }
      best_dist    <- .cand
      best_wt      <- weights[[.cand]]
      filtered_res <- .r
      break
    }

    if (is.null(best_dist)) {
      message("  [WARN] ", pathogen, ": no converged fit found — using pending panel.")
      panels[[pathogen]] <- .pending_panel(label, base_size)
      next
    }

    # Resolve per-panel x_max
    # NULL → pass NULL to .build_pathogen_panel which auto-derives from 99th pct
    if (is.null(x_max)) {
      panel_xmax <- NULL
    } else if (length(x_max) > 1L) {
      panel_xmax <- x_max[[pathogen]]
      if (is.null(panel_xmax)) panel_xmax <- x_max[["default"]]
      if (is.null(panel_xmax)) panel_xmax <- NULL   # fall back to auto
    } else {
      panel_xmax <- x_max
    }

    # Subgroup results using the same distribution
    sg_results <- list()
    if (show_subgroups) {
      for (key in setdiff(names(all_results[[pathogen]]), c("all", "filtered"))) {
        r <- all_results[[pathogen]][[key]][[best_dist]]
        if (!is.null(r) && !isTRUE(r$skipped) && !is.null(r$fit))
          sg_results[[key]] <- r
      }
    }

    panels[[pathogen]] <- tryCatch(
      .build_pathogen_panel(
        pathogen_label   = label,
        result_filtered  = filtered_res,
        subgroup_results = sg_results,
        best_dist        = best_dist,
        best_weight      = best_wt,
        x_max            = panel_xmax,
        n_draws          = n_draws,
        max_subgroups    = max_subgroups,
        base_size        = base_size,
        show_title       = TRUE,
        point_highlights = if (!is.null(pathogen_highlights)) pathogen_highlights[[pathogen]] else NULL
      ),
      error = function(e) {
        message("  [WARN] Panel build failed for '", pathogen,
                "': ", conditionMessage(e))
        .pending_panel(label, base_size)
      }
    )
  }

  # ── Assemble with patchwork ────────────────────────────────────────────────
  patchwork::wrap_plots(panels, ncol = ncol) +
    patchwork::plot_annotation(
      title = if (isTRUE(show_title)) "Incubation period distributions - main analysis" else NULL,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = base_size + 3L,
          hjust = 0.5, margin = ggplot2::margin(b = 8)
        )
      )
    )
}

#' Split main figure into three panels by data availability
#'
#' Calls [plot_main_figure()] three times — once per pathogen group — and
#' assembles the results into a single vertically-stacked figure labelled A,
#' B, and C.
#'
#' @param all_results Nested results list (same as [plot_main_figure()]).
#' @param group_a Character vector of pathogen keys for panel A (high data).
#' @param group_b Character vector of pathogen keys for panel B (moderate data).
#' @param group_c Character vector of pathogen keys for panel C (low data).
#' @param panel_labels Length-3 character vector of panel labels.
#' @param model_weights Pre-computed model weights from
#'   [compute_pathogen_model_bayes_factors()].  If `NULL`, computed once and
#'   shared across all three panels.
#' @param ncol Number of columns within each panel group.
#' @param ... Additional arguments forwarded to [plot_main_figure()].
#' @return A patchwork figure.
#' @export
plot_main_figure_split <- function(
  all_results,
  group_a       = c("SARS", "COVID_19", "Flu", "Measles"),
  group_b       = c("Nipah", "EVD", "MERS", "CCHF",
                    "Cholera", "Typhoid", "Dengue", "YFV", "Mpox", "Smallpox"),
  group_c       = c("MVD", "Lassa", "Zika", "RVF"),
  panel_labels  = c("A", "B", "C"),
  model_weights = NULL,
  ncol          = 4L,
  ...
) {
  # Pre-compute model weights once to avoid 3x recomputation
  if (is.null(model_weights))
    model_weights <- compute_pathogen_model_bayes_factors(all_results)

  fig_a <- plot_main_figure(
    all_results, pathogens = group_a,
    model_weights = model_weights, ncol = ncol, show_title = FALSE, ...
  )
  fig_b <- plot_main_figure(
    all_results, pathogens = group_b,
    model_weights = model_weights, ncol = ncol, show_title = FALSE, ...
  )
  fig_c <- plot_main_figure(
    all_results, pathogens = group_c,
    model_weights = model_weights, ncol = ncol, show_title = FALSE, ...
  )

  # Heights proportional to number of rows in each group
  n_rows  <- function(g) ceiling(length(g) / ncol)
  heights <- c(n_rows(group_a), n_rows(group_b), n_rows(group_c))

  # wrap_elements() prevents patchwork from looking inside each group when
  # applying tags — ensures A/B/C label the groups, not individual sub-panels
  (patchwork::wrap_elements(fig_a) /
   patchwork::wrap_elements(fig_b) /
   patchwork::wrap_elements(fig_c)) +
    patchwork::plot_layout(heights = heights) +
    patchwork::plot_annotation(
      tag_levels = list(panel_labels),
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold")
      )
    )
}
