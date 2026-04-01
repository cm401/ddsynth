
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
  res_tmp <- res_out %>%
    left_join(dplyr::select(scenarios, scenario_name, summary_type, n_obs_mean,
                            n_obs_sd, n_obs_min, n_obs_max),
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
