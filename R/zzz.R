# Suppress R CMD check NOTEs for column names used in dplyr/ggplot NSE.
# These variables are valid column names in the data frames passed to the
# plotting functions; they are not unbound global variables.
utils::globalVariables(c(
  # shared grouping / faceting variables
  "dist_type", "summary_type", "summary_type_label",
  "n_datasets", "n_datasets_bucket",
  "n_obs", "n_obs_bucket",
  "scenario_name", "scenario_group", "scenario_idx",
  # converged / diagnostic columns
  "converged", "max_rhat",
  # coverage / bias / mae columns
  "coverage", "coverage_mu0", "coverage_tau", "coverage_phi",
  "bias", "bias_mu0", "bias_tau", "bias_phi",
  "mae", "mae_mu0", "mae_tau", "mae_phi",
  # iqd / wis columns
  "iqd", "mean_iqd", "median_iqd", "sd_iqd",
  "wis", "rel_wis", "mean_wis", "median_wis", "sd_wis",
  # kappa coverage / bias columns
  "coverage_kappa", "bias_kappa", "rel_bias_kappa", "true_kappa",
  # predictive quantile coverage / bias columns (current names)
  "coverage_pred_median", "coverage_pred_q95",
  "bias_pred_median", "bias_pred_q95",
  # skipped / feasibility flags
  "skipped_reason", "scenario_feasibility_note",
  # parameter label after pivot_longer
  "parameter",
  # convergence rate (create_convergence_plot)
  "convergence_rate",
  # generate_cdf_plot_simulation_study_input
  "x", "y", "dist_label", "label",
  # generate_matched_moments_plot
  "density",
  # simulation / analysis result columns
  "analysis_", "cdf", "ci_95", "dataset", "datasets", "ESS",
  "high", "high_97.5", "i", "implied_phi", "is_outlier",
  "is_outlier_median", "is_outlier_spread",
  "log_prob", "low", "low_2.5",
  "median_bias_mu0", "median_bias_phi", "median_bias_tau",
  "model", "mu0",
  "n_obs_max", "n_obs_mean", "n_obs_min", "n_obs_sd",
  "newModel", "outside_prior_pi",
  "phi", "phi_hi", "phi_lo", "Rhat", "sd_est",
  "summaries_fun", "summary_names_used",
  "vary_n", "x_high", "x_low",
  # summary columns that share names with base/stats functions
  "median", "n", "sd"
))
