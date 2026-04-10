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
  # iqd
  "iqd", "mean_iqd",
  # parameter label after pivot_longer
  "parameter",
  # convergence rate (create_convergence_plot)
  "convergence_rate",
  # generate_cdf_plot_simulation_study_input
  "x", "y", "dist_label", "label",
  # generate_matched_moments_plot
  "density"
))
