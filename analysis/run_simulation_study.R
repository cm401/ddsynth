# ============================================================================
# analysis/run_simulation_study.R
# Hierarchical delay distribution synthesis — simulation study
#
# Distributions covered:
#   1 = lognormal   2 = gamma   3 = Weibull   4 = Burr XII   5 = Gen. Gamma
#
# Storage strategy — file-per-scenario:
#   results/by_scenario/<scenario_name>.rds  (one file, all n_sim rows)
#
#   * Interrupted runs resume automatically: the file is written after every
#     simulation replicate, so only the incomplete sim is re-run on restart.
#   * Adding new scenarios: define them in Section 3 and re-run the script;
#     already-complete scenarios are skipped automatically.
#   * Rerunning a single scenario: delete its .rds file, or set
#     FORCE_RERUN <- c("ScenarioName1", "ScenarioName2") (or TRUE for all).
#   * Final combined dataset: run collect_results() defined in Section 4.
# ============================================================================


# ── 0. Configuration ──────────────────────────────────────────────────────────

MAIN_FOLDER <- Sys.getenv("DDSYNTH_ROOT")
if (!nzchar(MAIN_FOLDER)) {
  # Default: project root = parent of the analysis/ directory this script lives in
  MAIN_FOLDER <- normalizePath(
    file.path(dirname(sys.frame(1)$ofile), ".."),
    winslash = "/", mustWork = FALSE
  )
}

RESULTS_DIR <- file.path(MAIN_FOLDER, "results", "by_scenario")
N_SIM       <- 100    # replications per scenario
SEED        <- 123    # global RNG seed (per-sim seeds derived from this)

# FORCE_RERUN: FALSE = skip completed scenarios (default resume behaviour).
#              TRUE  = rerun every scenario from scratch.
#              character vector = rerun only the named scenarios.
FORCE_RERUN <- FALSE

# RUN_ONLY: NULL = run all scenarios.
#           character vector = run only the named scenarios (useful for
#           testing a single scenario before committing to the full run).
RUN_ONLY <- NULL


# ── 1. Setup ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(purrr)
  library(tibble)
})

rstan_options(auto_write = TRUE)
options(mc.cores = 1L)   # Stan chains run sequentially within each scenario

source(file.path(MAIN_FOLDER, "R", "utils.R"))

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(MAIN_FOLDER, "results"), recursive = TRUE, showWarnings = FALSE)

cat("Compiling Stan model...\n")
stan_model <- rstan::stan_model(
  file.path(MAIN_FOLDER, "inst", "stan",
            "hierarchical_data_synthesis_summary_stats.stan")
)
cat("Done.\n\n")


# ── 2. Scenario-building helpers ──────────────────────────────────────────────

# Reconstruct the summary_type argument for generate_hierarchical_data_mixed
# from the summary_type_*_prop columns of a single scenario row.
summary_type_arg <- function(sc) {
  p <- c(sc$summary_type_1_prop, sc$summary_type_2_prop,
         sc$summary_type_3_prop, sc$summary_type_4_prop)
  idx <- which(p == 1)
  if (length(idx) == 1L) return(as.integer(idx))  # homogeneous
  p   # probability vector for mixed types
}

# n_obs for one scenario replicate (either fixed or randomly varied).
draw_n_obs <- function(sc, seed_offset) {
  if (!isTRUE(sc$vary_n)) return(sc$n_obs_mean)
  set.seed(seed_offset)
  n <- round(rnorm(sc$n_datasets, mean = sc$n_obs_mean, sd = sc$n_obs_sd))
  pmax(pmin(n, sc$n_obs_max), sc$n_obs_min)
}


# ── 3. Scenario definitions ────────────────────────────────────────────────────

# All scenarios come from generate_scenario_library(), which now covers
# dist_type 1–5.  Section 3 previously defined Burr XII and GG scenarios
# manually; those definitions now live inside the package function so that the
# simulation_study vignette and this script always use the same grid.

cat("Generating scenario library...\n")
all_scenarios <- generate_scenario_library(
  include_homogeneous   = TRUE,
  include_mixed         = TRUE,
  include_varied_n      = TRUE,
  include_freq_table    = TRUE,
  include_burr12        = TRUE,
  include_gengamma      = TRUE,
  include_gg_limitation = TRUE
)

cat(sprintf("  %d base (dist 1–3) scenarios\n",
            sum(startsWith(all_scenarios$scenario_group, "base_"))))
cat(sprintf("  %d Burr XII scenarios\n",
            sum(all_scenarios$scenario_group == "burr12")))
cat(sprintf("  %d GG standard scenarios\n",
            sum(all_scenarios$scenario_group == "gg_standard")))
cat(sprintf("  %d GG limitation scenarios\n",
            sum(all_scenarios$scenario_group == "gg_limitation")))
cat(sprintf("\nTotal scenarios: %d\n", nrow(all_scenarios)))

# Apply RUN_ONLY filter if specified
if (!is.null(RUN_ONLY)) {
  all_scenarios <- filter(all_scenarios, scenario_name %in% RUN_ONLY)
  cat(sprintf("Filtered to %d scenarios (RUN_ONLY)\n", nrow(all_scenarios)))
}


# ── 4. Runner functions ────────────────────────────────────────────────────────

# Run a single simulation replicate and return a one-row data frame whose
# columns are compatible with create_results_summary() / create_convergence_plot()
# in R/ploting_utils.R.
run_one_sim <- function(sc, sim_idx, stan_model, seed) {
  set.seed(seed + sc$scenario_idx * 10000L + sim_idx)
  n_obs_vec <- draw_n_obs(sc, seed + sc$scenario_idx * 10000L + sim_idx + 99999L)

  na_row <- tibble(
    scenario_name        = sc$scenario_name,
    scenario_group       = sc$scenario_group,
    scenario_idx         = sc$scenario_idx,
    sim                  = sim_idx,
    dist_type            = sc$dist_type,
    n_datasets           = sc$n_datasets,
    n_obs_mean           = sc$n_obs_mean,
    n_obs_sd             = sc$n_obs_sd,
    n_obs_min            = sc$n_obs_min,
    n_obs_max            = sc$n_obs_max,
    prop_summary_type_1  = sc$summary_type_1_prop,
    prop_summary_type_2  = sc$summary_type_2_prop,
    prop_summary_type_3  = sc$summary_type_3_prop,
    prop_summary_type_4  = sc$summary_type_4_prop,
    summary_type         = sc$summary_type,
    true_mu0 = sc$mu0, true_tau = sc$tau, true_phi = sc$phi, kappa = sc$kappa,
    max_rhat = NA_real_, min_neff = NA_real_, converged = FALSE,
    coverage_mu0 = NA_real_, bias_mu0 = NA_real_,
    coverage_tau = NA_real_, bias_tau = NA_real_,
    coverage_phi = NA_real_, bias_phi = NA_real_,
    coverage_kappa = NA_real_, bias_kappa = NA_real_, relbias_kappa = NA_real_,
    iqd = NA_real_
  )

  tryCatch({
    sim_data <- generate_hierarchical_data_mixed(
      n_datasets   = sc$n_datasets,
      n_obs        = n_obs_vec,
      dist_type    = sc$dist_type,
      mu0          = sc$mu0,
      tau          = sc$tau,
      phi          = sc$phi,
      kappa        = sc$kappa,
      summary_type = summary_type_arg(sc)
    )

    fit <- fit_model(sim_data, stan_model,
                     chains  = 4L,
                     refresh = 0L,
                     control = list(adapt_delta = 0.95, max_treedepth = 12L))

    fs       <- rstan::summary(fit)$summary
    max_rhat <- max(fs[, "Rhat"],  na.rm = TRUE)
    min_neff <- min(fs[, "n_eff"], na.rm = TRUE)

    tp        <- sim_data$true_params
    has_kappa <- sc$dist_type %in% c("burr12", "gengamma")

    bk <- if (has_kappa) compute_median_bias(fit, "kappa", tp$kappa) else NA_real_

    na_row |>
      mutate(
        max_rhat        = max_rhat,
        min_neff        = min_neff,
        converged       = (max_rhat <= 1.1) & (min_neff >= 100),
        coverage_mu0    = check_coverage(fit, "mu0", tp$mu0),
        bias_mu0        = compute_median_bias(fit, "mu0", tp$mu0),
        coverage_tau    = check_coverage(fit, "tau", tp$tau),
        bias_tau        = compute_median_bias(fit, "tau", tp$tau),
        coverage_phi    = check_coverage(fit, "phi", tp$phi),
        bias_phi        = compute_median_bias(fit, "phi", tp$phi),
        coverage_kappa  = if (has_kappa) check_coverage(fit, "kappa", tp$kappa)
                          else NA_real_,
        bias_kappa      = bk,
        relbias_kappa   = if (has_kappa) bk / tp$kappa else NA_real_,
        iqd             = NA_real_   # reserved for future compute_iqd() extension
      )
  }, error = function(e) {
    warning(sprintf("[scenario %s | sim %d] %s", sc$scenario_name, sim_idx, e$message))
    na_row
  })
}

# Run all n_sim replicates for one scenario row, saving incrementally.
# Returns invisibly; results are on disk at out_file.
run_scenario <- function(sc, stan_model, n_sim, seed, results_dir, force_rerun) {
  out_file   <- file.path(results_dir,
                          paste0(gsub("[^A-Za-z0-9_.-]", "_", sc$scenario_name), ".rds"))
  do_force   <- isTRUE(force_rerun) ||
                (is.character(force_rerun) && sc$scenario_name %in% force_rerun)

  # Determine which sims still need to run
  prev <- if (file.exists(out_file) && !do_force) readRDS(out_file) else tibble()
  done <- if (nrow(prev) > 0) max(prev$sim) else 0L

  if (done >= n_sim) {
    cat(sprintf("  [skip]   %s  (%d/%d done)\n", sc$scenario_name, done, n_sim))
    return(invisible(NULL))
  }
  if (done > 0L) {
    cat(sprintf("  [resume] %s from sim %d/%d\n", sc$scenario_name, done + 1L, n_sim))
  } else {
    cat(sprintf("  [start]  %s\n", sc$scenario_name))
    cat(sprintf("           dist=%s D=%d N=%.0f p=[%.2f %.2f %.2f %.2f] phi=%.2f kappa=%.2f\n",
                sc$dist_type, sc$n_datasets, sc$n_obs_mean,
                sc$summary_type_1_prop, sc$summary_type_2_prop,
                sc$summary_type_3_prop, sc$summary_type_4_prop,
                sc$phi, sc$kappa))
  }

  accum <- if (nrow(prev) > 0) list(prev) else list()

  for (sim in (done + 1L):n_sim) {
    cat(sprintf("    sim %d/%d\n", sim, n_sim))
    row      <- run_one_sim(sc, sim, stan_model, seed)
    accum    <- c(accum, list(row))
    saveRDS(bind_rows(accum), out_file)   # incremental save
  }

  invisible(NULL)
}

# Merge all per-scenario files into a single tibble.
collect_results <- function(results_dir) {
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    message("No result files found in ", results_dir)
    return(tibble())
  }
  purrr::map_dfr(files, readRDS)
}


# ── 5. Execution ──────────────────────────────────────────────────────────────

cat("\n========================================================\n")
cat("Running simulation study\n")
cat(sprintf("  Scenarios : %d\n", nrow(all_scenarios)))
cat(sprintf("  Sims/scen : %d\n", N_SIM))
cat(sprintf("  Results   : %s\n", RESULTS_DIR))
cat("========================================================\n\n")

t0 <- proc.time()

for (i in seq_len(nrow(all_scenarios))) {
  sc <- all_scenarios[i, ]
  cat(sprintf("\n[%d/%d] ", i, nrow(all_scenarios)))
  run_scenario(sc, stan_model,
               n_sim       = N_SIM,
               seed        = SEED,
               results_dir = RESULTS_DIR,
               force_rerun = FORCE_RERUN)
}

elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("\nAll scenarios finished in %.1f min\n", elapsed / 60))


# ── 6. Collect and save combined results ──────────────────────────────────────

cat("\nCollecting results...\n")
results_all <- collect_results(RESULTS_DIR)

out_combined <- file.path(MAIN_FOLDER, "results", "simulation_results_all.rds")
saveRDS(results_all, out_combined)
write.csv(results_all,
          sub("\\.rds$", ".csv", out_combined),
          row.names = FALSE)

cat(sprintf("Saved %d rows to %s\n", nrow(results_all), out_combined))

# Quick summary
cat("\nConvergence rates by group:\n")
results_all |>
  group_by(scenario_group) |>
  summarise(
    n_sims    = n(),
    pct_conv  = round(mean(converged, na.rm = TRUE) * 100, 1),
    med_rhat  = round(median(max_rhat, na.rm = TRUE), 3),
    cov_mu0   = round(mean(coverage_mu0,  na.rm = TRUE) * 100, 1),
    cov_phi   = round(mean(coverage_phi,  na.rm = TRUE) * 100, 1),
    cov_kappa = round(mean(coverage_kappa, na.rm = TRUE) * 100, 1)
  ) |>
  print()
