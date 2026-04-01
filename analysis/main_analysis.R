# =============================================================================
# main_analysis.R
# -----------------------------------------------------------------------------
# Primary analysis script for the incubation period synthesis paper.
#
# For every pathogen the script runs:
#   (a) "all"      — full dataset, no filtering
#   (b) "filtered" — datasets flagged by pre_inference_checks() removed
#   (c) Any subgroup analyses defined in `subgroup_config` below
#
# Each of these is fitted with all three parametric distributions (log-normal,
# gamma, Weibull), giving 6+ fits per pathogen.
#
# Results are written to `results/main_results.rds` as a nested list:
#
#   all_results[[pathogen]][[analysis]][[distribution]]
#     $fit        rstan stanfit object
#     $stan_data  list passed to Stan
#     $datasets   the dataset list used (possibly filtered / subsetted)
#     $checks     output of pre_inference_checks() (NULL for subgroup runs)
#
# The file is saved after every completed pathogen so a crash does not
# lose earlier work.  Re-running the script skips any (pathogen, analysis,
# distribution) triple that already has a non-NULL entry.
# =============================================================================

library(ddsynth)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# ── 1. Sampling settings ──────────────────────────────────────────────────────

CHAINS  <- 4
ITER    <- 12000
WARMUP  <- 2000
THIN    <- 1
SEED    <- 123
CONTROL <- list(adapt_delta = 0.999, max_treedepth = 12, stepsize = 0.01)

DIST_CODES <- c(lognormal = 1L, gamma = 2L, weibull = 3L, burr = 4L)


# ── 2. Output path ────────────────────────────────────────────────────────────

OUTPUT_DIR  <- here::here("results")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "main_results.rds")
dir.create(OUTPUT_DIR, showWarnings = FALSE)


# ── 3. Pathogen registry ──────────────────────────────────────────────────────
# Add or remove pathogens here.  The name on the left is used as the key in
# the results list and in all output labels.

pathogen_registry <- list(
  Nipah          = datasets_Nipah,
  MVD            = datasets_MVD,
  EVD            = datasets_EVD,
  Lassa          = datasets_Lassa,
  SARS           = datasets_SARS,
  MERS           = datasets_MERS,
  Zika           = datasets_Zika,
  Measles        = datasets_Measles,
  Mpox           = datasets_Mpox,
  Cholera        = datasets_Cholera,
  RVF            = datasets_RVF,
  CCHF           = datasets_CCHF,
  CCHF_extended  = datasets_CCHF_extended
)


# ── 4. Subgroup analysis config ───────────────────────────────────────────────
# For each pathogen, define zero or more named subgroup analyses.
# Each entry is a list of arguments forwarded to filter_datasets():
#   subgroup  character vector of subgroup values to retain  (or NULL)
#   location  character vector of location values to retain  (or NULL)
#
# Example — if datasets_CCHF entries carry subgroup/location fields:
#
#   CCHF = list(
#     tick_bite  = list(subgroup = "tick-bite"),
#     nosocomial = list(subgroup = "nosocomial"),
#     turkey     = list(location = "Turkey")
#   )
#
# Leave empty list() for pathogens with no subgroup analyses.

subgroup_config <- list(
  Nipah         = list(),
  MVD           = list(),
  EVD           = list(),
  Lassa         = list(),
  SARS          = list(),
  MERS          = list(),
  Zika          = list(),
  Measles       = list(),
  Mpox          = list(),
  Cholera       = list(),
  RVF           = list(),
  CCHF          = list(
    # Uncomment and adjust once subgroup/location fields are populated:
    # tick_bite  = list(subgroup = "tick-bite"),
    # nosocomial = list(subgroup = "nosocomial")
  ),
  CCHF_extended = list()
)


# ── 5. Compile Stan model (once) ──────────────────────────────────────────────

stan_file  <- system.file("stan", "hierarchical_data_synthesis_summary_stats.stan",
                          package = "ddsynth")
stan_model <- rstan::stan_model(stan_file)


# ── 6. Helper: fit one (dataset list, distribution) combination ──────────────

.fit_one <- function(datasets, dist_name, pathogen, label) {
  dist_type <- DIST_CODES[[dist_name]]

  stan_data <- tryCatch(
    prepare_stan_data_from_datasets(datasets, dist_type = dist_type),
    error = function(e) {
      message("  [SKIP] prepare_stan_data failed for ",
              pathogen, "/", label, "/", dist_name, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(stan_data)) return(NULL)

  stan_data <- update_phi_prior(stan_data, datasets)

  fit <- tryCatch(
    rstan::sampling(
      stan_model,
      data    = stan_data,
      chains  = CHAINS,
      iter    = ITER,
      warmup  = WARMUP,
      thin    = THIN,
      control = CONTROL,
      seed    = SEED,
      refresh = 200
    ),
    error = function(e) {
      message("  [FAIL] sampling failed for ",
              pathogen, "/", label, "/", dist_name, ": ", conditionMessage(e))
      NULL
    }
  )

  list(fit = fit, stan_data = stan_data, datasets = datasets)
}


# ── 7. Main loop ──────────────────────────────────────────────────────────────

# Load existing results if resuming
all_results <- if (file.exists(OUTPUT_FILE)) {
  message("Resuming from existing results file: ", OUTPUT_FILE)
  readRDS(OUTPUT_FILE)
} else {
  list()
}

for (pathogen in names(pathogen_registry)) {

  message("\n", strrep("=", 70))
  message("PATHOGEN: ", pathogen)
  message(strrep("=", 70))

  datasets_full <- pathogen_registry[[pathogen]]
  if (is.null(all_results[[pathogen]])) all_results[[pathogen]] <- list()

  # ── Build the analysis table for this pathogen ──────────────────────────────
  # Each row: analysis label, dataset list, run pre_inference_checks?
  #
  # "all"      uses the full dataset list, no filtering
  # "filtered" uses the output of pre_inference_checks(filter = TRUE)
  # subgroups  use filter_datasets() then no additional check-filtering

  # Run pre_inference_checks once per pathogen (distribution-agnostic checks
  # are used only for the filter; distribution-specific checks are re-run
  # inside the loop below).
  analyses <- list()

  # "all" analysis
  analyses[["all"]] <- list(datasets = datasets_full, checks = NULL)

  # "filtered" analysis — run checks with lognormal as a proxy distribution
  # (the filter flags are based on moment estimates and LOO fits, which are
  # robust to the specific distribution choice for this purpose).
  message("\n  Running pre_inference_checks for filter (lognormal proxy)...")
  checks_proxy <- tryCatch(
    suppressWarnings(
      pre_inference_checks(
        datasets_full, stan_model,
        dist_type = 1L,   # lognormal as proxy
        verbose   = FALSE,
        filter    = TRUE
      )
    ),
    error = function(e) {
      message("  [WARN] pre_inference_checks failed: ", conditionMessage(e),
              " — 'filtered' analysis will use full dataset.")
      list(datasets = datasets_full)
    }
  )
  analyses[["filtered"]] <- list(datasets = checks_proxy$datasets,
                                  checks   = checks_proxy)

  n_removed <- length(datasets_full) - length(checks_proxy$datasets)
  if (n_removed > 0) {
    message("  Filter removed ", n_removed, " dataset(s): ",
            paste(setdiff(names(datasets_full),
                          names(checks_proxy$datasets)), collapse = ", "))
  } else {
    message("  Filter: no datasets removed.")
  }

  # Subgroup analyses
  sg_cfg <- subgroup_config[[pathogen]]
  if (!is.null(sg_cfg) && length(sg_cfg) > 0) {
    for (sg_name in names(sg_cfg)) {
      sg_datasets <- tryCatch(
        do.call(filter_datasets, c(list(datasets = datasets_full), sg_cfg[[sg_name]])),
        error = function(e) {
          message("  [WARN] filter_datasets failed for subgroup '", sg_name,
                  "': ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(sg_datasets) && length(sg_datasets) > 0) {
        analyses[[sg_name]] <- list(datasets = sg_datasets, checks = NULL)
        message("  Subgroup '", sg_name, "': ", length(sg_datasets), " dataset(s).")
      } else {
        message("  [SKIP] Subgroup '", sg_name, "' yielded 0 datasets.")
      }
    }
  }

  # ── Fit all (analysis × distribution) combinations ─────────────────────────
  for (analysis_label in names(analyses)) {

    if (is.null(all_results[[pathogen]][[analysis_label]]))
      all_results[[pathogen]][[analysis_label]] <- list()

    analysis_datasets <- analyses[[analysis_label]]$datasets

    if (length(analysis_datasets) == 0) {
      message("\n  [SKIP] ", analysis_label, ": 0 datasets remaining.")
      next
    }

    for (dist_name in names(DIST_CODES)) {

      # Skip if already done
      if (!is.null(all_results[[pathogen]][[analysis_label]][[dist_name]])) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label, " / ", dist_name,
                " — already in results.")
        next
      }

      message("\n  Fitting: ", pathogen, " | ", analysis_label,
              " | ", dist_name,
              "  (n_datasets = ", length(analysis_datasets), ")")

      result <- .fit_one(analysis_datasets, dist_name, pathogen, analysis_label)

      # Attach check output to result (only for the analyses that had one)
      if (!is.null(result)) {
        result$checks <- analyses[[analysis_label]]$checks
      }

      all_results[[pathogen]][[analysis_label]][[dist_name]] <- result
    }
  }

  # Save after each pathogen
  saveRDS(all_results, OUTPUT_FILE)
  message("\n  Saved results to: ", OUTPUT_FILE)
}

message("\n", strrep("=", 70))
message("All done.  Results written to: ", OUTPUT_FILE)
message(strrep("=", 70))
