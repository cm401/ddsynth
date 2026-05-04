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
# Each of these is fitted with up to five parametric distributions (log-normal,
# gamma, Weibull, Burr, Generalised Gamma), giving 8+ fits per pathogen.
# The Generalised Gamma is skipped automatically when the identifiability
# heuristic (should_attempt_gg()) determines the data are insufficient.
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

# load_all() instead of library() so that unexported package objects
# (datasets_Nipah, datasets_SARS, etc.) are always loaded fresh from source.
# library(ddsynth) only puts *exported* symbols on the search path, so those
# dataset objects would resolve to whatever stale .GlobalEnv copy exists.
devtools::load_all(here::here(), quiet = TRUE)

# Remove any stale datasets_* objects that may be sitting in .GlobalEnv from a
# previous source() or load_all() call.  Because .GlobalEnv is searched before
# the package namespace, stale copies shadow the freshly-loaded versions and
# cause subgroup analyses to return 0 datasets (the old copies lack subgroup
# fields added since).
local({
  pkg_datasets <- grep("^datasets_",
                       ls(asNamespace("ddsynth")),
                       value = TRUE)
  stale <- intersect(ls(envir = .GlobalEnv), pkg_datasets)
  if (length(stale) > 0L) {
    message("Removing ", length(stale), " stale dataset object(s) from .GlobalEnv: ",
            paste(stale, collapse = ", "))
    rm(list = stale, envir = .GlobalEnv)
  }
})

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

DIST_CODES <- c(lognormal = 1L, gamma = 2L, weibull = 3L, burr = 4L,
                gengamma  = 5L)

# Pathogens listed here have their stored results wiped before the loop runs,
# forcing a full re-fit from scratch even if results already exist.
# Leave as character(0) for a normal resume run.
# Example: FORCE_RERUN <- c("SARS", "MERS")
FORCE_RERUN <- c()

# Distribution fits to skip unconditionally, regardless of data.
# Specified as a named list: pathogen key -> character vector of dist names.
# Example: FORCE_SKIP <- list(SARS = c("burr"), Nipah = c("gamma", "weibull"))
FORCE_SKIP <- list(
  Flu = c("gamma")
)

# ── 2. Output path ────────────────────────────────────────────────────────────

OUTPUT_DIR  <- here::here("results")
OUTPUT_FILE <- file.path(OUTPUT_DIR, "main_results.rds")
dir.create(OUTPUT_DIR, showWarnings = FALSE)


# ── 3. Pathogen registry ──────────────────────────────────────────────────────
# Add or remove pathogens here.  The name on the left is used as the key in
# the results list and in all output labels.

pathogen_registry <- list(
  Nipah          = datasets_Nipah,
  MVD            = datasets_MVD,        # really not enough data for a meaningful analysis, but included for completeness
  EVD            = datasets_EVD,
  Lassa          = datasets_Lassa,
  SARS           = datasets_SARS,
  MERS           = datasets_MERS,
  Zika           = datasets_Zika,
  Measles        = datasets_Measles,
  Mpox           = datasets_Mpox,
  Cholera        = datasets_Cholera,
  RVF            = datasets_RVF,
  #CCHF           = datasets_CCHF,
  CCHF           = datasets_CCHF_extended,
  COVID_19       = datasets_COVID_19,
  Dengue         = datasets_Dengue,
  YFV            = datasets_YFV,
  Typhoid        = datasets_typhoid,
  Smallpox       = datasets_Smallpox,
  Flu            = datasets_flu
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
  Nipah         = list(
    Bangladesh = list(subgroup = "Bangladesh")
  ),
  #MVD           = list(),   # insufficient data
  EVD           = list(
    WestAfrica    = list(subgroup = "West Africa"),
    nonWestAfrica = list(subgroup = "non-West Africa")
  ),
  Lassa         = list(),
  SARS          = list(
    HongKong  = list(subgroup = "Hong Kong"),
    Canada    = list(subgroup = "Canada"),
    Taiwan    = list(subgroup = "Taiwan"),
    China     = list(subgroup = "China"),
    Singapore = list(subgroup = "Singapore")
  ),
  MERS          = list(
    SaudiArabia = list(subgroup = "Saudi Arabia"),
    SouthKorea  = list(subgroup = "South Korea")
  ),
  Zika          = list(),
  Measles       = list(),
  Mpox          = list(),
  Cholera       = list(
    O1_El_Tor    = list(subgroup = "O1 El Tor"),
    O1_El_Tor_Ogawa = list(subgroup = "O1 El Tor Ogawa"),
    O139         = list(subgroup = "O139"),
    O1_Classical = list(subgroup = "O1 Classical")
  ),
  RVF           = list(),
  #CCHF          = list(),
  CCHF          = list(
    tick_bite  = list(subgroup = "tick-bite"),
    nosocomial = list(subgroup = "nosocomial")
  ),
  COVID_19      = list(
    Wildtype = list( subgroup = "Wildtype" ),
    Alpha    = list( subgroup = "Alpha" ),
    Delta    = list( subgroup = "Delta" ),
    Omicron  = list( subgroup = "Omicron" )
  ),
  Dengue        = list(
    inoculation = list(subgroup='inoculation'),
    mosquito    = list(subgroup='mosquito bite')
  ),
  YFV           = list(),
  Typhoid       = list(
    CateredMeal  = list(subgroup = 'Catered meal'),
    Experimental = list(subgroup = 'Experimental' ),
    Other        = list(subgroup = 'Other'),
    Water        = list(subgroup = 'Water')
  ),
  Smallpox      = list(),
  Flu           = list(
    H1N1 = list(subgroup = 'H1N1'),
    H3N2 = list(subgroup = 'H3N2'),
    H5N1 = list(subgroup = 'H5N1'),
    InfluenzaA = list(subgroup = 'Influenza A'),
    InfluenzaB = list(subgroup = 'Influenza B')
  )
)


# ── 5. Compile Stan model (once) ──────────────────────────────────────────────
# Set to "factorised" (default) or "joint" to switch between Stan models.
STAN_MODEL <- "factorised"
if (!STAN_MODEL %in% c("factorised", "joint"))
  stop('STAN_MODEL must be "factorised" or "joint", got: "', STAN_MODEL, '"')

stan_model_file <- switch(STAN_MODEL,
  factorised = "hierarchical_data_synthesis_summary_stats.stan",
  joint      = "hierarchical_data_synthesis_summary_stats_joint.stan"
)
stan_file  <- system.file("stan", stan_model_file, package = "ddsynth")
if (!nzchar(stan_file))
  stop("Stan file '", stan_model_file, "' not found in ddsynth installation.")
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


# ── 7. Load results and apply force-rerun ────────────────────────────────────

all_results <- if (file.exists(OUTPUT_FILE)) {
  message("Resuming from existing results file: ", OUTPUT_FILE)
  readRDS(OUTPUT_FILE)
} else {
  list()
}

# Wipe stored results for any pathogen listed in FORCE_RERUN so it is
# re-processed from scratch.
if (length(FORCE_RERUN) > 0) {
  unknown <- setdiff(FORCE_RERUN, names(pathogen_registry))
  if (length(unknown) > 0)
    warning("FORCE_RERUN contains unrecognised pathogen(s): ",
            paste(unknown, collapse = ", "), call. = FALSE)
  for (p in intersect(FORCE_RERUN, names(all_results))) {
    message("Force-rerun: clearing stored results for '", p, "'.")
    all_results[[p]] <- NULL
  }
}


# ── 8. Main loop ──────────────────────────────────────────────────────────────

for (pathogen in names(pathogen_registry)) {

  message("\n", strrep("=", 70))
  message("PATHOGEN: ", pathogen)
  message(strrep("=", 70))

  datasets_full <- pathogen_registry[[pathogen]]
  if (is.null(all_results[[pathogen]])) all_results[[pathogen]] <- list()

  # ── Build the analysis table for this pathogen ──────────────────────────────
  # For each analysis (all, filtered, subgroups) we need a dataset list.
  #
  # On a resume run the filtered dataset list and pre_inference_checks() output
  # are recovered from the previously stored results, avoiding the expensive
  # Stan optimisation that pre_inference_checks() calls internally.  The checks
  # are only re-run when no prior filtered result exists (first run, or after a
  # force-rerun).
  #
  # Subgroup dataset lists are similarly recovered from existing results, with
  # filter_datasets() called only when no prior result is available.

  analyses <- list()

  # "all" analysis — always the full dataset, nothing to recover
  analyses[["all"]] <- list(datasets = datasets_full, checks = NULL)

  # "filtered" analysis — recover or run fresh
  prior_filtered <- all_results[[pathogen]][["filtered"]]
  prior_filtered <- if (!is.null(prior_filtered)) Filter(Negate(is.null), prior_filtered) else list()

  if (length(prior_filtered) > 0) {
    ref <- prior_filtered[[1]]
    analyses[["filtered"]] <- list(datasets = ref$datasets, checks = ref$checks)
    message("\n  Recovered filtered datasets from existing results",
            " (skipping pre_inference_checks).")
  } else {
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
  }

  # Subgroup analyses — recover from existing results or build fresh
  sg_cfg <- subgroup_config[[pathogen]]
  if (!is.null(sg_cfg) && length(sg_cfg) > 0) {
    for (sg_name in names(sg_cfg)) {
      prior_sg <- all_results[[pathogen]][[sg_name]]
      prior_sg <- if (!is.null(prior_sg)) Filter(Negate(is.null), prior_sg) else list()

      if (length(prior_sg) > 0) {
        analyses[[sg_name]] <- list(datasets = prior_sg[[1]]$datasets, checks = NULL)
        message("  Subgroup '", sg_name, "': recovered from existing results.")
      } else {
        sg_datasets <- tryCatch(
          do.call(filter_datasets,
                  c(list(datasets = datasets_full), sg_cfg[[sg_name]])),
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
  }

  # ── Fast-skip if this pathogen is entirely complete ─────────────────────────
  # A slot is considered complete when it holds a non-NULL value — including the
  # sentinel list(skipped = TRUE, ...) stored below when the GG heuristic fires.
  all_done <- all(vapply(names(analyses), function(label) {
    all(vapply(names(DIST_CODES), function(dist_name) {
      !is.null(all_results[[pathogen]][[label]][[dist_name]])
    }, logical(1)))
  }, logical(1)))

  if (all_done) {
    message("  All fits already complete — nothing new to do.")
    next
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

      # Skip if already done (includes the skipped-GG sentinel)
      if (!is.null(all_results[[pathogen]][[analysis_label]][[dist_name]])) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label, " / ", dist_name,
                " — already in results.")
        next
      }

      # Manual force-skip: skip distributions listed in FORCE_SKIP for this
      # pathogen.  Stores the same sentinel structure as the automatic skips
      # so the slot is recognised as complete on future re-runs.
      if (dist_name %in% FORCE_SKIP[[pathogen]]) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label,
                " / ", dist_name, " — listed in FORCE_SKIP.")
        all_results[[pathogen]][[analysis_label]][[dist_name]] <-
          list(skipped  = TRUE,
               reason   = "force_skip",
               datasets = analysis_datasets)
        next
      }

      # GG identifiability check: skip if data cannot reliably identify the
      # extra shape parameter (kappa/Q).  Store a sentinel so this slot is
      # recognised as complete on future re-runs.
      if (dist_name == "gengamma" &&
          !should_attempt_gg(analysis_datasets, verbose = FALSE)) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label,
                " / gengamma — failed identifiability heuristic ",
                "(CV spread or information richness check).")
        all_results[[pathogen]][[analysis_label]][["gengamma"]] <-
          list(skipped  = TRUE,
               reason   = "identifiability_heuristic",
               datasets = analysis_datasets)
        next
      }

      # Gamma shape check: skip if type-2 (median+IQR) data imply a
      # concentrated distribution where central quantiles carry little
      # shape information, causing flat likelihoods and slow sampling.
      if (dist_name == "gamma" &&
          !gamma_type2_reliable(analysis_datasets, verbose = FALSE)) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label,
                " / gamma — type-2 datasets have high implied shape;",
                " central quantiles insufficient for reliable gamma shape",
                " inference (gamma_type2_reliable() check failed).")
        all_results[[pathogen]][[analysis_label]][["gamma"]] <-
          list(skipped  = TRUE,
               reason   = "gamma_type2_reliable",
               datasets = analysis_datasets)
        next
      }

      message("\n  Fitting: ", pathogen, " | ", analysis_label,
              " | ", dist_name,
              "  (n_datasets = ", length(analysis_datasets), ")")

      result <- .fit_one(analysis_datasets, dist_name, pathogen, analysis_label)

      if (!is.null(result))
        result$checks <- analyses[[analysis_label]]$checks

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
