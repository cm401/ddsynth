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
#     $tier, $max_rhat, $divergences, $converged, $runtime_secs
#                 diagnostics from fit_with_escalation() (R/utils.R)
#
# Fitting itself is queued across every pathogen first, then dispatched in
# one parallel batch (fit_corpus_parallel(), reserving RESERVE_CORES cores -
# see section 1) rather than one fit at a time. Each individual fit is
# additionally saved to its own file the moment it completes
# (results/main_analysis_tasks/), independent of the merged
# results/main_results.rds written once at the end - so a crash partway
# through only loses fits that hadn't completed yet, not the whole run, and
# re-running the script skips any (pathogen, analysis, distribution) triple
# whose task file (or, from before this restructuring, whose merged-results
# entry) already exists.
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
#
# Fitting now goes through fit_with_escalation() (R/utils.R), which tries
# increasingly conservative settings until Rhat<RHAT_TARGET and divergence
# rate <1%, falling back to the original fixed settings below (kept exactly
# as they were, not relaxed - they were originally needed for some chains to
# converge at all) as the final tier. There is no longer a single CHAINS/
# ITER/CONTROL to set here; SEED is still passed through directly.

SEED        <- 123
RHAT_TARGET <- 1.05

# Hard requirement, not a tuning default: leave this many cores free so the
# machine stays usable while a run is in progress. Do not reduce without
# being told to.
RESERVE_CORES <- 4

# Per-task wall-clock ceiling. A stuck fit never errors on its own; past
# this the task is killed and recorded as timed out so the run keeps moving.
TASK_TIMEOUT_SECS <- 90 * 60

# Per-dataset ceiling for pre_inference_checks()'s Check 5 LOO fits (a much
# smaller, fixed-settings single-dataset fit, not a full escalation-ladder
# corpus fit) - measured normal cost is ~4-5s/dataset, so 2 min gives ample
# headroom while still catching a fit stuck in a pathological region.
CHECKS_LOO_TIMEOUT_SECS <- 120

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

# Each (pathogen, analysis, distribution) fit is saved to its own file the
# moment it completes (fit_corpus_parallel(), R/utils.R), independent of
# whether/when the merged all_results list below is next written. This is
# the resume unit: a task whose file already exists is skipped even if a
# crash happened before the merge step ran.
TASK_DIR <- file.path(OUTPUT_DIR, "main_analysis_tasks")
dir.create(TASK_DIR, showWarnings = FALSE)

# pre_inference_checks() results (one per pathogen x hierarchical family),
# cached independently of TASK_DIR/main_results.rds so recovery of the
# "filtered" dataset list doesn't depend on whether the corresponding fit -
# or any fit at all - has completed yet.
CHECKS_DIR <- file.path(OUTPUT_DIR, "main_analysis_checks")
dir.create(CHECKS_DIR, showWarnings = FALSE)


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


# ── 5. Compile Stan models (once) ─────────────────────────────────────────────
# Set to "factorised" (default) or "joint" to switch between Stan models.
#
# Two models are compiled per STAN_MODEL choice: the hierarchical per-study
# `phi_d` model (used for lognormal, Weibull, Burr XII, and gengamma when
# should_attempt_gg() passes), and the pre-point-4 shared-scalar-`phi` model
# (used for gamma). See POINT4_LIKELIHOOD_MATHS.md Part E.5: gamma's shape
# parameter is a concentration-type parameter whose per-study information
# vanishes as phi grows, creating a genuine hierarchical funnel (not just a
# prior-calibration problem) rather than a family that can safely take the
# phi_d/omega reparameterization.
STAN_MODEL <- "factorised"
if (!STAN_MODEL %in% c("factorised", "joint"))
  stop('STAN_MODEL must be "factorised" or "joint", got: "', STAN_MODEL, '"')

stan_model_file <- switch(STAN_MODEL,
  factorised = "hierarchical_data_synthesis_summary_stats.stan",
  joint      = "hierarchical_data_synthesis_summary_stats_joint.stan"
)
stan_model_shared_phi_file <- switch(STAN_MODEL,
  factorised = "hierarchical_data_synthesis_summary_stats_shared_phi.stan",
  joint      = "hierarchical_data_synthesis_summary_stats_joint_shared_phi.stan"
)

stan_file  <- system.file("stan", stan_model_file, package = "ddsynth")
if (!nzchar(stan_file))
  stop("Stan file '", stan_model_file, "' not found in ddsynth installation.")
stan_model <- rstan::stan_model(stan_file)

stan_file_shared_phi <- system.file("stan", stan_model_shared_phi_file, package = "ddsynth")
if (!nzchar(stan_file_shared_phi))
  stop("Stan file '", stan_model_shared_phi_file, "' not found in ddsynth installation.")
stan_model_shared_phi <- rstan::stan_model(stan_file_shared_phi)

# Families using the hierarchical phi_d model vs. the shared-scalar-phi model
# (POINT4_LIKELIHOOD_MATHS.md Part E.5). Gen. gamma is hierarchical whenever
# it is attempted at all (the should_attempt_gg() gate below already decides
# whether to attempt it, per Part E.5's "two sequential gates" requirement).
HIERARCHICAL_FAMILIES <- c("lognormal", "weibull", "burr", "gengamma")
SHARED_PHI_FAMILIES    <- c("gamma")


# ── 6. Helper: resolve which compiled model a family uses ────────────────────
# (data prep, phi-prior update, and the escalation ladder itself are now
# handled inside fit_corpus_parallel()/fit_with_escalation(), R/utils.R)

.model_for_family <- function(dist_name) {
  if (dist_name %in% SHARED_PHI_FAMILIES) stan_model_shared_phi else stan_model
}

.task_file <- function(pathogen, analysis_label, dist_name) {
  file.path(TASK_DIR, paste0(pathogen, "__", analysis_label, "__", dist_name, ".rds"))
}

.checks_file <- function(pathogen, dist_name) {
  file.path(CHECKS_DIR, paste0(pathogen, "__", dist_name, ".rds"))
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

  # A force-rerun should also recompute pre_inference_checks(), not just the
  # fits - otherwise the cache below would keep serving the old filtered
  # dataset list straight through the "cleared" run.
  for (p in FORCE_RERUN) {
    stale_checks <- Sys.glob(file.path(CHECKS_DIR, paste0(p, "__*.rds")))
    if (length(stale_checks) > 0) {
      message("Force-rerun: clearing ", length(stale_checks),
              " cached pre_inference_checks() result(s) for '", p, "'.")
      unlink(stale_checks)
    }
  }
}


# ── 8. Build the analysis table for every pathogen, and one global task list ─
#
# Building `analyses` (all/filtered/subgroups) stays per-pathogen and serial,
# as before - it's comparatively cheap and each pathogen's "filtered" variant
# depends on that pathogen's own pre_inference_checks() result. What changes
# is the fitting step: rather than fitting each (pathogen, analysis, family)
# combination immediately and serially, every combination that still needs
# fitting is appended to one task list across ALL pathogens, dispatched in a
# single parallel batch at the end (section 9) so the concurrency-limited
# core budget stays continuously busy instead of idling between pathogens.

tasks <- list()

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

  # "filtered" analysis - family-specific (see POINT4_LIKELIHOOD_MATHS.md
  # Part E.10 for why a single shared lognormal proxy isn't enough). One
  # filtered list per HIERARCHICAL_FAMILIES member. gamma reuses "all"
  # unfiltered - it's protected instead by gamma_phi_extreme_safe().
  analyses[["filtered"]] <- list(by_family = list(gamma = analyses[["all"]]))

  for (dist_name in HIERARCHICAL_FAMILIES) {

    checks_file <- .checks_file(pathogen, dist_name)

    if (file.exists(checks_file)) {
      checks_res <- readRDS(checks_file)
      analyses[["filtered"]]$by_family[[dist_name]] <- list(datasets = checks_res$datasets, checks = checks_res)
      message("\n  Recovered filtered datasets for ", dist_name, " from cached checks file.")
      next
    }

    prior <- all_results[[pathogen]][["filtered"]][[dist_name]]
    if (!is.null(prior$datasets)) {
      analyses[["filtered"]]$by_family[[dist_name]] <- list(datasets = prior$datasets, checks = prior$checks)
      message("\n  Recovered filtered datasets for ", dist_name, " from existing results.")
      if (!is.null(prior$checks)) saveRDS(prior$checks, checks_file)
      next
    }

    message("\n  Running pre_inference_checks for filter (", dist_name, ")...")
    checks_res <- tryCatch(
      suppressWarnings(
        pre_inference_checks(
          datasets_full, stan_model,
          dist_type = DIST_CODES[[dist_name]],
          verbose   = FALSE,
          filter    = TRUE,
          loo_timeout_secs = CHECKS_LOO_TIMEOUT_SECS
        )
      ),
      error = function(e) {
        message("  [WARN] pre_inference_checks failed for ", dist_name, ": ", conditionMessage(e),
                " — 'filtered' will use the full dataset for this family.")
        list(datasets = datasets_full)
      }
    )
    saveRDS(checks_res, checks_file)
    analyses[["filtered"]]$by_family[[dist_name]] <- list(datasets = checks_res$datasets, checks = checks_res)
    n_removed <- length(datasets_full) - length(checks_res$datasets)
    if (n_removed > 0) {
      message("  Filter (", dist_name, ") removed ", n_removed, " dataset(s): ",
              paste(setdiff(names(datasets_full),
                            names(checks_res$datasets)), collapse = ", "))
    } else {
      message("  Filter (", dist_name, "): no datasets removed.")
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

  # ── Decide what needs fitting; queue it rather than fit immediately ────────
  for (analysis_label in names(analyses)) {

    if (is.null(all_results[[pathogen]][[analysis_label]]))
      all_results[[pathogen]][[analysis_label]] <- list()

    for (dist_name in names(DIST_CODES)) {

      # "filtered" is family-specific; other analyses aren't.
      analysis_datasets <- if (analysis_label == "filtered") {
        analyses[["filtered"]]$by_family[[dist_name]]$datasets
      } else {
        analyses[[analysis_label]]$datasets
      }

      if (length(analysis_datasets) == 0) {
        message("\n  [SKIP] ", analysis_label, " / ", dist_name, ": 0 datasets remaining.")
        next
      }

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

      # Gamma numerical-fragility check (see gamma_phi_extreme_safe() docs).
      if (dist_name == "gamma" &&
          !gamma_phi_extreme_safe(analysis_datasets, verbose = FALSE)) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label,
                " / gamma — a dataset's implied shape is extreme enough to",
                " risk grad_reg_lower_inc_gamma's slow-series regime",
                " (gamma_phi_extreme_safe() check failed).")
        all_results[[pathogen]][[analysis_label]][["gamma"]] <-
          list(skipped  = TRUE,
               reason   = "gamma_phi_extreme_safe",
               datasets = analysis_datasets)
        next
      }

      # Gamma shared-(mu0,tau,phi) weak-identifiability check (see
      # gamma_shared_phi_reliable() docs). Distinct from the two checks above:
      # not a hard boundary/slow-series issue, a genuinely flat/multimodal
      # surface diagnosed via Zika/all/gamma. Not predicted by sample size
      # alone (MVD/all/gamma has the same n=2 but converges cleanly), so this
      # runs unconditionally rather than being gated on n_datasets.
      if (dist_name == "gamma" &&
          !gamma_shared_phi_reliable(analysis_datasets, stan_model_shared_phi, verbose = FALSE)) {
        message("\n  [SKIP] ", pathogen, " / ", analysis_label,
                " / gamma — independent optimiser starts found no dominant",
                " mode in the shared (mu0, tau, phi) surface",
                " (gamma_shared_phi_reliable() check failed).")
        all_results[[pathogen]][[analysis_label]][["gamma"]] <-
          list(skipped  = TRUE,
               reason   = "gamma_shared_phi_reliable",
               datasets = analysis_datasets)
        next
      }

      # Burr XII per-study phi_d reliability check (see burr_phid_reliable() docs).
      task_stan_model <- .model_for_family(dist_name)
      if (dist_name == "burr" &&
          !burr_phid_reliable(analysis_datasets, stan_model, verbose = FALSE)) {
        message("\n  [INFO] ", pathogen, " / ", analysis_label,
                " / burr — falling back to shared-phi Burr XII ",
                "(burr_phid_reliable() check failed).")
        task_stan_model <- stan_model_shared_phi
      }

      out_file <- .task_file(pathogen, analysis_label, dist_name)
      task <- list(
        label       = paste(pathogen, analysis_label, dist_name, sep = "/"),
        pathogen    = pathogen,
        analysis_label = analysis_label,
        dist_name   = dist_name,
        datasets    = analysis_datasets,
        stan_model  = task_stan_model,
        output_file = out_file,
        checks      = if (analysis_label == "filtered") {
          analyses[["filtered"]]$by_family[[dist_name]]$checks
        } else {
          analyses[[analysis_label]]$checks
        }
      )

      if (file.exists(out_file)) {
        # Completed in an earlier (possibly interrupted) run of this script,
        # but not yet merged into all_results below - merge it, don't refit.
        message("\n  [RESUME] ", task$label, " — task file already exists, will merge without refitting.")
      } else {
        message("\n  [QUEUE] ", task$label, "  (n_datasets = ", length(analysis_datasets), ")")
      }
      tasks[[length(tasks) + 1]] <- task
    }
  }
}

# ── 9. Merge already-done tasks (from an earlier interrupted run) up front ───
# Tasks whose output file already exists (see [RESUME] above) are excluded
# from dispatch in section 10, but still need merging in - do that now so
# main_results.rds reflects everything on disk even before any new fitting
# happens this run.

for (task in tasks) {
  if (file.exists(task$output_file) && is.null(all_results[[task$pathogen]][[task$analysis_label]][[task$dist_name]])) {
    result <- readRDS(task$output_file)
    if (is.null(result$error)) result$checks <- task$checks
    all_results[[task$pathogen]][[task$analysis_label]][[task$dist_name]] <- result
  }
}
saveRDS(all_results, OUTPUT_FILE)

# ── 10. Dispatch every remaining task, saving as each one completes ─────────
# fit_corpus_parallel()'s on_complete callback fires in this (parent) process
# the moment each individual task finishes - not batched, not deferred to the
# end - so results/main_results.rds is never more than one fit's worth of
# work behind, and stopping the run part-way through loses at most whatever
# was still mid-fit, never anything already completed.

to_fit <- Filter(function(t) !file.exists(t$output_file), tasks)
message("\n", strrep("=", 70))
message(length(tasks), " total (pathogen, analysis, family) combinations queued; ",
        length(to_fit), " need fitting (", length(tasks) - length(to_fit), " already done).")
message(strrep("=", 70))

if (length(to_fit) > 0) {
  fit_corpus_parallel(
    to_fit, reserve_cores = RESERVE_CORES, rhat_target = RHAT_TARGET,
    task_timeout_secs = TASK_TIMEOUT_SECS,
    on_complete = function(task, result) {
      if (is.null(result$error)) result$checks <- task$checks
      all_results[[task$pathogen]][[task$analysis_label]][[task$dist_name]] <<- result
      saveRDS(all_results, OUTPUT_FILE)
      message("  Saved results to: ", OUTPUT_FILE, " (", task$label, ")")
    }
  )
}
message("\n", strrep("=", 70))
message("All done.  Results written to: ", OUTPUT_FILE)
message(strrep("=", 70))
