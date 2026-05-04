# ============================================================================
# analysis/new_simulation_study.R
# Hierarchical delay distribution synthesis — extended simulation study
#
# Distributions covered:
#   1 = lognormal   2 = gamma   3 = Weibull   4 = Burr XII   5 = Gen. Gamma
#
# Extensions over run_simulation_study.R:
#   * Coverage and bias of the marginal MEDIAN and 95th PERCENTILE
#   * Weighted Interval Score (WIS) for overall predictive accuracy
#   * GG identifiability heuristic (should_attempt_gg): per-replicate skip
#   * Gamma type-2 reliability heuristic (gamma_type2_reliable): per-replicate skip
#   * scenario_feasibility_note: pre-annotates scenarios where a heuristic is
#     expected to fire, so skips can be distinguished from accidental failures
#   * Extended scenario grid:
#       - tau sensitivity (0.2, 0.6) × all 5 distributions, mixed summary types
#       - mu0 sensitivity (log(3), log(14)) × all 5 distributions, mixed summary types
#       - Extended GG scenarios with freq-table data (≥50% type-4) to ensure
#         the should_attempt_gg richness check passes
#       - Extended Burr XII type-2 (median+IQR) and type-4 (freq table) scenarios
#   * Results written to results/new_scenarios/
#
# Storage: file-per-scenario, same resume logic as run_simulation_study.R.
# ============================================================================


# ── 0. Configuration ──────────────────────────────────────────────────────────

MAIN_FOLDER <- Sys.getenv("DDSYNTH_ROOT")
if (!nzchar(MAIN_FOLDER)) {
  MAIN_FOLDER <- normalizePath(
    file.path(dirname(sys.frame(1)$ofile), ".."),
    winslash = "/", mustWork = FALSE
  )
}

RESULTS_DIR <- file.path(MAIN_FOLDER, "results", "new_scenarios")
N_SIM       <- 100    # replications per scenario
SEED        <- 123    # global RNG seed (per-sim seeds derived from this)

# STAN_MODEL: "factorised" (default) or "joint".
STAN_MODEL  <- "factorised"

# FORCE_RERUN: FALSE = resume; TRUE = rerun all; character vector = rerun named.
FORCE_RERUN <- FALSE

# RUN_ONLY: NULL = all; character vector = only these scenario names.
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
options(mc.cores = 1L)

source(file.path(MAIN_FOLDER, "R", "utils.R"))

dir.create(RESULTS_DIR,                       recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(MAIN_FOLDER, "results"), recursive = TRUE, showWarnings = FALSE)

cat("Compiling Stan model...\n")
stan_model_file <- switch(STAN_MODEL,
  factorised = "hierarchical_data_synthesis_summary_stats.stan",
  joint      = "hierarchical_data_synthesis_summary_stats_joint.stan"
)
stan_model <- rstan::stan_model(
  file.path(MAIN_FOLDER, "inst", "stan", stan_model_file)
)
cat("Done.\n\n")


# ── 2. Scenario-building helpers ──────────────────────────────────────────────

summary_type_arg <- function(sc) {
  p <- c(sc$summary_type_1_prop, sc$summary_type_2_prop,
         sc$summary_type_3_prop, sc$summary_type_4_prop)
  idx <- which(p == 1)
  if (length(idx) == 1L) return(as.integer(idx))
  p
}

draw_n_obs <- function(sc) {
  if (!isTRUE(sc$vary_n)) return(sc$n_obs_mean)
  n <- round(rnorm(sc$n_datasets, mean = sc$n_obs_mean, sd = sc$n_obs_sd))
  pmax(pmin(n, sc$n_obs_max), sc$n_obs_min)
}

# Build a one-row scenario data frame from explicit arguments.
# summary_props: named vector c(p1=, p2=, p3=, p4=) must sum to 1.
make_row <- function(name, group, dist, D, N, mu0, tau, phi, kappa,
                     summary_props) {
  p  <- summary_props
  st <- if (p["p1"] == 1) 1L else if (p["p2"] == 1) 2L else
        if (p["p3"] == 1) 3L else if (p["p4"] == 1) 4L else 5L
  data.frame(
    scenario_name       = name,
    scenario_group      = group,
    dist_type           = dist,
    n_datasets          = as.integer(D),
    mu0                 = mu0,
    tau                 = tau,
    phi                 = phi,
    kappa               = kappa,
    n_obs_mean          = N,
    n_obs_sd            = 0,
    n_obs_min           = N,
    n_obs_max           = N,
    summary_type_1_prop = unname(p["p1"]),
    summary_type_2_prop = unname(p["p2"]),
    summary_type_3_prop = unname(p["p3"]),
    summary_type_4_prop = unname(p["p4"]),
    summary_type        = st,
    vary_n              = FALSE,
    stringsAsFactors    = FALSE
  )
}

# Shared mixed-summary proportions used for sensitivity scenarios:
#   50% median+range, 10% median+IQR, 20% mean+SD, 20% freq table
SENS_PROPS <- c(p1 = 0.50, p2 = 0.10, p3 = 0.20, p4 = 0.20)

# GG sensitivity scenarios need ≥30% freq-table to pass should_attempt_gg.
# Use a 50/50 median+range / freq-table split.
GG_SENS_PROPS <- c(p1 = 0.50, p2 = 0.00, p3 = 0.00, p4 = 0.50)

# Default phi / kappa per distribution for sensitivity scenarios.
# These are representative "typical" values: moderate dispersion.
SENS_PHI   <- c(lognormal = 0.5,  gamma = 2.0,  weibull = 2.0,
                burr12    = 2.5,  gengamma = 0.5)
SENS_KAPPA <- c(lognormal = 1.0,  gamma = 1.0,  weibull = 1.0,
                burr12    = 3.0,  gengamma = 1.0)


# ── 3. Scenario definitions ────────────────────────────────────────────────────

# ── 3a. Base library (existing homogeneous / mixed / varied-N / Burr / GG) ───

cat("Generating base scenario library...\n")
all_scenarios <- generate_scenario_library(
  include_homogeneous   = TRUE,
  include_mixed         = TRUE,
  include_varied_n      = TRUE,
  include_freq_table    = TRUE,
  include_burr12        = TRUE,
  include_gengamma      = TRUE,
  include_gg_limitation = TRUE
)

cat(sprintf("  %d base (dist 1-3)\n",
            sum(startsWith(all_scenarios$scenario_group, "base_"))))
cat(sprintf("  %d Burr XII\n",   sum(all_scenarios$scenario_group == "burr12")))
cat(sprintf("  %d GG standard\n", sum(all_scenarios$scenario_group == "gg_standard")))
cat(sprintf("  %d GG limitation\n", sum(all_scenarios$scenario_group == "gg_limitation")))


# ── 3b. tau sensitivity scenarios ─────────────────────────────────────────────
#
# All 5 distributions × tau ∈ {0.2, 0.6} × mixed summary proportions.
# tau=0.4 is the baseline used throughout the base library.
# mu0 = log(7), D = 20, N = 40.
#
# Note: GG scenarios use GG_SENS_PROPS (50% type-1, 50% type-4) so that the
# should_attempt_gg richness check (≥30% freq table) passes.
# The base distributions (lognormal, gamma, weibull) and Burr XII use
# SENS_PROPS (20% type-4), which is fine for those distributions.

tau_rows <- do.call(rbind, lapply(
  c("lognormal", "gamma", "weibull", "burr12", "gengamma"),
  function(dist) {
    props <- if (dist == "gengamma") GG_SENS_PROPS else SENS_PROPS
    do.call(rbind, lapply(c(0.2, 0.6), function(tau_val) {
      make_row(
        name    = sprintf("TauSens_%s_tau%.1f_D20_N40", dist, tau_val),
        group   = "tau_sensitivity",
        dist    = dist,
        D       = 20, N = 40,
        mu0     = log(7),
        tau     = tau_val,
        phi     = SENS_PHI[dist],
        kappa   = SENS_KAPPA[dist],
        summary_props = props
      )
    }))
  }
))

cat(sprintf("  %d tau sensitivity\n", nrow(tau_rows)))
all_scenarios <- bind_rows(all_scenarios, tau_rows)


# ── 3c. mu0 sensitivity scenarios ─────────────────────────────────────────────
#
# All 5 distributions × mu0 ∈ {log(3), log(14)} × mixed summary proportions.
# log(3)  ≈ 1.10  →  mean delay ≈ 3 days   (short)
# log(14) ≈ 2.64  →  mean delay ≈ 14 days  (long)
# tau = 0.4, D = 20, N = 40.

mu0_rows <- do.call(rbind, lapply(
  c("lognormal", "gamma", "weibull", "burr12", "gengamma"),
  function(dist) {
    props <- if (dist == "gengamma") GG_SENS_PROPS else SENS_PROPS
    do.call(rbind, lapply(
      list(list(val = log(3),  lab = "Short"),
           list(val = log(14), lab = "Long")),
      function(m) {
        make_row(
          name    = sprintf("Mu0Sens_%s_%s_D20_N40", dist, m$lab),
          group   = "mu0_sensitivity",
          dist    = dist,
          D       = 20, N = 40,
          mu0     = m$val,
          tau     = 0.4,
          phi     = SENS_PHI[dist],
          kappa   = SENS_KAPPA[dist],
          summary_props = props
        )
      }
    ))
  }
))

cat(sprintf("  %d mu0 sensitivity\n", nrow(mu0_rows)))
all_scenarios <- bind_rows(all_scenarios, mu0_rows)


# ── 3d. Extended GG scenarios ─────────────────────────────────────────────────
#
# The base GG scenarios (from generate_scenario_library) all use summary type 1
# only (0% freq table), so they will always be skipped by should_attempt_gg's
# richness check (requires ≥30% freq table). They are kept in the grid to
# document that GG is infeasible without freq-table data.
#
# These new scenarios add freq-table data (≥50% type-4) to confirm that
# richer data makes GG identifiable:
#
#   Mix1_4: 50% median+range + 50% freq table
#   ST4:    100% freq table
#
# kappa values: 0.5 (near Weibull), 1.0 (lognormal), 2.0 (high right skew)
# sigma (phi) values: 0.3, 0.5, 0.7
# D: 5, 10, 20;  N: 10, 30, 100
# Limitation cases (kappa ≈ 0): also tested with freq tables

gg_extended <- rbind(

  # --- Core kappa/sigma grid with 50/50 mix ---------------------------------
  make_row("GG_Mix1_4_s0.5_Q0.5_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=0.5,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.5_Q1.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.5_Q2.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=2.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),

  # --- Varying sigma (phi) with kappa=1.0, 50/50 mix -----------------------
  make_row("GG_Mix1_4_s0.3_Q1.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.3, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.7_Q1.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.7, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),

  # --- Varying D and N with kappa=1.0, 50/50 mix ----------------------------
  make_row("GG_Mix1_4_s0.5_Q1.0_D5_N30",   "gg_extended", "gengamma",
           D= 5, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.5_Q1.0_D20_N30",  "gg_extended", "gengamma",
           D=20, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.5_Q1.0_D10_N10",  "gg_extended", "gengamma",
           D=10, N=10, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Mix1_4_s0.5_Q1.0_D10_N100", "gg_extended", "gengamma",
           D=10, N=100, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),

  # --- All-freq-table scenarios (strongest possible data) -------------------
  make_row("GG_ST4_s0.5_Q0.5_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=0.5,
           c(p1=0, p2=0, p3=0, p4=1)),
  make_row("GG_ST4_s0.5_Q1.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=1.0,
           c(p1=0, p2=0, p3=0, p4=1)),
  make_row("GG_ST4_s0.5_Q2.0_D10_N30", "gg_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=2.0,
           c(p1=0, p2=0, p3=0, p4=1)),

  # --- Limitation cases with freq-table data (does richer data help?) -------
  # kappa≈0 (near log-normal): hardest case; does 50% freq-table rescue it?
  make_row("GG_Lim_Mix1_4_Q0.1_D10_N30", "gg_limitation_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=0.1,
           c(p1=0.5, p2=0, p3=0, p4=0.5)),
  make_row("GG_Lim_ST4_Q0.1_D10_N30",    "gg_limitation_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=0.1,
           c(p1=0, p2=0, p3=0, p4=1)),
  # kappa=3.0 (high Q): does freq-table data confirm identifiability?
  make_row("GG_Lim_Mix1_4_Q3.0_D10_N30", "gg_limitation_extended", "gengamma",
           D=10, N=30, mu0=log(7), tau=0.4, phi=0.5, kappa=3.0,
           c(p1=0.5, p2=0, p3=0, p4=0.5))
)

cat(sprintf("  %d extended GG\n", nrow(gg_extended)))
all_scenarios <- bind_rows(all_scenarios, gg_extended)


# ── 3e. Extended Burr XII scenarios (type-2 and type-4) ───────────────────────
#
# The base Burr XII grid is almost entirely type-1 (median+range).
# These scenarios add type-2 (median+IQR) and type-4 (freq table) coverage
# across the key (phi=c, kappa=k) parameter combinations.

burr_ext <- rbind(

  # --- Type-2 (median + IQR) ------------------------------------------------
  make_row("Burr12_c2.5_k3_D10_N30_ST2",  "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=2.5, kappa=3.0,
           c(p1=0, p2=1, p3=0, p4=0)),
  make_row("Burr12_c2.5_k3_D20_N50_ST2",  "burr12", "burr12",
           D=20, N=50, mu0=log(7), tau=0.4, phi=2.5, kappa=3.0,
           c(p1=0, p2=1, p3=0, p4=0)),
  make_row("Burr12_c2_k2_D10_N30_ST2",    "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=2.0, kappa=2.0,
           c(p1=0, p2=1, p3=0, p4=0)),
  make_row("Burr12_c3_k5_D10_N30_ST2",    "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=3.0, kappa=5.0,
           c(p1=0, p2=1, p3=0, p4=0)),

  # --- Type-4 (freq table) --------------------------------------------------
  make_row("Burr12_c2.5_k3_D10_N30_ST4",  "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=2.5, kappa=3.0,
           c(p1=0, p2=0, p3=0, p4=1)),
  make_row("Burr12_c2.5_k3_D20_N50_ST4",  "burr12", "burr12",
           D=20, N=50, mu0=log(7), tau=0.4, phi=2.5, kappa=3.0,
           c(p1=0, p2=0, p3=0, p4=1)),
  make_row("Burr12_c2_k2_D10_N30_ST4",    "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=2.0, kappa=2.0,
           c(p1=0, p2=0, p3=0, p4=1)),
  make_row("Burr12_c3_k5_D10_N30_ST4",    "burr12", "burr12",
           D=10, N=30, mu0=log(7), tau=0.4, phi=3.0, kappa=5.0,
           c(p1=0, p2=0, p3=0, p4=1))
)

cat(sprintf("  %d extended Burr XII (type-2 / type-4)\n", nrow(burr_ext)))
all_scenarios <- bind_rows(all_scenarios, burr_ext)


# ── 3f. Finalise scenario table ───────────────────────────────────────────────

# Re-assign contiguous scenario indices (used as seed offsets in run_one_sim).
all_scenarios$scenario_idx <- seq_len(nrow(all_scenarios))

# ── scenario_feasibility_note ─────────────────────────────────────────────────
#
# Two-tier documentation of expected heuristic skips:
#
#   scenario_feasibility_note  (set here, before any run):
#     Explains WHY a scenario is expected to be skipped most or all of the time.
#     Recorded in every output row so results can be filtered/labelled cleanly.
#
#   skipped_reason  (set inside run_one_sim, after data generation):
#     Records which heuristic fired for that specific replicate.
#     NA means the replicate was fitted normally (or failed with an error).
#
# The two columns together tell the full story:
#   scenario_feasibility_note == "gg_insufficient_freq_table"
#     → should_attempt_gg richness check expected to fire every replicate
#   scenario_feasibility_note == "gamma_type2_may_be_slow"
#     → gamma_type2_reliable may fire depending on the observed IQR/median;
#       check skipped_reason to see which replicates were actually skipped
#   NA  → no known structural issue; fitted normally

all_scenarios$scenario_feasibility_note <- NA_character_

# GG scenarios with <30% freq-table data: richness check will always fire
gg_no_freq <- all_scenarios$dist_type == "gengamma" &
              all_scenarios$summary_type_4_prop < 0.30
all_scenarios$scenario_feasibility_note[gg_no_freq] <- "gg_insufficient_freq_table"

# Gamma scenarios that include type-2 data: may trigger gamma_type2_reliable
gamma_has_t2 <- all_scenarios$dist_type == "gamma" &
                all_scenarios$summary_type_2_prop > 0
all_scenarios$scenario_feasibility_note[gamma_has_t2] <- "gamma_type2_may_be_slow"

cat(sprintf("\nScenario feasibility notes:\n"))
table(all_scenarios$scenario_feasibility_note, useNA = "always") |>
  as.data.frame() |>
  setNames(c("note", "n")) |>
  print()

cat(sprintf("\nTotal scenarios: %d\n", nrow(all_scenarios)))

if (!is.null(RUN_ONLY)) {
  all_scenarios <- filter(all_scenarios, scenario_name %in% RUN_ONLY)
  cat(sprintf("Filtered to %d scenarios (RUN_ONLY)\n", nrow(all_scenarios)))
}


# ── 4. Runner functions ────────────────────────────────────────────────────────

# Convert the flat Stan-format obs_data list returned by
# generate_hierarchical_data_mixed() into a named list of per-dataset entries
# that match the structure expected by should_attempt_gg() and
# gamma_type2_reliable() (which were designed for the main analysis pipeline
# where each dataset is a separate named list with fields like $median, $Q1,
# $Q3, $mean, $sd, $freq_value, etc.).
#
# generate_hierarchical_data_mixed() returns a SINGLE flat list with vector
# fields (obs_stat1, obs_stat2, obs_stat3, summary_type, ...).  Passing that
# flat list directly to the heuristics causes vapply() to iterate over the
# top-level field names, none of which have $median/$Q1/etc., so every dataset
# appears as "not type-2" / "not rich" and the heuristics always return TRUE.
obs_data_to_dataset_list <- function(od) {
  n <- od$n_datasets
  lapply(seq_len(n), function(i) {
    st  <- od$summary_type[i]
    n_i <- od$n_obs[i]
    if (st == 1L) {          # median + range
      list(median = od$obs_stat1[i],
           min    = od$obs_stat2[i],
           max    = od$obs_stat3[i],
           n      = n_i)
    } else if (st == 2L) {  # median + IQR
      list(median = od$obs_stat1[i],
           Q1     = od$obs_stat2[i],
           Q3     = od$obs_stat3[i],
           n      = n_i)
    } else if (st == 3L) {  # mean + SD
      list(mean = od$obs_stat1[i],
           sd   = od$obs_stat2[i],
           n    = n_i)
    } else if (st == 4L) {  # frequency table
      start <- od$freq_start[i]
      len   <- od$freq_len[i]
      idx   <- seq_len(len) + start - 1L
      list(freq_value = od$freq_value[idx],
           freq_count = od$freq_count[idx],
           n          = n_i)
    } else {
      list(n = n_i)
    }
  })
}

# WIS probability grid (shared between WIS and posterior quantile CI).
WIS_ALPHA_LEVELS <- c(0.50, 0.80, 0.90, 0.95)
WIS_ALL_PROBS    <- sort(unique(c(
  (1 - WIS_ALPHA_LEVELS) / 2,        # 0.025 0.05 0.10 0.25
  0.5,
  1 - (1 - WIS_ALPHA_LEVELS) / 2    # 0.75  0.90 0.95 0.975
)))
# -> c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
IDX_MEDIAN <- which(abs(WIS_ALL_PROBS - 0.50) < 1e-9)  # 5
IDX_Q95    <- which(abs(WIS_ALL_PROBS - 0.95) < 1e-9)  # 8

# Run a single simulation replicate.
run_one_sim <- function(sc, sim_idx, stan_model, seed) {
  set.seed(seed + sc$scenario_idx * 10000L + sim_idx)
  n_obs_vec <- draw_n_obs(sc)

  # Base row with NAs for all computed columns.
  # scenario_feasibility_note is carried through unchanged from the scenario
  # table, providing static documentation in every output row.
  na_row <- tibble(
    scenario_name             = sc$scenario_name,
    scenario_group            = sc$scenario_group,
    scenario_idx              = sc$scenario_idx,
    scenario_feasibility_note = sc$scenario_feasibility_note,
    sim                       = sim_idx,
    dist_type                 = sc$dist_type,
    n_datasets                = sc$n_datasets,
    n_obs_mean                = sc$n_obs_mean,
    n_obs_sd                  = sc$n_obs_sd,
    n_obs_min                 = sc$n_obs_min,
    n_obs_max                 = sc$n_obs_max,
    prop_summary_type_1       = sc$summary_type_1_prop,
    prop_summary_type_2       = sc$summary_type_2_prop,
    prop_summary_type_3       = sc$summary_type_3_prop,
    prop_summary_type_4       = sc$summary_type_4_prop,
    summary_type              = sc$summary_type,
    true_mu0 = sc$mu0, true_tau = sc$tau, true_phi = sc$phi, kappa = sc$kappa,
    # skipped_reason: set when a heuristic fires; NA = normally fitted
    skipped_reason    = NA_character_,
    max_rhat          = NA_real_, min_neff  = NA_real_, converged = FALSE,
    # Hyperparameter coverage / bias
    coverage_mu0      = NA_real_, bias_mu0      = NA_real_,
    coverage_tau      = NA_real_, bias_tau      = NA_real_,
    coverage_phi      = NA_real_, bias_phi      = NA_real_,
    coverage_kappa    = NA_real_, bias_kappa    = NA_real_,
    relbias_kappa     = NA_real_,
    # Predictive quantile coverage / bias
    coverage_pred_median = NA_real_, bias_pred_median = NA_real_,
    coverage_pred_q95    = NA_real_, bias_pred_q95    = NA_real_,
    # Density / scoring
    iqd     = NA_real_,
    wis     = NA_real_,
    rel_wis = NA_real_
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

    obs_datasets <- sim_data$obs_data

    # Convert flat Stan-format obs_data to per-dataset list so that
    # should_attempt_gg() and gamma_type2_reliable() can inspect individual
    # dataset fields ($median, $Q1, $Q3, $freq_value, etc.).
    ds_list <- obs_data_to_dataset_list(obs_datasets)

    # ── Heuristic checks ────────────────────────────────────────────────────
    # Skips are per-replicate because they depend on the observed data.
    # skipped_reason records which heuristic fired; scenario_feasibility_note
    # records whether this was expected (see Section 3f above).

    if (sc$dist_type == "gengamma" &&
        !should_attempt_gg(ds_list, verbose = FALSE)) {
      return(mutate(na_row, skipped_reason = "gg_heuristic"))
    }

    if (sc$dist_type == "gamma" &&
        !gamma_type2_reliable(ds_list, verbose = FALSE)) {
      return(mutate(na_row, skipped_reason = "gamma_type2_heuristic"))
    }

    # ── Fit ─────────────────────────────────────────────────────────────────

    fit <- fit_model(sim_data, stan_model,
                     chains  = 4L,
                     refresh = 0L,
                     control = list(adapt_delta = 0.95, max_treedepth = 12L))

    fs       <- rstan::summary(fit)$summary
    max_rhat <- max(fs[, "Rhat"],  na.rm = TRUE)
    min_neff <- min(fs[, "n_eff"], na.rm = TRUE)

    tp        <- sim_data$true_params
    has_kappa <- sc$dist_type %in% c("burr12", "gengamma")

    # ── True marginal quantiles (50th and 95th) ──────────────────────────────

    true_qs <- compute_true_marginal_quantile(
      dist_type = sc$dist_type,
      mu0       = tp$mu0,
      tau       = tp$tau,
      phi       = tp$phi,
      kappa     = if (!is.null(tp$kappa)) tp$kappa else 1.0,
      probs     = c(0.5, 0.95)
    )

    # ── Posterior predictive quantile CI (reused by both metrics and WIS) ───

    post_ci <- compute_posterior_predictive_quantile_ci(
      fit           = fit,
      dist_type     = sc$dist_type,
      n_datasets    = sc$n_datasets,
      probs         = WIS_ALL_PROBS,
      n_post_draws  = 500L,
      n_mc_per_draw = 1000L
    )

    # ── WIS (pred_q and true_median passed in to avoid redundant MC draws) ──

    wis_result <- compute_wis(
      fit          = fit,
      dist_type    = sc$dist_type,
      n_datasets   = sc$n_datasets,
      true_params  = tp,
      alpha_levels = WIS_ALPHA_LEVELS,
      pred_q       = post_ci$median,
      true_median  = true_qs["50%"]
    )

    bk <- if (has_kappa) compute_median_bias(fit, "kappa", tp$kappa) else NA_real_

    mutate(na_row,
      max_rhat          = .env$max_rhat,
      min_neff          = .env$min_neff,
      converged         = (.env$max_rhat <= 1.05) & (.env$min_neff >= 100),
      # Hyperparameter
      coverage_mu0      = check_coverage(fit, "mu0", tp$mu0),
      bias_mu0          = compute_median_bias(fit, "mu0", tp$mu0),
      coverage_tau      = check_coverage(fit, "tau", tp$tau),
      bias_tau          = compute_median_bias(fit, "tau", tp$tau),
      coverage_phi      = check_coverage(fit, "phi", tp$phi),
      bias_phi          = compute_median_bias(fit, "phi", tp$phi),
      coverage_kappa    = if (has_kappa) check_coverage(fit, "kappa", tp$kappa)
                         else NA_real_,
      bias_kappa        = bk,
      relbias_kappa     = if (has_kappa) bk / tp$kappa else NA_real_,
      # Predictive quantile (median)
      coverage_pred_median = true_qs["50%"] >= post_ci$lower[IDX_MEDIAN] &
                             true_qs["50%"] <= post_ci$upper[IDX_MEDIAN],
      bias_pred_median     = post_ci$median[IDX_MEDIAN] - true_qs["50%"],
      # Predictive quantile (95th percentile)
      coverage_pred_q95 = true_qs["95%"] >= post_ci$lower[IDX_Q95] &
                          true_qs["95%"] <= post_ci$upper[IDX_Q95],
      bias_pred_q95     = post_ci$median[IDX_Q95] - true_qs["95%"],
      # Scoring
      iqd     = compute_iqd(fit, tp, sc$dist_type),
      wis     = wis_result$wis,
      rel_wis = wis_result$rel_wis
    )

  }, error = function(e) {
    warning(sprintf("[scenario %s | sim %d] %s", sc$scenario_name, sim_idx, e$message))
    na_row
  })
}

# Run all n_sim replicates for one scenario, saving incrementally.
run_scenario <- function(sc, stan_model, n_sim, seed, results_dir, force_rerun) {
  out_file <- file.path(results_dir,
                        paste0(gsub("[^A-Za-z0-9_.-]", "_", sc$scenario_name), ".rds"))
  do_force <- isTRUE(force_rerun) ||
              (is.character(force_rerun) && sc$scenario_name %in% force_rerun)

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
    note <- if (!is.na(sc$scenario_feasibility_note))
              paste0(" [note: ", sc$scenario_feasibility_note, "]") else ""
    cat(sprintf("           dist=%s D=%d N=%.0f p=[%.2f %.2f %.2f %.2f] phi=%.2f kappa=%.2f%s\n",
                sc$dist_type, sc$n_datasets, sc$n_obs_mean,
                sc$summary_type_1_prop, sc$summary_type_2_prop,
                sc$summary_type_3_prop, sc$summary_type_4_prop,
                sc$phi, sc$kappa, note))
  }

  accum <- if (nrow(prev) > 0) list(prev) else list()

  for (sim in (done + 1L):n_sim) {
    cat(sprintf("    sim %d/%d\n", sim, n_sim))
    row   <- run_one_sim(sc, sim, stan_model, seed)
    accum <- c(accum, list(row))
    saveRDS(bind_rows(accum), out_file)
  }

  invisible(NULL)
}

# Merge all per-scenario files into a single tibble.
#
# Completeness rule: a scenario file is complete when it contains exactly
# n_sim rows — one per replicate — where every row is accounted for by one of:
#   (a) a fitted replicate  (MCMC ran; skipped_reason is NA, metrics are present)
#   (b) an explicit skip    (skipped_reason is set by a heuristic; metrics are NA)
#   (c) an error replicate  (all metrics NA, no skipped_reason; replicate was
#                            attempted but the model threw an error)
#
# run_one_sim() always writes exactly one row per replicate regardless of
# outcome, so nrow(d) == number of replicates attempted.  Using nrow() rather
# than max(sim) guards against files that have a high max(sim) but missing
# intermediate rows (e.g. from a mid-write corruption).
#
# A scenario where every replicate was skipped by a heuristic (all rows have
# skipped_reason set) is treated as complete: 100/100 replicates are
# accounted for, just none were fitted.
collect_results <- function(results_dir, n_sim = N_SIM) {
  files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    message("No result files found in ", results_dir)
    return(tibble())
  }
  data_list <- purrr::map(files, readRDS)
  names(data_list) <- basename(files)

  # For each file, count accounted-for replicates: fitted + skipped + errored.
  n_done <- purrr::map_int(data_list, function(d) as.integer(nrow(d)))

  is_complete  <- n_done >= n_sim
  n_incomplete <- sum(!is_complete)

  if (n_incomplete > 0L) {
    message(sprintf(
      "%d incomplete scenario file(s) excluded (%d accounted-for replicates required):",
      n_incomplete, n_sim
    ))
    for (i in which(!is_complete))
      message(sprintf("  [%d/%d rows]  %s", n_done[i], n_sim, names(data_list)[i]))
  }

  message(sprintf("Loaded %d complete / %d total scenario files.",
                  sum(is_complete), length(files)))
  purrr::map_dfr(data_list[is_complete], identity)
}


# ── 5. Execution ──────────────────────────────────────────────────────────────

cat("\n========================================================\n")
cat("Running extended simulation study\n")
cat(sprintf("  Scenarios   : %d\n", nrow(all_scenarios)))
cat(sprintf("  Sims/scen   : %d\n", N_SIM))
cat(sprintf("  Results dir : %s\n", RESULTS_DIR))
cat("  New metrics : coverage/bias of pred. median and 95th pct, WIS\n")
cat("  Heuristics  : GG (should_attempt_gg), Gamma+type-2 (gamma_type2_reliable)\n")
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

out_combined <- file.path(MAIN_FOLDER, "results", "new_simulation_results_all.rds")
saveRDS(results_all, out_combined)
write.csv(results_all, sub("\\.rds$", ".csv", out_combined), row.names = FALSE)
cat(sprintf("Saved %d rows to %s\n", nrow(results_all), out_combined))


# ── 7. Summary ────────────────────────────────────────────────────────────────

# 7a. Expected vs actual skip rates
# scenario_feasibility_note = what we predicted before running
# skipped_reason            = what actually happened per replicate
cat("\n--- Skip summary (expected note × actual reason) ---\n")
results_all |>
  mutate(
    note   = tidyr::replace_na(scenario_feasibility_note, "(none)"),
    reason = tidyr::replace_na(skipped_reason,            "(fitted)")
  ) |>
  count(note, reason) |>
  print()

# 7b. Performance summary for fitted replicates
cat("\n--- Convergence and coverage by scenario group (fitted replicates only) ---\n")
results_all |>
  filter(is.na(skipped_reason)) |>
  group_by(scenario_group) |>
  summarise(
    n_sims           = n(),
    pct_conv         = round(mean(converged,            na.rm = TRUE) * 100, 1),
    med_rhat         = round(median(max_rhat,           na.rm = TRUE), 3),
    cov_mu0          = round(mean(coverage_mu0,         na.rm = TRUE) * 100, 1),
    cov_phi          = round(mean(coverage_phi,         na.rm = TRUE) * 100, 1),
    cov_kappa        = round(mean(coverage_kappa,       na.rm = TRUE) * 100, 1),
    cov_pred_median  = round(mean(coverage_pred_median, na.rm = TRUE) * 100, 1),
    cov_pred_q95     = round(mean(coverage_pred_q95,    na.rm = TRUE) * 100, 1),
    med_wis          = round(median(wis,                na.rm = TRUE), 4),
    med_rel_wis      = round(median(rel_wis,            na.rm = TRUE), 4)
  ) |>
  print()
