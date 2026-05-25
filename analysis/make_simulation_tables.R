# analysis/make_simulation_tables.R
# =============================================================================
# Generates LaTeX supplementary tables for the simulation study.
#
# Table S1: Performance by distribution family x summary type
# Table S2: Performance by data availability (N datasets x n obs x summary type)
#
# Usage:
#   Rscript analysis/make_simulation_tables.R
#   -- or --
#   source("analysis/make_simulation_tables.R")
#
# Outputs written to results/
#   simulation_table1.tex  (Table Supplementary Information)
#   simulation_table2.tex  (Table Supplementary Information)
# =============================================================================

devtools::load_all(here::here(), quiet = TRUE)

# ── Load simulation results ────────────────────────────────────────────────────

rds_path <- system.file("extdata", "simulation_results_all.rds", package = "ddsynth")
if (!nzchar(rds_path))
  rds_path <- file.path(here::here(), "vignettes", "simulation_results_all.rds")

if (!file.exists(rds_path))
  stop("simulation_results_all.rds not found. Run analysis/new_simulation_study.R first.")

res_out <- readRDS(rds_path)
cat(sprintf("Loaded %d rows across %d scenarios.\n",
            nrow(res_out), length(unique(res_out$scenario_name))))

# ── Summarise ──────────────────────────────────────────────────────────────────

summary_res <- create_results_summary(res_out)
cat(sprintf("Summary: %d scenario-level rows.\n", nrow(summary_res)))

# ── Generate and write tables ─────────────────────────────────────────────────

out_dir <- file.path(here::here(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

tbl1 <- generate_simulation_table1(summary_res)
tbl2 <- generate_simulation_table2(summary_res)

writeLines(tbl1, file.path(out_dir, "simulation_table1.tex"))
writeLines(tbl2, file.path(out_dir, "simulation_table2.tex"))

cat("Written: results/simulation_table1.tex\n")
cat("Written: results/simulation_table2.tex\n")
