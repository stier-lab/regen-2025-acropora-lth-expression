# =============================================================================
# Purpose: Run the full LTH analysis pipeline in dependency order.
#          Source this once from the project root: `Rscript code/_run_all.R`
#          Total runtime ~3 minutes (Apex XML parse dominates first pass).
#
# What & why: the MASTER pipeline runner — the single entry point that
#   reproduces every table and figure from the raw data in one command. It
#   source()s each analysis script in turn. The ONLY thing that matters here is
#   the ORDER: many scripts read .rds/.csv files that an EARLIER script wrote, so
#   running them out of order would read stale (or missing) inputs. The list
#   below is therefore a dependency order, not an alphanumeric one. A few
#   examples (see also the inline notes):
#     - 01 (load/clean metadata) runs first: nearly everything downstream reads
#       the cleaned data it produces.
#     - 07 (wax dipping) runs BEFORE 05 (buoyant weight) because 05's areal
#       calcification needs the per-coral surface areas that 07 writes — despite
#       "05" sorting before "07" numerically.
#     - 31 (RNA-seq covariate handoff) runs late: it joins the symbiont table
#       (06) and the genet resilience summary (19), so both must exist first.
#     - 32 (preliminary SNP integration) runs after 31: it joins Rachael's
#       provisional cluster file to the RNA-seq covariates and phenotype
#       summaries. It is exploratory and exits cleanly if the file is absent.
#     - 33-35 (trade-off screens) run late: they reuse phenotype summary tables
#       to ask whether source C's lower heat sensitivity appears to come with
#       growth or regeneration costs.
#     - 30 (manuscript audit) runs DEAD LAST: it checks the manuscript against
#       the freshly regenerated tables, so every table must already be rebuilt.
#   Edit this list with care: reordering can silently feed a script stale
#   outputs.
# =============================================================================

# Source from project root (resolve via here::here, which uses .Rproj/.here).
# setwd() so the relative file.path("code", s) calls below resolve regardless of
# where Rscript was launched.
project_root <- here::here()
setwd(project_root)

# Pipeline steps, in dependency order. Read top-to-bottom = the order they run.
scripts <- c(
  "01_load_clean_metadata.R",
  "02_pam_analysis.R",
  "03_color_card_analysis.R",
  "04_physio_morphology.R",
  "07_wax_dipping.R",          # must run before 05 — 05 reads wax_clean.rds (surface areas)
  "05_buoyant_weight.R",
  "06_symbiont_chl.R",
  "08_apex_temperature.R",
  "09_ysi_water_chem.R",
  "10_worms.R",
  "11_microscope_physio.R",    # separate photo-only microscope cohort; not pooled with main morphology
  "11c_morphology_dataset_plots.R", # comprehensive morphology plots across both separate cohorts
  "12_models.R",               # primary models + color-CLMM & morphology-blme robustness (was 12/12b/12c)
  "13_genet_interaction.R",
  "14_morphology_kaplan.R",
  "15_multivariate.R",
  "16_main_figure.R",
  "17_figure_audit.R",
  "18_data_validation.R",
  "19_genet_dashboard.R",
  "sensitivity/22_sensitivity_flagged.R",
  "sensitivity/23_timeseries_diagnostics.R",
  "sensitivity/24_headline_model_comparison.R",
  "sensitivity/26_thermal_context.R",
  "sensitivity/27_variance_partitioning.R",
  "sensitivity/28_multiple_testing.R",
  "sensitivity/29_morphology_prob_contrasts.R",
  "diagnostics/A_continuous_lmm.R",
  "diagnostics/B_morphology_glmm.R",
  "diagnostics/C_cox_diagnostics.R",
  "diagnostics/D_pca_lrt.R",
  "diagnostics/E_design_alignment.R",
  "diagnostics/F_model_reproducibility.R",
  "diagnostics/G_diagnostic_plots.R",
  "20_master_results_table.R",
  "diagnostics/H_spreadsheet_coverage.R",
  "31_rnaseq_covariate_table.R",    # handoff: per-library phenotype covariates for the RNA-seq
  "32_prelim_snp_phenotype_integration.R", # exploratory: preliminary SNP clusters x RNA-seq covariates
  "33_growth_tradeoff_screen.R",    # exploratory: ambient growth vs growth lost under heat
  "34_regeneration_tradeoff_screen.R", # exploratory: growth/regeneration and heat-tolerance/regeneration screens
  "35_tradeoff_summary.R",          # team-summary four-panel trade-off check
  "sensitivity/25_model_diagnostic_coverage.R",
  "30_manuscript_audit.R"          # advisory phenotype reproducibility check — warns (never fails) if phenotype numbers drift
)

# Run each script in turn, printing a banner and per-script timing. A failure in
# any script throws and stops the loop here (so you see exactly which step broke)
# — except 30's audit, which only WARNS and lets the run finish.
t0 <- Sys.time()
for (s in scripts) {
  cat("\n============================================================\n",
      "Running:", s, "\n",
      "============================================================\n", sep = "")
  tic <- Sys.time()
  source(file.path("code", s))
  toc <- Sys.time()
  cat("  -> done in", round(as.numeric(toc - tic, units = "secs"), 1), "s\n")
}
t1 <- Sys.time()
cat("\n\nFull pipeline complete in",
    round(as.numeric(t1 - t0, units = "mins"), 2), "min.\n")

# Save session info (R version + exact package versions) for reproducibility
# provenance — so a future reader knows the software state this run was built on.
dir.create("output", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
