# =============================================================================
# Purpose: Parent inventory of every analysis family and its diagnostic coverage.
#
# What & why: this is the repo-level answer to "have we checked diagnostics and
#   goodness of fit for all the analyses?" It groups repeated model fits into
#   analysis families, records the scripts and outputs that support each family,
#   and states whether the relevant quality check is model goodness-of-fit,
#   data/provenance checking, exact-test estimability, or "not applicable".
#   The table is intentionally coarser than output/tables/20_master_results.csv:
#   one row per analysis family, not one row per coefficient or contrast.
# Input:   output/tables/*.csv, output/diagnostics/*.csv, figures/*
# Output:  output/tables/38_analysis_diagnostic_parent_list.csv
#          output/tables/38_analysis_diagnostic_inventory.csv
#          output/tables/38_diagnostic_artifact_status.csv
#          docs/analysis_diagnostic_inventory.md
# =============================================================================

source(here::here("code", "00_setup.R"))

DIAG_OUT <- file.path(OUT_DIR, "diagnostics")
DOC_PATH <- here::here("docs", "analysis_diagnostic_inventory.md")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(DOC_PATH), recursive = TRUE, showWarnings = FALSE)

rel <- function(...) file.path(...)

path_exists <- function(path) {
  expanded <- Sys.glob(here::here(path))
  if (length(expanded) == 0) {
    file.exists(here::here(path))
  } else {
    any(file.exists(expanded))
  }
}

split_paths <- function(x) {
  if (is.na(x) || !nzchar(x) || x %in% c("none", "n/a", "not applicable")) {
    return(character())
  }
  str_split(x, "\\s*;\\s*")[[1]] |>
    str_trim() |>
    purrr::discard(\(z) !nzchar(z) || z %in% c("none", "n/a", "not applicable"))
}

status_counts <- function(path, status_candidates = c("status", "test_status",
                                                       "conclusion")) {
  full <- here::here(path)
  if (!file.exists(full)) return("missing")
  dat <- suppressMessages(readr::read_csv(full, show_col_types = FALSE))
  status_col <- intersect(status_candidates, names(dat))[1]
  if (is.na(status_col)) return(sprintf("present, %d rows", nrow(dat)))
  dat |>
    count(.data[[status_col]], name = "n") |>
    arrange(.data[[status_col]]) |>
    mutate(part = paste0(.data[[status_col]], "=", n)) |>
    pull(part) |>
    paste(collapse = "; ")
}

coverage_summary <- function() {
  p <- rel("output", "tables", "25_model_diagnostic_coverage.csv")
  if (!file.exists(here::here(p))) return("model coverage table missing")
  x <- suppressMessages(readr::read_csv(here::here(p), show_col_types = FALSE))
  n_gap <- sum(str_detect(str_to_lower(x$status), "gap|missing|failed|fail"))
  sprintf("Script 25 snapshot: %d/%d fitted models covered; %d coverage gaps (new script 42 fits audited separately)",
          nrow(x) - n_gap, nrow(x), n_gap)
}

fail_count <- function(path, fail_levels = c("FAIL", "failed", "gap", "missing")) {
  full <- here::here(path)
  if (!file.exists(full)) return(NA_integer_)
  dat <- suppressMessages(readr::read_csv(full, show_col_types = FALSE))
  status_col <- intersect(c("status", "test_status", "conclusion"), names(dat))[1]
  if (is.na(status_col)) return(0L)
  sum(str_detect(str_to_lower(dat[[status_col]]),
                 paste(str_to_lower(fail_levels), collapse = "|")),
      na.rm = TRUE)
}

diag_counts <- list(
  validation = status_counts(rel("output", "tables", "18_validation_summary.csv")),
  cross_dataset_morph = status_counts(rel("output", "diagnostics",
                                          "11f_morphology_dataset_diagnostic_checks.csv")),
  continuous = status_counts(rel("output", "diagnostics",
                                 "A_continuous_diagnostics.csv")),
  morph_glmm = status_counts(rel("output", "diagnostics",
                                 "B_morphology_glmm_diagnostics.csv")),
  cox = status_counts(rel("output", "diagnostics", "C_cox_diagnostics.csv")),
  pca_lrt = status_counts(rel("output", "diagnostics",
                              "D_pca_lrt_diagnostics.csv")),
  design = status_counts(rel("output", "diagnostics", "E_design_alignment.csv")),
  reproducibility = status_counts(rel("output", "diagnostics",
                                      "F_model_reproducibility.csv")),
  plot_inventory = status_counts(rel("output", "diagnostics",
                                     "G_plot_inventory.csv")),
  spreadsheet = status_counts(rel("output", "diagnostics",
                                  "H_coverage_by_source.csv")),
  timeseries = "9 checks: AR(1), random-slope, nonlinear-time, and headline-effect robustness; reported heat effects unchanged",
  multiple = {
    p <- rel("output", "tables", "28_multiple_testing.csv")
    if (file.exists(here::here(p))) {
      x <- suppressMessages(readr::read_csv(here::here(p), show_col_types = FALSE))
      sprintf("%d tests; %d reported significant after confirmatory/raw or exploratory/BH rule",
              nrow(x), sum(x$significant, na.rm = TRUE))
    } else {
      "multiple-testing table missing"
    }
  },
  tradeoff = status_counts(rel("output", "tables",
                               "37_tradeoff_model_diagnostics.csv"))
)

model_coverage <- coverage_summary()

parent_list <- tibble::tribble(
  ~parent_id, ~parent_analysis, ~plain_question, ~analysis_ids,
  "P1", "Data, design, and environment",
  "Did the experiment and environmental data match the planned design?",
  "data_integrity; temperature_context",
  "P2", "Whole-fragment physiology",
  "Did heat change photosynthesis score, color, growth, or symbionts?",
  "pam_fvfm; color_score; skeletal_growth; symbiont_density",
  "P3", "Wound healing and regeneration",
  "Did heat change wound closure, tip regrowth, and new corallites?",
  "main_morphology_states; morphology_timing; microscope_validation",
  "P4", "Source-patch response",
  "Did source patches A, C, and D respond differently?",
  "source_interactions; multivariate_pca; source_resilience_composite",
  "P5", "External genotype and RNA-seq bridge",
  "What can the preliminary DNA-marker and external Hauru data tell us without overclaiming identity?",
  "preliminary_snp; trinity_context; rnaseq_covariates",
  "P6", "Trade-off tests",
  "Is resilience associated with a cost in growth or regeneration?",
  "tradeoff_formal; tradeoff_screens",
  "P7", "Repo-wide reproducibility",
  "Do figures, tables, model diagnostics, and manuscript numbers stay in sync?",
  "multiple_testing; model_diagnostic_coverage; figure_and_manuscript_audits"
) |>
  mutate(
    diagnostic_summary = case_when(
      parent_id == "P1" ~ paste("Validation:", diag_counts$validation),
      parent_id == "P2" ~ paste("Continuous models:", diag_counts$continuous),
      parent_id == "P3" ~ paste("Morphology GLMM:", diag_counts$morph_glmm,
                                "| Cox:", diag_counts$cox),
      parent_id == "P4" ~ paste("PCA/LRT:", diag_counts$pca_lrt),
      parent_id == "P5" ~ "Join/design checks done; SNP-symbiont model diagnostics are partial because the SNP file is preliminary and source-patch identity is unresolved.",
      parent_id == "P6" ~ paste("Trade-off models:", diag_counts$tradeoff),
      TRUE ~ paste(model_coverage, "| Figure inventory:", diag_counts$plot_inventory)
    )
  )

inventory <- tibble::tribble(
  ~parent_id, ~analysis_id, ~analysis_family, ~role, ~scripts, ~primary_outputs,
  ~display_items, ~model_or_test, ~diagnostics_or_qc, ~diagnostics_done,
  ~goodness_of_fit_done, ~diagnostic_outcome, ~remaining_gap, ~manuscript_use,

  "P1", "data_integrity", "Cleaned data and experimental design",
  "support",
  "code/01_load_clean_metadata.R; code/18_data_validation.R; code/diagnostics/H_spreadsheet_coverage.R",
  "output/tables/01_metadata_summary.csv; output/tables/18_validation_summary.csv; output/diagnostics/H_coverage_by_source.csv",
  "none",
  "Data validation, ID checks, tank-treatment checks, and spreadsheet coverage",
  "Sample-size checks, factor-level checks, orphan-ID checks, tank assignment checks, duplicate checks",
  "yes", "not applicable",
  paste("Data validation:", diag_counts$validation),
  "No model fit here; keep this as the first stop-gate before biological interpretation.",
  "Use as provenance and design QC, not as biological inference.",

  "P1", "temperature_context", "Apex logger, YSI, and heat-treatment context",
  "support",
  "code/08_apex_temperature.R; code/09_ysi_water_chem.R; code/36_molly_followup_checks.R",
  "output/tables/36_molly_temperature_coverage_days0_16.csv",
  "figures/08_apex_temperature.pdf; figures/08b_apex_temperature_full.pdf; figures/09_ysi_water_chem.pdf",
  "Descriptive environmental summaries",
  "Coverage audit for days 0-16; tank-treatment map checks through the shared setup",
  "yes", "not applicable",
  "Recovered Apex exports cover all eight tanks every ten minutes through Days 0-15; compiled sheet matches XML after duplicate checks. YSI coverage remains separate.",
  "Confirm the logger clock setting before exact clock-hour exposure calculations; see the recovery provenance note.",
  "Daily tank means differ by 3.09 C through Day 15. Post-experiment records are excluded from treatment means.",

  "P2", "pam_fvfm", "Photosynthesis-efficiency score (PAM Fv/Fm)",
  "primary",
  "code/02_pam_analysis.R; code/12_models.R",
  "output/tables/02_pam_treatment_contrasts.csv; output/tables/12_anova_summary.csv; output/tables/12_emmeans_contrasts.csv",
  "figures/02_pam_fvfm_trajectory.pdf",
  "Linear mixed model with treatment x wound x day x source patch plus tank and coral ID random intercepts",
  "DHARMa residuals, convergence/singularity, influence checks, time-series AR(1)/random-slope/quadratic sensitivity, headline model comparison",
  "yes", "yes",
  paste(diag_counts$continuous, "|", diag_counts$timeseries),
  "One Cook's-distance warning remains, but the heat direction check and upgraded-model comparison preserve the conclusion.",
  "Primary whole-fragment physiology result.",

  "P2", "pam_location", "Paired near-tip versus near-base photosynthesis",
  "support",
  "code/42_pam_location_and_score_audit.R",
  "output/tables/42_pam_location_heat_tests.csv; output/tables/42_pam_location_diagnostics.csv; output/tables/42_pam_location_trajectory_test.csv",
  "figures/42_pam_location_comparison.pdf; figures/42_pam_location_diagnostics.pdf",
  "Tank-level exact label permutations; paired-difference repeated-measures LMM; correlated unequal-variance sensitivity",
  "Complete-pair and canonical-average agreement; residual/Q-Q inspection; convergence; singularity; residual serial correlation; leave-one-tank-out; tank medians",
  "yes", "partial",
  "LMM converged and is not singular, but residual tails and unequal spread remain. Tank permutations give mean heat-gap p=0.20 and trajectory p=0.143; correlated model gives mean-gap p=0.204.",
  "Only eight tanks. Do not interpret the LMM-only treatment-by-day p=0.030 as decisive; no prespecified equivalence margin. See September 22 follow-up note.",
  "Both sites retain the heat decline; averaging is a two-location summary, not proof of equivalence.",

  "P1", "morphology_score_review", "Tip reversals and missing pigment scores",
  "support",
  "code/42_pam_location_and_score_audit.R",
  "output/tables/42_morphology_fragment_review.csv; output/tables/42_morphology_missing_calls.csv",
  "none",
  "Descriptive score audit, no imputation or model refit",
  "Unique ID/day keys; all wounded-fragment tip sequences and missing calls",
  "yes", "not applicable",
  "Tip reversals in 116, 121, 143; 23/24 Day 9 pigment calls missing. Original observations preserved.",
  "Field-note review pending. Separate microscope cohort cannot validate individual main-experiment scores.",
  "Carry score uncertainty in captions and first-observation interpretation.",

  "P2", "color_score", "Color-card paling score",
  "primary",
  "code/03_color_card_analysis.R; code/12_models.R",
  "output/tables/03_color_end_proportions.csv; output/tables/12_anova_summary.csv; output/tables/12b_color_clmm.csv",
  "figures/03_color_trajectory.pdf",
  "Linear mixed model plus ordinal cumulative-link mixed-model robustness check",
  "DHARMa/residual checks for the LMM, observed-vs-fitted CLMM check, influence checks, time-shape sensitivity",
  "yes", "yes",
  paste(diag_counts$continuous, "| CLMM diagnostic plot covered in", model_coverage),
  "Gaussian model on an ordinal scale is explicitly handled by the CLMM robustness model.",
  "Primary whole-fragment physiology result, with ordinal caveat disclosed.",

  "P2", "skeletal_growth", "Buoyant-weight skeletal growth",
  "primary",
  "code/05_buoyant_weight.R; code/12_models.R",
  "output/tables/05_buoyant_weight_lm.csv; output/tables/05_buoyant_weight_tank_test.csv; output/tables/12_anova_summary.csv",
  "figures/05_buoyant_weight_growth.pdf",
  "Endpoint LM/LMM and tank-level randomization check",
  "Residual plots, Shapiro-Wilk, Cook's distance, tank-level randomization, model coverage",
  "yes", "yes",
  paste(diag_counts$continuous, "|", model_coverage),
  "Surface area caveat remains: growth is percent skeletal mass change, not areal calcification.",
  "Primary whole-fragment physiology result.",

  "P2", "symbiont_density", "Symbiont density",
  "primary",
  "code/06_symbiont_chl.R; code/12_models.R",
  "output/tables/06_symbiont_chl_summary.csv; output/tables/12_anova_summary.csv; output/tables/12_genet_treatment_effects.csv",
  "figures/06_symbiont_density_by_day.pdf",
  "Linear mixed model on log symbiont density",
  "DHARMa residuals, outlier sensitivity, influence checks, model coverage",
  "yes", "yes",
  paste(diag_counts$continuous, "|", model_coverage),
  "Some residual/outlier checks are handled by top-residual sensitivity; conclusion stays the same.",
  "Primary whole-fragment physiology result and current direct SNP-linked response.",

  "P3", "main_morphology_states", "Daily wound-healing and regeneration states",
  "primary",
  "code/04_physio_morphology.R; code/12_models.R; code/sensitivity/29_morphology_prob_contrasts.R",
  "output/tables/04_morphology_trait_glmm_summaries.csv; output/tables/12c_morph_blme_anova.csv; output/tables/29_morphology_prob_contrasts.csv",
  "figures/04_morphology_trajectories.pdf; figures/04b_morphology_trajectories_by_genet.pdf",
  "Binomial GLMMs and weakly regularized binomial GLMMs for separated/saturated traits",
  "DHARMa residual screens, convergence/singularity checks, separation-aware penalized fits, probability-scale contrasts",
  "yes", "yes",
  paste(diag_counts$morph_glmm, "|", model_coverage),
  "Several traits saturate or separate; these are handled, and final regeneration inference leans on timing models.",
  "Primary morphology trajectory result.",

  "P3", "morphology_timing", "Time to wound-healing and regeneration milestones",
  "primary",
  "code/14_morphology_kaplan.R; code/diagnostics/C_cox_diagnostics.R",
  "output/tables/14_interval_survreg.csv; output/tables/14_cox_hazard_ratios.csv; output/tables/14_cox_ph_tests.csv",
  "figures/14_morphology_KM.pdf; figures/14b_morphology_KM_by_genet.pdf",
  "Interval-censored Weibull AFT models, Kaplan-Meier curves, Cox proportional-hazards checks",
  "Cox Schoenfeld proportional-hazards tests, event-count checks, interval-censoring sensitivity",
  "yes", "yes",
  paste(diag_counts$cox, "| PH checks have no failures"),
  "Some low-event subsets are handled/skipped; use overall timing models for the manuscript headline.",
  "Primary regeneration timing result.",

  "P3", "microscope_validation", "Separate microscope/photo morphology cohort",
  "support",
  "code/11_microscope_physio.R; code/11c_morphology_dataset_plots.R",
  "output/tables/11_microscope_event_tests.csv; output/tables/11f_morphology_dataset_endpoint_tests.csv; output/tables/11f_morphology_dataset_logrank_tests.csv",
  "figures/11_microscope_trait_trajectories.pdf; figures/11c_morphology_dataset_all_trait_trajectories.pdf; figures/11e_morphology_dataset_event_summary.pdf",
  "Fisher exact endpoint tests and log-rank first-observed timing tests",
  "Scoring coverage, endpoint balance, and test estimability checks",
  "yes", "not applicable",
  paste("Cross-dataset morphology diagnostics:", diag_counts$cross_dataset_morph),
  "Small wounded-only photo cohort; do not pool with the main morphology models.",
  "Supporting visual-validation result.",

  "P4", "source_interactions", "Source-patch dependence of physiology and morphology",
  "primary/support",
  "code/13_genet_interaction.R; code/14_morphology_kaplan.R; code/diagnostics/D_pca_lrt.R",
  "output/tables/13_genet_anova.csv; output/tables/13_genet_emmeans.csv; output/tables/14_cox_genet_LRT.csv",
  "figures/13_genet_response_panel.pdf; figures/14b_morphology_KM_by_genet.pdf",
  "Likelihood-ratio tests comparing models with and without source-patch terms",
  "Design-alignment audit, model reproducibility, generated DHARMa plots, PCA/LRT reproduction checks",
  "yes", "yes",
  paste(diag_counts$pca_lrt, "|", diag_counts$reproducibility),
  "Source patch is a fixed label with only three levels; final genetic identity is pending.",
  "Supports source-patch wording and the C-resilience pattern.",

  "P4", "multivariate_pca", "End-of-experiment multivariate physiology",
  "support/synthesis",
  "code/15_multivariate.R; code/diagnostics/D_pca_lrt.R",
  "output/tables/15_pca_loadings.csv; output/tables/15_genet_pca_displacement.csv",
  "figures/15_physio_PCA_biplot.pdf; figures/15b_physio_PCA_by_genet.pdf",
  "PCA on centered and scaled endpoint physiology variables",
  "Sample-size, complete-case, centering/scaling, variance-explained, loading, and interpretation checks",
  "yes", "yes",
  paste("PCA/LRT diagnostics:", diag_counts$pca_lrt),
  "PCA is a synthesis, not a causal model; missing values are handled by complete-case filtering.",
  "Supports the multivariate source-patch displacement claim.",

  "P4", "source_resilience_composite", "Source-patch heat-sensitivity composite",
  "synthesis/exploratory",
  "code/19_genet_dashboard.R; code/37_tradeoff_formal_tests.R",
  "output/tables/19_genet_resilience_summary.csv; output/tables/19c_resilience_decomp_by_scope.csv; output/tables/37_tradeoff_heat_penalty_leave_one_trait_out.csv",
  "figures/19_genet_dashboard.pdf; figures/19b_genet_resilience_ranking.pdf; figures/19c_decomposed_resilience.pdf",
  "Derived row-max standardized heat-penalty score",
  "Input decomposition and leave-one-trait-out sensitivity; no residual goodness-of-fit because this is a derived index",
  "yes", "not applicable",
  "Leave-one-trait-out scores keep source patch C as the lowest heat-penalty source.",
  "Do not present as a formal hypothesis test; use as a ranking/synthesis.",
  "Useful plain-language synthesis figure.",

  "P5", "rnaseq_covariates", "RNA-seq phenotype covariate handoff",
  "support",
  "code/31_rnaseq_covariate_table.R",
  "output/tables/31_rnaseq_phenotype_covariates.csv; output/tables/31_rnaseq_library_lookup_raw.csv",
  "none",
  "Sample-level covariate table, no fitted response model",
  "Completeness and library lookup checks printed by the script",
  "yes", "not applicable",
  "All 144 libraries have a symbiont-density value in the current handoff.",
  "Expression counts are not yet included; this is a handoff table.",
  "Defines the RNA-seq analysis spine.",

  "P5", "preliminary_snp", "Preliminary SNP cluster and direct SNP-symbiont check",
  "exploratory",
  "code/32_prelim_snp_phenotype_integration.R",
  "output/tables/32_prelim_snp_join_audit.csv; output/tables/32_prelim_snp_response_joinability.csv; output/tables/32_prelim_snp_symbiont_model_summary.csv; output/tables/32_prelim_snp_symbiont_nested_tests.csv",
  "figures/32_prelim_snp_structure.pdf; figures/32_prelim_snp_design_balance.pdf; figures/32_prelim_snp_symbiont_heat_effects.pdf",
  "Chi-square design-balance tests and exploratory LM on same-fragment symbiont density",
  "Join audit, design-balance checks, and coverage flag; residual diagnostics for exploratory symbiont LM not yet added",
  "partial", "partial",
  paste("SNP join audit:", status_counts(rel("output", "tables",
                                             "32_prelim_snp_join_audit.csv"))),
  "Before any SNP-phenotype inference is promoted, add residual/influence diagnostics and rerun with final SNP calls, genetic PCs or kinship, and any inversion/clone-structure flags.",
  "Preliminary context only; keep A/C/D as source-patch labels, not confirmed genets.",

  "P5", "trinity_context", "Trinity/Ross Hauru coordinate and ED50 context",
  "support/external context",
  "code/sensitivity/26b_trinity_hauru_context.R",
  "output/tables/26b_trinity_hauru_table_summary.csv; output/tables/26b_trinity_hauru_nearest_source_patch_matches.csv",
  "figures/26b_trinity_hauru_context.pdf",
  "Nearest-coordinate descriptive screen, no fitted model",
  "File provenance, distance/candidate-match checks, clonal-group context, and ED50 range check",
  "yes", "not applicable",
  "C is closest to Apul-115 by coordinate, A and D are closest to Apul-111, and nearby colonies share clonal group 2; no source patch is assigned an external ED50.",
  "Needs DNA-marker comparison before assigning Trinity/Ross colony IDs or ED50 values to A/C/D; Molly's forwarded context also makes clear that clones can span thickets, so distance alone is not enough.",
  "Context only; useful for prioritizing the Apul-115/Apul-111 genetic comparison.",

  "P6", "tradeoff_screens", "Exploratory growth, heat, and regeneration trade-off screens",
  "exploratory",
  "code/33_growth_tradeoff_screen.R; code/34_regeneration_tradeoff_screen.R; code/35_tradeoff_summary.R",
  "output/tables/33_growth_heat_tradeoff.csv; output/tables/34_heat_tolerance_regeneration_screen.csv; output/tables/35_tradeoff_summary.csv",
  "figures/33_growth_heat_tradeoff.pdf; figures/34_regeneration_tradeoff_screens.pdf; figures/35_tradeoff_summary.pdf",
  "Descriptive source-patch and fragment-level screens",
  "Join audits for growth/morphology merges; formal checks now live in script 37",
  "yes", "not applicable",
  "Exploratory screens are paired with formal endpoint, robust, timing, and leave-one-out checks.",
  "Do not infer source-level slopes strongly from only three source patches.",
  "Useful for interpretation, not a standalone formal result.",

  "P6", "tradeoff_formal", "Formal trade-off model set",
  "support/sensitivity",
  "code/37_tradeoff_formal_tests.R",
  "output/tables/37_tradeoff_main_results.csv; output/tables/37_tradeoff_model_diagnostics.csv; output/tables/37_tradeoff_heat_penalty_leave_one_trait_out.csv",
  "figures/37_tradeoff_diagnostics/*.png; figures/35_tradeoff_summary.pdf",
  "Endpoint LM, robust LM, bias-reduced logistic models, Cox timing, and discrete-time transition models",
  "Residual, influence, robust-weight, DHARMa, calibration, leverage, and Cox proportional-hazards checks",
  "yes", "yes",
  paste("Trade-off diagnostics:", diag_counts$tradeoff),
  "Endpoint models are small; robust and daily transition checks keep the same no-cost conclusion.",
  "Supports the current no-strong-regeneration-cost wording.",

  "P7", "multiple_testing", "Confirmatory versus exploratory p-value handling",
  "diagnostic",
  "code/sensitivity/28_multiple_testing.R",
  "output/tables/28_multiple_testing.csv",
  "none",
  "Benjamini-Hochberg false-discovery-rate correction for exploratory family only",
  "Rules encoded in script and output table",
  "yes", "not applicable",
  paste("Multiple-testing table:", diag_counts$multiple),
  "Keep confirmatory/exploratory labels synchronized with manuscript wording.",
  "Supports statistical-reporting discipline.",

  "P7", "model_diagnostic_coverage", "Every saved model has a diagnostic artifact",
  "diagnostic",
  "code/sensitivity/25_model_diagnostic_coverage.R; code/diagnostics/G_diagnostic_plots.R",
  "output/tables/25_model_diagnostic_coverage.csv; output/diagnostics/K_model_coverage_report.md; output/diagnostics/G_plot_inventory.csv",
  "figures/diagnostics/*.png; figures/37_tradeoff_diagnostics/*.png",
  "Model inventory coverage check",
  "Checks each saved model against residual, PH, ordinal, or trade-off diagnostic artifacts",
  "yes", "yes",
  paste(model_coverage, "| plot inventory:", diag_counts$plot_inventory),
  "Coverage says a diagnostic exists; the family-specific rows above say whether that diagnostic passed cleanly or had handled caveats.",
  "Repo-wide diagnostic gate.",

  "P7", "figure_and_manuscript_audits", "Figure freshness and manuscript-number checks",
  "diagnostic",
  "code/17_figure_audit.R; code/30_manuscript_audit.R",
  "output/tables/17_figure_audit.csv; output/tables/30_manuscript_audit.csv",
  "figures/FIGURE_INDEX.md",
  "File freshness checks and regenerated-number audit",
  "Figure existence/freshness and manuscript-number cross-checks",
  "yes", "not applicable",
  paste("Figure audit:", status_counts(rel("output", "tables",
                                           "17_figure_audit.csv")),
        "| Manuscript audit:",
        status_counts(rel("output", "tables", "30_manuscript_audit.csv"),
                      status_candidates = c("found"))),
  "This catches stale outputs, but it does not evaluate statistical model fit.",
  "Prevents prose/figure/table drift."
) |>
  left_join(parent_list |> select(parent_id, parent_analysis), by = "parent_id") |>
  relocate(parent_analysis, .after = parent_id)

artifact_files <- tibble::tribble(
  ~artifact, ~scope, ~status_counts,
  rel("output", "tables", "18_validation_summary.csv"),
  "data-integrity checks", diag_counts$validation,
  rel("output", "diagnostics", "A_continuous_diagnostics.csv"),
  "continuous-response mixed-model diagnostics", diag_counts$continuous,
  rel("output", "diagnostics", "B_morphology_glmm_diagnostics.csv"),
  "morphology GLMM diagnostics", diag_counts$morph_glmm,
  rel("output", "diagnostics", "C_cox_diagnostics.csv"),
  "Cox proportional-hazards diagnostics", diag_counts$cox,
  rel("output", "diagnostics", "D_pca_lrt_diagnostics.csv"),
  "PCA and source-patch LRT diagnostics", diag_counts$pca_lrt,
  rel("output", "diagnostics", "E_design_alignment.csv"),
  "model/design structure audit", diag_counts$design,
  rel("output", "diagnostics", "F_model_reproducibility.csv"),
  "model refit reproducibility", diag_counts$reproducibility,
  rel("output", "diagnostics", "G_plot_inventory.csv"),
  "diagnostic plot inventory", diag_counts$plot_inventory,
  rel("output", "diagnostics", "11f_morphology_dataset_diagnostic_checks.csv"),
  "cross-dataset morphology support checks", diag_counts$cross_dataset_morph,
  rel("output", "tables", "23_timeseries_diagnostics.csv"),
  "time-series sensitivity", diag_counts$timeseries,
  rel("output", "tables", "25_model_diagnostic_coverage.csv"),
  "all saved model diagnostic coverage", model_coverage,
  rel("output", "tables", "28_multiple_testing.csv"),
  "confirmatory/exploratory p-value handling", diag_counts$multiple,
  rel("output", "tables", "32_prelim_snp_join_audit.csv"),
  "preliminary SNP join audit",
  status_counts(rel("output", "tables", "32_prelim_snp_join_audit.csv")),
  rel("output", "tables", "37_tradeoff_model_diagnostics.csv"),
  "formal trade-off model diagnostics", diag_counts$tradeoff,
  rel("output", "tables", "30_manuscript_audit.csv"),
  "manuscript-number reproducibility",
  status_counts(rel("output", "tables", "30_manuscript_audit.csv"),
                status_candidates = c("found"))
) |>
  mutate(
    exists = map_lgl(artifact, path_exists),
    fail_rows = map_int(artifact, \(p) coalesce(fail_count(p), 0L))
  )

inventory_checked <- inventory |>
  rowwise() |>
  mutate(
    script_paths_checked = paste(split_paths(scripts), collapse = "; "),
    missing_scripts = {
      s <- split_paths(scripts)
      miss <- s[!map_lgl(s, path_exists)]
      if (length(miss)) paste(miss, collapse = "; ") else "none"
    },
    missing_primary_outputs = {
      p <- split_paths(primary_outputs)
      miss <- p[!map_lgl(p, path_exists)]
      if (length(miss)) paste(miss, collapse = "; ") else "none"
    }
  ) |>
  ungroup()

parent_checked <- parent_list |>
  mutate(
    n_analysis_rows = map_int(parent_id, \(id) sum(inventory_checked$parent_id == id)),
    diagnostic_readiness = case_when(
      parent_id == "P5" ~ "mixed: SNP and Trinity/Ross layers remain preliminary",
      parent_id == "P2" ~ "primary condition analyses covered; location-model fit checks remain partial",
      parent_id == "P7" ~ model_coverage,
      TRUE ~ "covered for current manuscript use"
    )
  )

write_csv(parent_checked,
          file.path(TBL_DIR, "38_analysis_diagnostic_parent_list.csv"))
write_csv(inventory_checked,
          file.path(TBL_DIR, "38_analysis_diagnostic_inventory.csv"))
write_csv(artifact_files,
          file.path(TBL_DIR, "38_diagnostic_artifact_status.csv"))

md_table <- function(dat) {
  knitr::kable(dat, format = "pipe")
}

doc_inventory <- inventory_checked |>
  select(parent_analysis, analysis_family, role, model_or_test,
         diagnostics_done, goodness_of_fit_done, diagnostic_outcome,
         remaining_gap)

doc_artifacts <- artifact_files |>
  transmute(artifact, scope, exists, status_counts)

doc <- c(
  "# Analysis and diagnostic inventory",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This is the parent list for the project. It answers a simple question:",
  "for each analysis family, have we checked the data, the model fit, or the",
  "right alternative quality check?",
  "",
  "A `not applicable` goodness-of-fit entry is not a gap. It means the row is a",
  "data audit, exact test, descriptive screen, or derived index rather than a",
  "fitted model with residuals. A `partial` entry means the analysis is useful",
  "as exploratory context but needs more diagnostics before becoming a formal",
  "manuscript claim.",
  "",
  "## Parent list",
  "",
  md_table(parent_checked |>
             select(parent_id, parent_analysis, plain_question,
                    n_analysis_rows, diagnostic_readiness)),
  "",
  "## Analysis-level table",
  "",
  md_table(doc_inventory),
  "",
  "## Diagnostic artifacts",
  "",
  md_table(doc_artifacts),
  "",
  "## Current bottom line",
  "",
  paste0("- Saved model coverage: ", model_coverage, "."),
  paste0("- Continuous mixed models: ", diag_counts$continuous, "."),
  paste0("- Morphology GLMMs: ", diag_counts$morph_glmm, "."),
  paste0("- Cox/time-to-event models: ", diag_counts$cox, "."),
  paste0("- PCA/source-patch LRT checks: ", diag_counts$pca_lrt, "."),
  paste0("- Formal trade-off models: ", diag_counts$tradeoff, "."),
  "- Main remaining diagnostic gap: the preliminary SNP-by-symbiont linear models need residual/influence checks before any SNP-phenotype result is treated as inferential. They are currently marked as exploratory.",
  "- Molly's recent LTH emails are incorporated here: A/C/D remain source-patch labels, Trinity/Ross ED50 values remain external context, and exact identity requires DNA-marker comparison rather than coordinate matching.",
  "- Environmental update: recovered Apex files fill Days 0-15 for all tanks. Exact clock-hour exposure still needs the logger time-zone convention confirmed; later YSI spot checks remain unavailable."
)

writeLines(doc, DOC_PATH)

if (any(inventory_checked$missing_scripts != "none")) {
  warning("Some inventory script paths are missing: ",
          paste(inventory_checked$missing_scripts[
            inventory_checked$missing_scripts != "none"
          ], collapse = " | "))
}
if (any(inventory_checked$missing_primary_outputs != "none")) {
  warning("Some inventory primary outputs are missing: ",
          paste(inventory_checked$missing_primary_outputs[
            inventory_checked$missing_primary_outputs != "none"
          ], collapse = " | "))
}

cat("\n=== Analysis diagnostic inventory ===\n")
cat(sprintf("Parent rows: %d | Analysis rows: %d | Diagnostic artifacts: %d\n",
            nrow(parent_checked), nrow(inventory_checked), nrow(artifact_files)))
cat(model_coverage, "\n")
cat("Wrote:\n")
cat(" - output/tables/38_analysis_diagnostic_parent_list.csv\n")
cat(" - output/tables/38_analysis_diagnostic_inventory.csv\n")
cat(" - output/tables/38_diagnostic_artifact_status.csv\n")
cat(" - docs/analysis_diagnostic_inventory.md\n")
