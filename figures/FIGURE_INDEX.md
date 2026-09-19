# Figure Index

Updated: 2026-09-09

This is the current figure catalog for the LTH heat x wound repo. The older
`FIGURE_INDEX.docx` is retained for provenance but may lag the current pipeline
and terminology. Reader-facing text should use **source patch** for A, C, and D;
some filenames still use legacy `genet` wording.

## Start Here

| Figure | Use |
|---|---|
| `08_apex_temperature.png` | Confirms the tanks held the 28 C and 31 C temperature treatments. |
| `14_morphology_KM.png` | Shows that heat mostly blocked new-corallite budding rather than wound sealing. |
| `36b_molly_stage_violin.png` | Violin-style healing/regeneration stage plot with points, non-reached fragments, and reached/scored counts. |
| `36_molly_stage_timing.png` | Supporting exact-day healing/regeneration stage plot with stacked points and X marks for stages not reached. |
| `19e_source_heat_penalty_summary.png` | Plain-language summary of source-patch heat penalties across response types. |
| `19f_source_physiology_heat_penalties.png` | Easier view of the whole-coral physiology and growth heat penalties. |
| `19_genet_dashboard.png` | Legacy-name detailed source-patch heat-penalty plot across measurements. |
| `32_prelim_snp_structure.png` | Preliminary RNA-seq SNP structure: C is coherent, A and D are mixed. |
| `26b_trinity_hauru_context.png` | Preliminary Trinity/Ross Hauru coordinate context; candidate matches only, not confirmed identities. |
| `26_thermal_context.png` | Shows that 31 C is below the acute CBASS ED50 range. |
| `35_tradeoff_summary.png` | Exploratory growth, heat-tolerance, and regeneration trade-off screen. |

## Current Pipeline Figures

| File stem | What it shows | Main script | Current use |
|---|---|---|---|
| `02_pam_fvfm_trajectory` | Photosynthesis-efficiency score over time | `code/02_pam_analysis.R` | Physiology result |
| `03_color_trajectory` | Color-card paling over time | `code/03_color_card_analysis.R` | Physiology result |
| `04_morphology_trajectories` | Eight wound-healing and regrowth traits over time | `code/04_physio_morphology.R` | Supporting morphology |
| `04b_morphology_trajectories_by_genet` | Same traits by legacy source-patch/genet label | `code/04_physio_morphology.R` | Supporting source-patch view |
| `05_buoyant_weight_growth` | Skeletal mass change | `code/05_buoyant_weight.R` | Growth result |
| `06_symbiont_density_by_day` | Algal symbiont density by day | `code/06_symbiont_chl.R` | Symbiont result |
| `07_wax_standard_curve` | Wax-dipping surface-area calibration | `code/07_wax_dipping.R` | Calibration |
| `08_apex_temperature` | Main Apex tank-temperature record | `code/08_apex_temperature.R` | Treatment validation |
| `08b_apex_temperature_full` | Full Apex temperature record | `code/08_apex_temperature.R` | Supporting treatment validation |
| `09_ysi_water_chem` | Daily YSI water chemistry checks | `code/09_ysi_water_chem.R` | Husbandry context |
| `10_worm_presence` | Flatworm presence by treatment/tank | `code/10_worms.R` | Contamination context |
| `11_microscope_trait_trajectories` | Microscope-photo healing and regrowth trajectories | `code/11_microscope_physio.R` | Supporting photo cohort |
| `11b_microscope_trait_trajectories_by_genet` | Microscope-photo trajectories by source patch | `code/11_microscope_physio.R` | Supporting photo cohort |
| `11c_morphology_dataset_all_trait_trajectories` | Main and microscope morphology datasets kept separate | `code/11c_morphology_dataset_plots.R` | Cross-dataset check |
| `11d_morphology_dataset_shared_trait_comparison` | Shared morphology traits across datasets | `code/11c_morphology_dataset_plots.R` | Cross-dataset check |
| `11e_morphology_dataset_event_summary` | Morphology event summary across datasets | `code/11c_morphology_dataset_plots.R` | Cross-dataset check |
| `13_genet_response_panel` | Legacy-name response panel by source-patch label | `code/13_genet_interaction.R` | Supporting source-patch view |
| `14_morphology_KM` | Kaplan-Meier time-to-onset curves | `code/14_morphology_kaplan.R` | Main morphology result |
| `14b_morphology_KM_by_genet` | Kaplan-Meier curves by source-patch label | `code/14_morphology_kaplan.R` | Supporting source-patch view |
| `15_physio_PCA_biplot` | End-of-study physiology PCA | `code/15_multivariate.R` | Multivariate condition |
| `15b_physio_PCA_by_genet` | Physiology PCA by source-patch label | `code/15_multivariate.R` | Supporting source-patch view |
| `16_manuscript_fig1` | Standalone four-panel overview | `code/16_main_figure.R` | Overview figure |
| `19_genet_dashboard` | Legacy-name detailed heat-penalty dashboard | `code/19_genet_dashboard.R` | Detailed source-patch summary |
| `19b_genet_resilience_ranking` | Composite heat-penalty ranking | `code/19_genet_dashboard.R` | Supporting source-patch summary |
| `19c_decomposed_resilience` | Heat penalties split by wounded/unwounded scope | `code/19_genet_dashboard.R` | Supporting source-patch summary |
| `19d_wound_healing_heat_penalties` | Healing and regrowth heat penalties | `code/19_genet_dashboard.R` | Wound-specific source-patch summary |
| `19e_source_heat_penalty_summary` | Plain-language heat-penalty summary | `code/19_genet_dashboard.R` | Team-summary Section 4 |
| `19f_source_physiology_heat_penalties` | Whole-coral heat penalties by response | `code/19_genet_dashboard.R` | Team-summary Section 4 |
| `26_thermal_context` | LTH temperatures versus acute CBASS ED50 | `code/sensitivity/26_thermal_context.R` | Thermal context |
| `26b_trinity_hauru_context` | Preliminary Trinity/Ross Hauru coordinate screen | `code/sensitivity/26b_trinity_hauru_context.R` | Candidate-match context |
| `32_prelim_snp_structure` | Preliminary SNP clusters and genetic PCs | `code/32_prelim_snp_phenotype_integration.R` | Preliminary genetics |
| `32_prelim_snp_design_balance` | Preliminary SNP clusters by treatment design | `code/32_prelim_snp_phenotype_integration.R` | Preliminary genetics |
| `32_prelim_snp_symbiont_heat_effects` | Direct SNP-by-symbiont heat-penalty check | `code/32_prelim_snp_phenotype_integration.R` | Preliminary genetics |
| `33_growth_heat_tradeoff` | Growth versus heat-penalty screen | `code/33_growth_heat_tradeoff.R` | Exploratory trade-off |
| `34_regeneration_tradeoff_screens` | Regeneration trade-off screens | `code/34_regeneration_tradeoff_screens.R` | Exploratory trade-off |
| `35_tradeoff_summary` | Four-panel trade-off summary | `code/35_tradeoff_summary.R` | Team-summary trade-off check |
| `36_molly_stage_timing` | Exact-day healing/regeneration stage timing, with unreached late stages shown explicitly | `code/36_molly_followup_checks.R` | Supporting output for Molly follow-up |
| `36b_molly_stage_violin` | Violin-style healing/regeneration stage timing, with individual points, non-reached stages, and reached/scored labels | `code/36_molly_followup_checks.R` | Team-summary Section 3 and Molly follow-up |
| `41_stage_progression` | **Recommended stage-timing display.** Per-fragment first-observation day, an explicit never-reached zone on the same axis, and the proportion reached with binomial intervals. Suppresses a median wherever <80% of a group reached the milestone. | `code/41_stage_progression_figure.R` | Milestone timing — prefer over 36/39a/39c/40 |

## Notes

- **Stage-timing figures:** use `41_stage_progression`. The earlier options are kept
  as the design comparison they were commissioned as, but `40_stage_mean_review`
  plots a mean first-observation day among fragments that *reached* the milestone,
  which renders the headline result as a 0.9-day gap (12.6 vs 13.5 d) when the
  finding is 12/12 versus 4/12; and `39a_stage_timing_violin_option` draws kernel
  densities from at most 12 reachers. Do not put either in the manuscript.
- Most figure stems have both `.png` and `.pdf` outputs.
- `figures/12_diagnostics/` contains diagnostic plots for model checks and is
  not part of the reader-facing figure set.
- The Trinity/Ross Hauru figure is preliminary context only. It should not be used to assign
  ED50 values to LTH source patches until DNA-marker matching is complete.
