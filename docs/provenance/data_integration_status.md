# Data Integration Status

Last checked: 2026-09-01

Status categories:

- `full_analysis`: cleaned from raw data, written to processed/output artifacts,
  and included in the master pipeline.
- `supporting_analysis`: cleaned and included in the master pipeline, but mainly
  as design, calibration, environmental, or quality-control context rather than
  as a primary response model.
- `handoff_only`: used to build a downstream handoff table; no statistical
  analysis is run in this repo.
- `provenance_only`: imported and characterized, but not read by the current
  analysis pipeline.

| source group | status | current use |
|---|---|---|
| `metadata` | `supporting_analysis` | Cleaned by `code/01_load_clean_metadata.R`; writes `data/processed/coral_metadata.rds` and the metadata summary. This is the design spine for downstream joins and validation. |
| `pam` | `full_analysis` | Cleaned/analyzed by `code/02_pam_analysis.R`; writes `data/processed/pam_clean.rds`, PAM figures, and treatment contrasts. |
| `color_card` | `full_analysis` | Cleaned/analyzed by `code/03_color_card_analysis.R`; writes `data/processed/color_clean.rds`, color figures, end proportions, and contributes to the main models. |
| `physio_morphology` | `full_analysis` | Cleaned/analyzed by `code/04_physio_morphology.R`; writes `data/processed/physio_clean.rds`, morphology trajectory figures, GLMM summaries, and feeds survival/KM/model scripts. Also feeds the cross-dataset morphology plots, endpoint/timing tests, and diagnostic coverage checks in `code/11c_morphology_dataset_plots.R`. |
| `microscope_physio` | `full_analysis` | Cleaned/analyzed by `code/11_microscope_physio.R` as a separate photo-only cohort; writes cleaned data, design/event/trajectory tables, and microscope figures. Also feeds the cross-dataset morphology plots, endpoint/timing tests, and diagnostic coverage checks in `code/11c_morphology_dataset_plots.R`. |
| `wax_dipping` | `supporting_analysis` | Cleaned/analyzed by `code/07_wax_dipping.R`; writes `data/processed/wax_clean.rds`, the wax standard-curve figure, and surface areas used by growth/symbiont calculations. |
| `buoyant_weight` | `full_analysis` | Cleaned/analyzed by `code/05_buoyant_weight.R`; writes `data/processed/buoyant_weight_clean.rds`, growth figures, growth model tables, tank tests, and metric comparisons. |
| `symbiont_counts` | `full_analysis` | Cleaned/analyzed by `code/06_symbiont_chl.R`; writes `data/processed/symbiont_chl_clean.rds`, symbiont-density figure, summaries, and RNA-seq covariates. Chlorophyll-a is explicitly handled as not run. |
| `apex` | `supporting_analysis` | Parsed/analyzed by `code/08_apex_temperature.R`; writes hourly/daily temperature RDS files and thermal-context figures. |
| `ysi` | `supporting_analysis` | Cleaned/analyzed by `code/09_ysi_water_chem.R`; writes `data/processed/ysi_clean.rds` and water-chemistry figures. |
| `worm_presence` | `supporting_analysis` | Cleaned/summarized by `code/10_worms.R`; writes `data/processed/worm_clean.rds`, worm-presence figure, and summary table. Used as contamination/QC context. |
| `plate_layout` | `handoff_only` | Read by `code/31_rnaseq_covariate_table.R`; writes RNA-seq library lookup and phenotype-covariate handoff tables. No expression model is fit in this repo. |
| `shipping` | `provenance_only` | Imported/codebooked sample handling and storage metadata. The current pipeline does not read these sheets directly. |
| `daily_health_log` | `provenance_only` | Imported/codebooked for completeness; mostly a filled design/template log and not used by the current pipeline. |
| `printable_data_sheets` | `provenance_only` | Imported/codebooked collection templates for YSI, PAM, color-card, and daily-health logs; not analysis inputs. |
| `wax_dipping_standard_curve` | `provenance_only` | Standalone standard-curve copy retained for traceability. The current pipeline uses `data/raw/wax_dipping/Standard_curve.csv`. |
| `buoyant_weight_archive` | `provenance_only` | Duplicate/working buoyant-weight calculation workbook retained for traceability. The current pipeline uses `data/raw/buoyant_weight/` tab CSVs. |
| `chlorophyll` | `provenance_only` | Example calculation workbook only. Chlorophyll-a was planned but not populated/run for this project. |
| Drive method docs | `provenance_only` | Exported to `docs/provenance/drive_methods/`; used to document design, methods, photo-source paths, and interpretation. |
| NAS image sources | `provenance_only` | Local NAS paths and file counts are documented in `docs/provenance/nas_image_sources.md`; raw images are not imported or pixel-analyzed in this repo. |

Validation status after this classification:

`39 PASS, 3 HANDLED, 0 WARN, 0 FAIL`
