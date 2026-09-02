# =============================================================================
# Purpose: Symbiont density (cells/cm^2); chlorophyll-a was planned but not run.
#          Endpoint physiological response variables sampled destructively at
#          each biopsy timepoint.
#          Joins symbiont counts with metadata. The chl-a column is retained only
#          as provenance for the planned-but-not-collected assay.
#
# What & why: corals get most of their energy from symbiotic algae
#   (zooxanthellae) living in their tissue; bleaching is the loss of these
#   symbionts under heat stress. We measured symbiont density as cells per cm^2
#   of coral surface. To get it, tissue was removed from a known surface area into
#   a slurry of known volume, and cells in that slurry were counted on a
#   hemocytometer (a gridded counting chamber). This sampling is DESTRUCTIVE, so
#   there is exactly one observation per coral — no repeated measures within a
#   coral. This script reconstructs cells/cm^2 from the raw quadrant counts (and
#   falls back to the pre-computed sheet average when raw counts are missing),
#   then summarizes and plots density by timepoint and temperature. Densities are
#   large and right-skewed, so downstream models log-transform them; here we keep
#   the natural cells/cm^2 scale for the cleaned data and figure. Chlorophyll-a
#   was planned but never assayed, so that column is carried only for provenance
#   and the figure shows a placeholder panel for it.
# Input:   data/raw/symbiont_counts/Raw_counts.csv  (4 hemocytometer reps per coral)
#          data/raw/symbiont_counts/metadata_ordered_merge.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/symbiont_chl_clean.rds
#          figures/06_symbiont_density_by_day.{pdf,png}
#          output/tables/06_symbiont_chl_summary.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_RAW, DATA_PROC, TBL_DIR),
# theme_pub(), and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# Three inputs: raw hemocytometer counts, the per-coral count metadata, and the
# master coral metadata (for the chl-a provenance column). suppressWarnings +
# guess_max keeps readr from choking on the mixed/blank cells in these exports;
# clean_names() standardises every header to snake_case.
raw_counts <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts", "Raw_counts.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names()
zoox_meta  <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts",
                     "metadata_ordered_merge.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names()
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))

# ---- Reps per coral -------------------------------------------------------
# Raw_counts has 4 hemocytometer quadrant counts (Q1, Q2, Q3, Q4) per
# replicate count per coral, in a "long" tidy layout.
# Columns: day, (blank), coral_id, count, q1, q2, q3, q4, average, coral_id_2
# Some rows have Q values, some don't (partial cells). Aggregate per coral_id:
reps <- raw_counts |>
  rename(coral_id = coral_id_3) |>
  filter(!is.na(coral_id)) |>
  mutate(coral_id = as.integer(coral_id)) |>
  mutate(across(c(q1, q2, q3, q4), as.numeric)) |>
  rowwise() |>     # rowwise so mean() averages the 4 quadrants WITHIN each row,
  mutate(quad_mean = mean(c(q1, q2, q3, q4), na.rm = TRUE)) |>  # not down a column
  ungroup() |>
  group_by(coral_id) |>
  summarise(
    mean_quad_count = mean(quad_mean, na.rm = TRUE),  # average the reps -> 1 value/coral
    n_reps          = n(),                            # how many count reps backed it
    .groups = "drop"
  ) |>
  filter(is.finite(mean_quad_count))

# ---- Build per-coral table & compute cells/cm^2 ----------------------------
# Join with per-coral metadata (treatment, day, SA, etc.)
# metadata_ordered_merge has columns: species, thicket, id_3, sample, wound,
# tank, treatment, biopsy_day, biopsy_date, average, coral_id, id_12,
# calculated_sa_standard_curve, slurry_volume_m_l
zoox <- zoox_meta |>
  rename(thicket = matches("^thicket"),
         id      = id_3) |>
  mutate(
    id              = as.integer(id),
    treatment       = factor(as.integer(treatment), levels = c(28, 31),
                              labels = c("28C", "31C")),     # baseline = 28C ambient
    biopsy_day      = as.integer(biopsy_day),
    thicket         = str_to_lower(str_squish(thicket)),     # normalize genet label
    wound           = factor(wound, levels = c("no", "yes")),  # baseline = unwounded
    sa_cm2          = as.numeric(calculated_sa_standard_curve),  # wax-curve SA (cm^2)
    slurry_ml       = as.numeric(slurry_volume_m_l),         # tissue slurry volume (mL)
    zoox_avg_sheet  = as.numeric(average),                   # spreadsheet's mean count
    # Hemocytometer formula: density = (cells/quadrant) * 10^4 (cells/mL)
    # (the 10^4 is the chamber's standard volume factor: a quadrant holds 1e-4 mL)
    # Total in slurry = density * slurry_mL ; then per cm^2 = / SA
    cells_per_cm2_sheet = (zoox_avg_sheet * 1e4 * slurry_ml) / sa_cm2
  ) |>
  filter(!is.na(id)) |>
  left_join(reps, by = c("id" = "coral_id")) |>   # attach the re-derived raw counts
  mutate(
    # Prefer the count recomputed from raw quadrants; fall back to the sheet's
    # average only when raw counts are missing for that coral.
    zoox_avg_hemo = coalesce(mean_quad_count, zoox_avg_sheet),
    count_source = if_else(!is.na(mean_quad_count), "raw_counts", "metadata_average"),
    # Hemocytometer formula: density = (cells/quadrant) * 10^4 (cells/mL)
    # Total in slurry = density * slurry_mL ; then per cm^2 = / SA
    cells_per_cm2 = (zoox_avg_hemo * 1e4 * slurry_ml) / sa_cm2
  ) |>
  mutate(tank = as.integer(tank)) |>
  select(id, treatment, wound, biopsy_day, thicket, tank, sa_cm2,
         cells_per_cm2, count_source, n_reps)

# Pull chlorophyll from master metadata (carried for provenance; assay not run,
# so this is empty for this project).
chl <- meta |>
  filter(!is.na(chlorophyll_ug_cm2)) |>
  select(id, chlorophyll_ug_cm2)

# Final analysis table: one row per coral, symbiont density + chl-a placeholder.
phys <- zoox |>
  left_join(chl, by = "id")

saveRDS(phys, file.path(DATA_PROC, "symbiont_chl_clean.rds"))

# ---- Summary --------------------------------------------------------------
# Cell means and standard errors per temperature x wound x timepoint, for the
# manuscript summary table. SE = sd / sqrt(n).
summary_tbl <- phys |>
  group_by(treatment, wound, biopsy_day) |>
  summarise(
    n             = n(),
    mean_cells    = mean(cells_per_cm2, na.rm = TRUE),
    se_cells      = sd(cells_per_cm2, na.rm = TRUE) / sqrt(n()),
    mean_chl      = mean(chlorophyll_ug_cm2, na.rm = TRUE),   # NaN when chl absent
    se_chl        = sd(chlorophyll_ug_cm2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(TBL_DIR, "06_symbiont_chl_summary.csv"))

# ---- Figure ----------------------------------------------------------------
# Express counts in millions for readable axes
p_zoox <- ggplot(phys, aes(factor(biopsy_day), cells_per_cm2 / 1e6,
                            fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.8) +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  labs(x = "Biopsy day",
       y = expression(Symbionts~(10^6~cells~cm^{-2})),
       title = "Symbiont density loss under heating",
       subtitle = expression(italic(A.~pulchra)~"biopsies, n = 192 across 4 timepoints")) +
  theme_pub(10)

# Chlorophyll-a was never run for this project (the chl-a column is retained
# only as a provenance placeholder), so we emit the single-panel symbiont-
# density figure. There is no chlorophyll panel.
save_fig(p_zoox, "06_symbiont_density_by_day", width = 140, height = 100)

cat("Wrote symbiont_chl_clean.rds, 06_symbiont_chl_summary.csv,",
    "06_symbiont_density_by_day.{pdf,png}\n")
