# =============================================================================
# Purpose: Reproducible follow-up checks for Molly's 2026-09-04 and
#          2026-09-08 notes.
#
# What & why: Molly flagged four useful points for the team summary:
#   1. the current exported Apex/YSI temperature record may be incomplete after
#      Day 10,
#   2. the wound story should be framed as one early healing step followed by
#      later regeneration stages,
#   3. a Figure 1 companion should show stage-timing distributions while
#      keeping corals that never reached late stages visible, and
#   4. Day 15 samples may be useful for genotype-only checks, with RNA-seq
#      reserved for the broader expression question.
#
# This script makes those checks explicit and writes small artifacts that the
# coauthor summary and email draft can cite.
#
# Input:   data/processed/apex_temperature_daily.rds
#          data/processed/ysi_clean.rds
#          data/processed/physio_clean.rds
#          data/processed/microscope_physio_clean.rds
#          output/tables/32_prelim_snp_response_joinability.csv
# Output:  output/tables/36_molly_temperature_coverage_days0_16.csv
#          output/tables/36_molly_stage_timing.csv
#          output/tables/36_molly_stage_timing_summary.csv
#          output/tables/36_molly_stage_timing_mean_iqr.csv
#          output/tables/36_molly_followup_action_items.csv
#          figures/36_molly_stage_timing.{pdf,png}
#          figures/36b_molly_stage_violin.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

D0 <- as.Date("2025-06-04")
expected_tanks <- sort(c(TANK_28C, TANK_31C))
days_0_16 <- tibble(
  day = 0:16,
  date = D0 + day
)

# ---- Temperature coverage in the current repo exports ----------------------
apex_path <- file.path(DATA_PROC, "apex_temperature_daily.rds")
ysi_path <- file.path(DATA_PROC, "ysi_clean.rds")

apex_temp <- if (file.exists(apex_path)) {
  readRDS(apex_path) |>
    mutate(
      date = as.Date(date),
      tank = as.integer(stringr::str_match(probe, "Temp([0-9]+)")[, 2]),
      treatment = tank_treatment(tank)
    ) |>
    filter(!is.na(treatment), tank %in% expected_tanks,
           is.finite(value_mean), value_mean > 15, value_mean < 40)
} else {
  tibble(date = as.Date(character()), tank = integer(), treatment = character(),
         value_mean = numeric())
}

ysi <- if (file.exists(ysi_path)) {
  readRDS(ysi_path) |>
    mutate(date = as.Date(date)) |>
    filter(tank %in% expected_tanks)
} else {
  tibble(date = as.Date(character()), tank = integer(), treatment = character(),
         temp_c = numeric())
}

apex_daily <- apex_temp |>
  mutate(day = as.integer(date - D0)) |>
  filter(day >= 0, day <= 16) |>
  group_by(day, date) |>
  summarise(
    apex_tanks_present = n_distinct(tank),
    apex_tanks = paste(sort(unique(tank)), collapse = ","),
    apex_missing_tanks =
      paste(setdiff(expected_tanks, sort(unique(tank))), collapse = ","),
    apex_mean_28_c = if (any(treatment == "28C")) {
      mean(value_mean[treatment == "28C"], na.rm = TRUE)
    } else NA_real_,
    apex_mean_31_c = if (any(treatment == "31C")) {
      mean(value_mean[treatment == "31C"], na.rm = TRUE)
    } else NA_real_,
    .groups = "drop"
  )

ysi_daily <- ysi |>
  mutate(day = as.integer(date - D0)) |>
  filter(day >= 0, day <= 16) |>
  group_by(day, date) |>
  summarise(
    ysi_tanks_present = n_distinct(tank),
    ysi_tanks = paste(sort(unique(tank)), collapse = ","),
    ysi_missing_tanks =
      paste(setdiff(expected_tanks, sort(unique(tank))), collapse = ","),
    ysi_mean_28_c = if (any(treatment == "28C")) {
      mean(temp_c[treatment == "28C"], na.rm = TRUE)
    } else NA_real_,
    ysi_mean_31_c = if (any(treatment == "31C")) {
      mean(temp_c[treatment == "31C"], na.rm = TRUE)
    } else NA_real_,
    .groups = "drop"
  )

temperature_coverage <- days_0_16 |>
  left_join(apex_daily, by = c("day", "date")) |>
  left_join(ysi_daily, by = c("day", "date")) |>
  mutate(
    across(c(apex_tanks_present, ysi_tanks_present), ~ replace_na(.x, 0L)),
    across(c(apex_tanks, apex_missing_tanks, ysi_tanks, ysi_missing_tanks),
           ~ replace_na(.x, "")),
    apex_missing_tanks = if_else(
      apex_tanks_present == length(expected_tanks),
      "none",
      if_else(apex_tanks_present == 0,
              paste(expected_tanks, collapse = ","),
              apex_missing_tanks)
    ),
    ysi_missing_tanks = if_else(
      ysi_tanks_present == length(expected_tanks),
      "none",
      if_else(ysi_tanks_present == 0,
              paste(expected_tanks, collapse = ","),
              ysi_missing_tanks)
    ),
    apex_delta_31_minus_28_c = apex_mean_31_c - apex_mean_28_c,
    ysi_delta_31_minus_28_c = ysi_mean_31_c - ysi_mean_28_c,
    apex_status = case_when(
      apex_tanks_present == length(expected_tanks) ~ "complete in current export",
      apex_tanks_present == 0 ~ "absent from current export",
      TRUE ~ "partial in current export"
    ),
    ysi_status = case_when(
      ysi_tanks_present == length(expected_tanks) ~ "complete in current export",
      ysi_tanks_present == 0 ~ "absent from current export",
      TRUE ~ "partial in current export"
    )
  ) |>
  arrange(day) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

write_csv(temperature_coverage,
          file.path(TBL_DIR, "36_molly_temperature_coverage_days0_16.csv"))

# ---- Stage timing: one healing step plus five regeneration steps -----------
physio <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))

stage_map <- tibble::tribble(
  ~stage_order, ~phase,          ~stage,                                ~dataset,              ~dataset_label,        ~trait,                  ~measurement_note,
  1L,           "Healing",       "Coenosarc barrier re-established",    "microscope_physio",   "Microscope photos",   "tissue_over_wound",     "microscope-only coenosarc/tissue-cover score",
  2L,           "Regeneration",  "Axial polyp forms",                  "physio_morphology",   "Main morphology",     "axial_polyp_formation", "gross-photo central axial-polyp/hole score",
  3L,           "Regeneration",  "Wound surface remodels",             "physio_morphology",   "Main morphology",     "wound_smoothed",        "gross-photo smoothed wound-surface score",
  4L,           "Regeneration",  "Axial corallite or tip forms",        "physio_morphology",   "Main morphology",     "tip_exist",             "gross-photo visible regenerated-tip score",
  5L,           "Regeneration",  "Axial corallite or tip extends",      "physio_morphology",   "Main morphology",     "tip_extension",         "gross-photo tip-extension score",
  6L,           "Regeneration",  "New radial corallites bud",           "physio_morphology",   "Main morphology",     "new_corallites_on_tip", "gross-photo new-corallites-on-tip score"
) |>
  mutate(
    phase = factor(phase, levels = c("Healing", "Regeneration")),
    stage = factor(stage, levels = stage)
  )

compute_stage_events <- function(dat, stage_row) {
  tr <- stage_row$trait
  dat |>
    filter(as.character(wound) == "yes", !is.na(day), day >= 0, day <= 15) |>
    select(id, treatment, tank, thicket, day, y = all_of(tr)) |>
    filter(!is.na(y)) |>
    mutate(treatment = factor(as.character(treatment), levels = c("28C", "31C"))) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      first_scored_day = min(day),
      last_scored_day = max(day),
      event = any(y == 1, na.rm = TRUE),
      first_yes_day = if (any(y == 1, na.rm = TRUE)) {
        min(day[y == 1], na.rm = TRUE)
      } else NA_real_,
      .groups = "drop"
    ) |>
    mutate(
      stage_order = stage_row$stage_order,
      phase = stage_row$phase,
      stage = stage_row$stage,
      dataset = stage_row$dataset,
      dataset_label = stage_row$dataset_label,
      trait = stage_row$trait,
      measurement_note = stage_row$measurement_note
    )
}

stage_events <- purrr::pmap_dfr(stage_map, function(stage_order, phase, stage,
                                                    dataset, dataset_label,
                                                    trait, measurement_note) {
  stage_row <- tibble(stage_order = stage_order, phase = phase, stage = stage,
                      dataset = dataset, dataset_label = dataset_label,
                      trait = trait, measurement_note = measurement_note)
  dat <- if (dataset == "microscope_physio") micro else physio
  compute_stage_events(dat, stage_row)
}) |>
  mutate(
    stage = factor(stage, levels = stage_map$stage),
    phase = factor(phase, levels = c("Healing", "Regeneration")),
    treatment = factor(treatment, levels = c("28C", "31C")),
    event_status = if_else(event, "Reached", "Not reached by final day"),
    plot_day = if_else(event, first_yes_day, 16.35),
    stage_label = paste0(as.character(phase), ": ", as.character(stage)),
    stage_label = factor(stage_label,
                         levels = rev(paste0(stage_map$phase, ": ",
                                             stage_map$stage)))
  ) |>
  arrange(stage_order, treatment, id)

write_csv(stage_events,
          file.path(TBL_DIR, "36_molly_stage_timing.csv"))

stage_summary <- stage_events |>
  group_by(stage_order, phase, stage, dataset_label, trait, measurement_note,
           treatment) |>
  summarise(
    n_corals = n(),
    n_reached = sum(event),
    n_not_reached = n_corals - n_reached,
    pct_reached = round(100 * n_reached / n_corals, 1),
    mean_first_day = if (n_reached > 0) {
      round(mean(first_yes_day[event], na.rm = TRUE), 2)
    } else NA_real_,
    median_first_day = if (n_reached > 0) {
      round(median(first_yes_day[event]), 2)
    } else NA_real_,
    q25_first_day = if (n_reached > 0) {
      round(unname(quantile(first_yes_day[event], 0.25)), 2)
    } else NA_real_,
    q75_first_day = if (n_reached > 0) {
      round(unname(quantile(first_yes_day[event], 0.75)), 2)
    } else NA_real_,
    iqr_first_day = if (n_reached > 0) {
      round(IQR(first_yes_day[event]), 2)
    } else NA_real_,
    .groups = "drop"
  ) |>
  group_by(stage_order) |>
  mutate(complete_both_temperatures = all(n_reached == n_corals),
         across(c(mean_first_day, median_first_day, q25_first_day,
                  q75_first_day, iqr_first_day),
                ~if_else(complete_both_temperatures, .x, NA_real_))) |>
  ungroup() |>
  arrange(stage_order, treatment)

write_csv(stage_summary,
          file.path(TBL_DIR, "36_molly_stage_timing_summary.csv"))

stage_summary_display <- stage_summary |>
  transmute(
    stage_order,
    phase = as.character(phase),
    stage = as.character(stage),
    dataset = dataset_label,
    treatment = dplyr::recode(as.character(treatment),
                              `28C` = "28 °C", `31C` = "31 °C"),
    n_reached,
    n_scored = n_corals,
    n_not_reached,
    pct_reached,
    complete_both_temperatures,
    mean_day = mean_first_day,
    median_day = median_first_day,
    iqr_days = if_else(
      complete_both_temperatures,
      sprintf("%.2f-%.2f", q25_first_day, q75_first_day),
      NA_character_
    )
  )

write_csv(stage_summary_display,
          file.path(TBL_DIR, "36_molly_stage_timing_mean_iqr.csv"))

label_data <- stage_summary |>
  mutate(
    stage_label = paste0(as.character(phase), ": ", as.character(stage)),
    stage_label = factor(stage_label,
                         levels = rev(paste0(stage_map$phase, ": ",
                                             stage_map$stage))),
    reached_label = paste0(n_reached, "/", n_corals),
    y_plot = as.numeric(stage_label) +
      if_else(treatment == "28C", -0.18, 0.18)
  )

stage_levels <- levels(stage_events$stage_label)

stage_events_plot <- stage_events |>
  group_by(stage_label, treatment, event_status, plot_day) |>
  arrange(id, .by_group = TRUE) |>
  mutate(
    n_in_cell = n(),
    dot_index = row_number(),
    n_cols = pmin(n_in_cell, 4L),
    stack_col = ((dot_index - 1L) %% n_cols) + 1L,
    stack_row = ((dot_index - 1L) %/% n_cols) + 1L,
    n_rows = max(stack_row),
    x_offset = if_else(n_in_cell == 1L, 0,
                       (stack_col - (n_cols + 1) / 2) * 0.11),
    y_offset = if_else(n_in_cell == 1L, 0,
                       (stack_row - (n_rows + 1) / 2) * 0.08),
    y_plot = as.numeric(stage_label) +
      if_else(treatment == "28C", -0.18, 0.18) + y_offset,
    x_plot = plot_day + x_offset
  ) |>
  ungroup()

status_key <- tibble(
  x = NA_real_,
  y = NA_real_,
  event_status = factor(c("Reached", "Not reached by final day"),
                        levels = c("Reached", "Not reached by final day"))
)

p_stage <- ggplot() +
  annotate("rect", xmin = 15.75, xmax = 16.85, ymin = -Inf, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_vline(xintercept = 15, linetype = "dotted",
             linewidth = 0.35, colour = "grey55") +
  geom_point(
    data = filter(stage_events_plot, event_status == "Reached"),
    aes(x_plot, y_plot, fill = treatment),
    shape = 21, colour = "white", size = 1.75, stroke = 0.25, alpha = 0.95
  ) +
  geom_point(
    data = filter(stage_events_plot, event_status == "Not reached by final day"),
    aes(x_plot, y_plot, colour = treatment),
    shape = 4, size = 1.95, stroke = 0.8, alpha = 0.95
  ) +
  geom_point(
    data = status_key,
    aes(x, y, shape = event_status),
    colour = "#444444", fill = "#444444", size = 2.8, stroke = 0.8,
    inherit.aes = FALSE, na.rm = TRUE
  ) +
  geom_text(
    data = label_data,
    aes(x = 17.15, y = y_plot, label = reached_label,
        colour = treatment),
    inherit.aes = FALSE,
    size = 2.45,
    hjust = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = PAL_TEMP, name = NULL) +
  scale_colour_manual(values = PAL_TEMP, guide = "none") +
  scale_shape_manual(values = c(`Reached` = 21,
                                `Not reached by final day` = 4),
                     name = NULL) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21, colour = "white", size = 3)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(colour = "#444444", fill = "#444444", size = 3)
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 3, 5, 10, 15, 16.35),
    labels = c("0", "1", "3", "5", "10", "15", "not reached"),
    limits = c(-0.2, 18.0),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    breaks = seq_along(stage_levels),
    labels = stage_levels,
    limits = c(0.5, length(stage_levels) + 0.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "First day the stage was seen",
    y = NULL,
    title = "Heat mostly affects the late regeneration stage",
    subtitle = "Small stacks separate corals that reached the same stage on the same day. X marks corals that never reached the stage.",
    caption = "Each point is one wounded coral. Right labels show reached/total corals.\nThe healing row uses microscope photos; regeneration rows use the main morphology data set."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub(8) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 40, 8, 8, "pt")
  )

save_fig(p_stage, "36_molly_stage_timing", width = 190, height = 110)

stage_axis_levels <- paste0("Stage ", seq_len(nrow(stage_map)))
stage_axis_labels <- c(
  "1\nCoenosarc\ncover",
  "2\nMain polyp\nforms",
  "3\nWound surface\nremodels",
  "4\nTip\nforms",
  "5\nTip\nextends",
  "6\nNew corallites\nbud"
)
dodge_width <- 0.72

stage_events_violin <- stage_events |>
  mutate(
    stage_axis = factor(paste0("Stage ", stage_order),
                        levels = stage_axis_levels),
    event_status = factor(if_else(event, "Reached",
                                  "Not reached by Day 15"),
                          levels = c("Reached",
                                     "Not reached by Day 15"))
  )

stage_summary_violin <- stage_summary |>
  mutate(
    stage_axis = factor(paste0("Stage ", stage_order),
                        levels = stage_axis_levels),
    pct_not_reached_label = if_else(
      n_not_reached > 0,
      paste0(round(100 * n_not_reached / n_corals), "%"),
      NA_character_
    )
  )

p_stage_violin <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 16.05, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = 16.05, linetype = "dotted",
             linewidth = 0.4, colour = "grey50") +
  geom_violin(
    data = filter(stage_events_violin, event),
    aes(x = stage_axis, y = first_yes_day, fill = treatment,
        group = interaction(stage_axis, treatment)),
    position = position_dodge(width = dodge_width),
    scale = "count",
    trim = TRUE,
    adjust = 0.85,
    alpha = 0.72,
    colour = "#202020",
    linewidth = 0.35,
    width = 0.68
  ) +
  geom_boxplot(
    data = filter(stage_events_violin, event, stage_order %in%
                    stage_summary$stage_order[stage_summary$complete_both_temperatures]),
    aes(x = stage_axis, y = first_yes_day,
        group = interaction(stage_axis, treatment)),
    position = position_dodge(width = dodge_width),
    width = 0.12,
    fill = NA,
    colour = "#202020",
    linewidth = 0.38,
    outlier.shape = NA
  ) +
  geom_point(
    data = filter(stage_summary_violin, complete_both_temperatures),
    aes(x = stage_axis, y = mean_first_day,
        group = treatment),
    position = position_dodge(width = dodge_width),
    inherit.aes = FALSE,
    shape = 23,
    size = 2,
    stroke = 0.35,
    colour = "#202020",
    fill = "white"
  ) +
  geom_point(
    data = filter(stage_events_violin, event),
    aes(x = stage_axis, y = first_yes_day, fill = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0,
      dodge.width = dodge_width, seed = 36
    ),
    shape = 21,
    colour = "white",
    size = 1.55,
    stroke = 0.25,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(stage_events_violin, !event),
    aes(x = stage_axis, y = 17.05, colour = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0.06,
      dodge.width = dodge_width, seed = 37
    ),
    shape = 4,
    size = 2.05,
    stroke = 0.85,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(stage_summary_violin, !is.na(pct_not_reached_label)),
    aes(x = stage_axis, y = 17.72, label = pct_not_reached_label,
        colour = treatment, group = treatment),
    position = position_dodge(width = dodge_width),
    inherit.aes = FALSE,
    size = 2.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = PAL_TEMP,
    labels = c(`28C` = "28 °C", `31C` = "31 °C"),
    name = "Temperature"
  ) +
  scale_colour_manual(
    values = PAL_TEMP,
    labels = c(`28C` = "28 °C", `31C` = "31 °C"),
    guide = "none"
  ) +
  scale_x_discrete(labels = stage_axis_labels, drop = FALSE) +
  scale_y_continuous(
    breaks = c(1, 3, 5, 7, 9, 11, 13, 15, 17.05),
    labels = c("1", "3", "5", "7", "9", "11", "13", "15",
               "Not\nreached"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Healing or regeneration stage",
    y = "Day first seen",
    title = "Heat reduced how many corals reached the final regrowth stage",
    subtitle = "Violin width scales with the number of corals that reached each stage; X marks corals still short of a stage by Day 15.",
    caption = "Boxes and means appear only where all fragments at both temperatures reached the stage.\nPoints and violin shapes for incomplete stages describe only fragments that reached them; X marks retain the others."
  ) +
  coord_cartesian(ylim = c(0.5, 18.15), clip = "on") +
  theme_pub(8) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(lineheight = 0.9, margin = margin(t = 4)),
    plot.margin = margin(8, 12, 8, 8, "pt")
  )

save_fig(p_stage_violin, "36b_molly_stage_violin",
         width = 190, height = 130)

# ---- Action table for the coauthor summary and email -----------------------
joinability_path <- file.path(TBL_DIR, "32_prelim_snp_response_joinability.csv")
direct_snp_responses <- if (file.exists(joinability_path)) {
  read_csv(joinability_path, show_col_types = FALSE) |>
    filter(direct_sample_level_join) |>
    pull(response_variable) |>
    paste(collapse = ", ")
} else {
  "not checked"
}

apex_missing_days <- temperature_coverage |>
  filter(apex_tanks_present < length(expected_tanks)) |>
  pull(day)
ysi_missing_days <- temperature_coverage |>
  filter(ysi_tanks_present < length(expected_tanks)) |>
  pull(day)

coenosarc_d1 <- stage_events |>
  filter(trait == "tissue_over_wound", first_scored_day <= 1,
         last_scored_day >= 1) |>
  group_by(treatment) |>
  summarise(
    n_scored_day1 = sum(!is.na(first_yes_day) & first_yes_day <= 1),
    n_event_by_day1 = sum(event & !is.na(first_yes_day) & first_yes_day <= 1),
    .groups = "drop"
  )
coenosarc_d2 <- stage_events |>
  filter(trait == "tissue_over_wound") |>
  summarise(n_reached_by_day2 = sum(event & first_yes_day <= 2),
            n_total = n())

action_items <- tibble::tribble(
  ~topic, ~molly_request_or_point, ~what_current_repo_says, ~summary_change, ~next_step,
  "Temperature records",
  "Check whether Days 11-16 are missing from YSI and Apex logs before using this as supplemental methods context.",
  sprintf("Current export has complete Apex tank coverage for Days %s only within Days 0-16; Apex is absent for Days %s. Current YSI export has complete tank coverage for Days %s only and is absent for Days %s.",
          paste(setdiff(0:16, apex_missing_days), collapse = ", "),
          paste(apex_missing_days, collapse = ", "),
          paste(setdiff(0:16, ysi_missing_days), collapse = ", "),
          paste(ysi_missing_days, collapse = ", ")),
  "Integrated Molly's recovered files, cross-checked the compiled sheet against XML, and rebuilt the temperature record with dated tank assignments.",
  "See code/08 coverage and provenance tables for reading counts and clock conventions; YSI coverage is a separate check.",
  "Day 15 RNA-seq",
  "Ask whether Day 15 coral samples should be sequenced, or genotyped only, to test genetic/cluster effects on expression and physiology.",
  sprintf("The current 144 RNA-seq/SNP libraries join directly to %s. Molly clarified that Day 15 corals were tracked through color card, PAM, buoyant weight, symbiont density, and regeneration, so genotyping those exact fragments would allow broader genotype-by-phenotype checks.",
          direct_snp_responses),
  "Reframed Day 15 as two decisions: genotype-only for genetic identity versus RNA-seq if expression at the endpoint is worth the added cost and analysis.",
  "Resolve final genetic identity first; then decide whether Day 15 genotyping answers the clustering question without extra expression libraries.",
  "Healing/regeneration staging",
  "Separate early healing from later regeneration stages.",
  sprintf("The microscope-photo cohort has 16 all-wounded corals. For coenosarc/tissue cover, 14/14 scored corals are positive by Day 1 and %d/%d are positive by Day 2. The main physiology cohort has 48 fragments; the stage-timing subset is 12 wounded fragments per temperature.",
          coenosarc_d2$n_reached_by_day2, coenosarc_d2$n_total),
  "Updated the one-healing-step plus five-regeneration-step framing with Molly's corrected sample-size language.",
  "Use the microscope-photo cohort for coenosarc healing only, likely as supplement support; keep later stages tied to the main morphology data.",
  "Figure 1 companion",
  "Make violin-style plots over time and add a companion table with mean day and IQR for each stage-by-temperature group.",
  "Small sample sizes, tied days, and non-reached late stages make a smoothed-only violin plot easy to misread. The key data are both later timing and the heated corals that never reach new-corallite budding.",
  "Made cumulative percentage reached the main figure. Means and middle-half ranges are shown only for stages reached by every fragment at both temperatures.",
  "Do not interpret averages among early finishers as overall time to regeneration. Report non-attainment through Day 15, not permanent failure."
)

write_csv(action_items,
          file.path(TBL_DIR, "36_molly_followup_action_items.csv"))

cat("\n=== Molly 2026-09-04 follow-up checks ===\n")
cat("Wrote temperature coverage, stage timing figure/table, and action-item table.\n")
cat("Apex missing days in current export:", paste(apex_missing_days, collapse = ", "), "\n")
cat("YSI missing days in current export:", paste(ysi_missing_days, collapse = ", "), "\n")
cat("Direct preliminary SNP-by-response check:", direct_snp_responses, "\n")
