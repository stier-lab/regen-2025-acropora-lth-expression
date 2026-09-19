# =============================================================================
# Purpose: Compare three visual designs for healing/regeneration stage timing.
#
# What & why: Molly's violin-plot suggestion is one good option, but the same
# data can be easier to read in other formats depending on the audience. This
# script makes three reproducible alternatives from the stage-timing table:
#   A. count-scaled violins, close to Molly's mock-up;
#   B. cumulative percent reached through time;
#   C. a treatment-by-stage-day heatmap.
#
# Input:   output/tables/36_molly_stage_timing.csv
# Output:  figures/39a_stage_timing_violin_option.{pdf,png}
#          figures/39b_stage_timing_cumulative_option.{pdf,png}
#          figures/39c_stage_timing_heatmap_option.{pdf,png}
#          figures/39_stage_timing_visual_options.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

stage_events <- readr::read_csv(
  file.path(TBL_DIR, "36_molly_stage_timing.csv"),
  show_col_types = FALSE
) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    treatment_label = factor(
      dplyr::recode(as.character(treatment),
                    `28C` = "28 °C", `31C` = "31 °C"),
      levels = c("28 °C", "31 °C")
    ),
    stage_order = as.integer(stage_order),
    stage_short = factor(
      dplyr::recode(
        as.character(stage),
        `Coenosarc barrier re-established` = "Living tissue\ncovers wound",
        `Axial polyp forms` = "Main polyp\nforms",
        `Wound surface remodels` = "Wound surface\nremodels",
        `Axial corallite or tip forms` = "Tip\nforms",
        `Axial corallite or tip extends` = "Tip\nextends",
        `New radial corallites bud` = "New skeletal\ncups bud"
      ),
      levels = c("Living tissue\ncovers wound", "Main polyp\nforms",
                 "Wound surface\nremodels", "Tip\nforms",
                 "Tip\nextends", "New skeletal\ncups bud")
    ),
    event_status = factor(if_else(event, "Reached", "Not reached by Day 15"),
                          levels = c("Reached", "Not reached by Day 15"))
  )

stage_summary <- stage_events |>
  group_by(stage_order, stage_short, treatment, treatment_label) |>
  summarise(
    n_scored = n(),
    n_reached = sum(event),
    n_not_reached = n_scored - n_reached,
    pct_not_reached = 100 * n_not_reached / n_scored,
    .groups = "drop"
  ) |>
  group_by(stage_order) |>
  mutate(complete_both_temperatures = all(n_not_reached == 0)) |>
  ungroup() |>
  mutate(
    pct_not_label = if_else(n_not_reached > 0,
                            paste0(round(pct_not_reached), "%"),
                            NA_character_)
  )

# ---- A. Count-scaled violin, close to Molly's mock-up ----------------------
p_violin <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 16.05, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = 16.05, linetype = "dotted",
             linewidth = 0.35, colour = "grey50") +
  geom_violin(
    data = filter(stage_events, event),
    aes(stage_short, first_yes_day, fill = treatment,
        group = interaction(stage_short, treatment)),
    position = position_dodge(width = 0.72),
    scale = "count",
    trim = TRUE,
    adjust = 0.85,
    alpha = 0.72,
    colour = "#222222",
    linewidth = 0.32,
    width = 0.68
  ) +
  geom_boxplot(
    data = filter(stage_events, event, stage_order %in%
                    stage_summary$stage_order[stage_summary$complete_both_temperatures]),
    aes(stage_short, first_yes_day,
        group = interaction(stage_short, treatment)),
    position = position_dodge(width = 0.72),
    width = 0.11,
    fill = NA,
    colour = "#222222",
    linewidth = 0.35,
    outlier.shape = NA
  ) +
  geom_point(
    data = filter(stage_events, event),
    aes(stage_short, first_yes_day, fill = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0,
      dodge.width = 0.72, seed = 39
    ),
    shape = 21,
    colour = "white",
    size = 1.25,
    stroke = 0.22,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(stage_events, !event),
    aes(stage_short, y = 17.05, colour = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0.06,
      dodge.width = 0.72, seed = 40
    ),
    shape = 4,
    size = 1.9,
    stroke = 0.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(stage_summary, !is.na(pct_not_label)),
    aes(stage_short, 17.7, label = pct_not_label, colour = treatment,
        group = treatment),
    position = position_dodge(width = 0.72),
    size = 2.45,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = PAL_TEMP, name = "Temperature",
                    labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_colour_manual(values = PAL_TEMP, guide = "none") +
  scale_y_continuous(
    breaks = c(1, 3, 5, 7, 9, 11, 13, 15, 17.05),
    labels = c("1", "3", "5", "7", "9", "11", "13", "15",
               "Not\nreached")
  ) +
  coord_cartesian(ylim = c(0.5, 18.15), clip = "on") +
  labs(
    title = "A. Count-scaled violins",
    subtitle = "Best for showing the distribution Molly asked for.",
    x = NULL,
    y = "Day first seen"
  ) +
  theme_pub(8) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(lineheight = 0.88, margin = margin(t = 4))
  )

save_fig(p_violin, "39a_stage_timing_violin_option", width = 180, height = 120)

# ---- B. Cumulative percent reached through time ----------------------------
score_days <- 0:15
stopifnot(!anyDuplicated(select(stage_events, stage_order, id)),
          all(stage_events$event | stage_events$last_scored_day == 15))

cumulative <- tidyr::expand_grid(
  stage_order = sort(unique(stage_events$stage_order)),
  treatment = levels(stage_events$treatment),
  day = score_days
) |>
  left_join(distinct(stage_events, stage_order, stage_short),
            by = "stage_order") |>
  left_join(
    stage_events |>
      select(stage_order, treatment, id, event, first_yes_day),
    by = c("stage_order", "treatment"),
    relationship = "many-to-many"
  ) |>
  group_by(stage_order, stage_short, treatment, day) |>
  summarise(
    n_scored = n(),
    n_reached = sum(event & first_yes_day <= day, na.rm = TRUE),
    pct_reached = 100 * n_reached / n_scored,
    .groups = "drop"
  ) |>
  mutate(treatment = factor(treatment, levels = c("28C", "31C")))
write_csv(cumulative, file.path(TBL_DIR, "39_stage_cumulative_percent.csv"))

# Keep missing daily scores distinct from observed absences. The main plot
# counts documented first observations, not imputed dates of biological onset.
score_audit <- map_dfr(unique(stage_events$stage_order), function(s) {
  meta <- filter(stage_events, stage_order == s)
  path <- if (s == 1) "microscope_physio_clean.rds" else "physio_clean.rds"
  d <- readRDS(file.path(DATA_PROC, path)) |>
    filter(wound == "yes", day %in% 0:15) |>
    select(id, treatment, day, y = all_of(unique(meta$trait)))
  stopifnot(!anyDuplicated(select(d, id, day)))
  d |> group_by(treatment, day) |>
    summarise(n_total = n(), n_scored_today = sum(!is.na(y)),
              n_missing_today = sum(is.na(y)),
              n_present_today = sum(y == 1, na.rm = TRUE), .groups = "drop") |>
    mutate(stage_order = s)
})
write_csv(score_audit, file.path(TBL_DIR, "39_stage_daily_score_audit.csv"))
stopifnot(all(cumulative$n_reached <= cumulative$n_scored))

p_cumulative <- ggplot(cumulative, aes(day, pct_reached, colour = treatment,
                                      linetype = treatment, shape = treatment)) +
  geom_step(linewidth = 0.75, direction = "hv") +
  geom_point(size = 1.5) +
  facet_wrap(~ stage_short, ncol = 3) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature",
                      labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_linetype_manual(values = c(`28C` = "solid", `31C` = "dashed"),
                        name = "Temperature", labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature", labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 10, 15),
                     limits = c(0, 15)) +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "Fewer heated fragments reached the final regrowth step by Day 15",
    subtitle = "All fragments stay in the denominator, including those that had not reached a stage",
    x = "Day",
    y = "Fragments first seen at this stage by each day",
    caption = "Tissue cover: 8 fragments per temperature in a separate photo experiment. Other panels: 12 per temperature.\nTwo heated fragments lacked a Day 1 tissue-cover score; both were first recorded as covered on Day 2."
  ) +
  theme_pub(10) +
  theme(legend.position = "bottom")

save_fig(p_cumulative, "39b_stage_timing_cumulative_option",
         width = 200, height = 145)

# ---- C. Stage-by-day heatmap -----------------------------------------------
heatmap_days <- c(1, 3, 5, 7, 10, 12, 15)

heatmap_dat <- tidyr::expand_grid(
  stage_order = sort(unique(stage_events$stage_order)),
  treatment = levels(stage_events$treatment),
  day = heatmap_days
) |>
  left_join(distinct(stage_events, stage_order, stage_short),
            by = "stage_order") |>
  left_join(
    stage_events |>
      select(stage_order, treatment, id, event, first_yes_day),
    by = c("stage_order", "treatment"),
    relationship = "many-to-many"
  ) |>
  group_by(stage_order, stage_short, treatment, day) |>
  summarise(
    n_scored = n(),
    n_reached = sum(event & first_yes_day <= day, na.rm = TRUE),
    pct_reached = round(100 * n_reached / n_scored),
    .groups = "drop"
  ) |>
  mutate(
    treatment_label = factor(
      dplyr::recode(treatment, `28C` = "28 °C", `31C` = "31 °C"),
      levels = c("28 °C", "31 °C")
    ),
    day = factor(day, levels = heatmap_days),
    stage_short = forcats::fct_rev(stage_short)
  )

p_heatmap <- ggplot(heatmap_dat, aes(day, stage_short, fill = pct_reached)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = paste0(pct_reached, "%")), size = 2.35,
            colour = "#1b1b1b") +
  facet_wrap(~ treatment_label, nrow = 1) +
  scale_fill_gradient(
    low = "#f4f4f4",
    high = "#0072B2",
    limits = c(0, 100),
    labels = \(x) paste0(x, "%"),
    name = "Reached"
  ) +
  labs(
    title = "C. Percent-reached heatmap",
    subtitle = "Best for a compact table-like summary.",
    x = "Day",
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )

save_fig(p_heatmap, "39c_stage_timing_heatmap_option",
         width = 180, height = 120)

# ---- Combined comparison sheet --------------------------------------------
comparison <- (p_violin / p_cumulative / p_heatmap) +
  plot_annotation(
    title = "Three ways to show healing and regeneration stage timing",
    subtitle = paste(
      "All panels use the same stage-timing table.",
      "\nChoose based on whether the priority is distribution, timing, or quick pattern recognition."
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )

save_fig(comparison, "39_stage_timing_visual_options",
         width = 190, height = 300)

cat("\n=== Stage-timing visual options ===\n")
cat("Wrote three individual options and one combined comparison sheet.\n")
