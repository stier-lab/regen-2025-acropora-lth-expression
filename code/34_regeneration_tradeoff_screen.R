# =============================================================================
# Purpose: Exploratory screens for growth-regeneration and
#          heat-tolerance-regeneration trade-offs.
#
#          These are descriptive plots, not formal trade-off tests. The
#          growth-regeneration panel uses wounded fragments with both growth and
#          Day-15 morphology. The heat-tolerance-regeneration panel compares
#          source-level physiology/growth heat penalties with the percentage of
#          heated wounded fragments that budded new corallites by Day 15.
#
# Input:   data/processed/buoyant_weight_clean.rds
#          data/processed/physio_clean.rds
#          output/tables/12_genet_treatment_effects.csv
# Output:  output/tables/34_growth_morphology_join_audit.csv
#          output/tables/34_growth_regeneration_summary.csv
#          output/tables/34_heat_tolerance_regeneration_screen.csv
#          figures/34_regeneration_tradeoff_screens.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Same-fragment growth and regeneration table --------------------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  select(id,
         treatment_growth = treatment,
         wound_growth = wound,
         thicket_growth = thicket,
         pct_growth)

phys15 <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(day == 15) |>
  select(id,
         treatment_morph = treatment,
         wound_morph = wound,
         thicket_morph = thicket,
         new_corallites_on_tip,
         tip_extension,
         tip_exist)

growth_morph <- inner_join(bw, phys15, by = "id") |>
  mutate(
    treatment_match = as.character(treatment_growth) ==
      as.character(treatment_morph),
    wound_match = as.character(wound_growth) == as.character(wound_morph),
    thicket_match = thicket_growth == thicket_morph
  )

join_audit <- tibble(
  check = c(
    "growth rows",
    "Day-15 morphology rows",
    "joined rows",
    "unique joined IDs",
    "treatment-label mismatches",
    "wound-label mismatches",
    "source-patch-label mismatches"
  ),
  value = c(
    nrow(bw),
    nrow(phys15),
    nrow(growth_morph),
    n_distinct(growth_morph$id),
    sum(!growth_morph$treatment_match),
    sum(!growth_morph$wound_match),
    sum(!growth_morph$thicket_match)
  ),
  notes = c(
    "One growth row per physiology fragment.",
    "One Day-15 morphology row per physiology fragment.",
    "Joined by coral fragment ID.",
    "Duplicate IDs would make same-fragment trade-off plots invalid.",
    "Should be zero before plotting.",
    "One known mismatch means morphology wound labels define the wounded-only regeneration screen.",
    "Should be zero before plotting."
  )
)

write_csv(join_audit,
          file.path(TBL_DIR, "34_growth_morphology_join_audit.csv"))

if (any(!growth_morph$treatment_match) || any(!growth_morph$thicket_match)) {
  stop("Growth and morphology disagree on treatment or source-patch labels.")
}

wounded <- growth_morph |>
  filter(wound_morph == "yes") |>
  mutate(
    treatment = factor(treatment_morph, levels = c("28C", "31C")),
    treatment_label = factor(dplyr::recode(as.character(treatment),
                                           `28C` = "28 °C",
                                           `31C` = "31 °C"),
                             levels = c("28 °C", "31 °C")),
    thicket = thicket_morph,
    source_label = str_to_upper(thicket),
    new_corallites_label = factor(
      if_else(new_corallites_on_tip == 1,
              "Budded by Day 15",
              "No new corallites"),
      levels = c("No new corallites", "Budded by Day 15")
    )
  )

growth_regen_summary <- wounded |>
  group_by(treatment, new_corallites_label) |>
  summarise(
    n = n(),
    growth_mean = mean(pct_growth, na.rm = TRUE),
    growth_se = sd(pct_growth, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  complete(treatment, new_corallites_label,
           fill = list(n = 0, growth_mean = NA_real_, growth_se = NA_real_))

source_regen_summary <- wounded |>
  group_by(treatment, thicket, source_label) |>
  summarise(
    n = n(),
    growth_mean = mean(pct_growth, na.rm = TRUE),
    new_corallites_pct = 100 * mean(new_corallites_on_tip, na.rm = TRUE),
    tip_extension_pct = 100 * mean(tip_extension, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(growth_regen_summary,
          file.path(TBL_DIR, "34_growth_regeneration_summary.csv"))

# ---- Independent source-level heat-tolerance axis -------------------------
# Use only whole-coral condition/growth responses here. Do not include
# regeneration in this heat-tolerance score, or the x and y axes would not be
# independent.
phys_growth_heat_penalty <- read_csv(
  file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
  show_col_types = FALSE
) |>
  filter(response %in% c("pam_fvfm", "color_dscale", "growth_pct",
                         "log_zoox_density"),
         is.finite(estimate)) |>
  group_by(response, thicket) |>
  summarise(estimate = mean(estimate, na.rm = TRUE),
            .groups = "drop") |>
  group_by(response) |>
  mutate(scaled_heat_penalty =
           estimate / max(abs(estimate), na.rm = TRUE)) |>
  ungroup() |>
  group_by(thicket) |>
  summarise(
    phys_growth_heat_penalty = mean(scaled_heat_penalty, na.rm = TRUE),
    n_responses = n(),
    .groups = "drop"
  )

heat_regen <- source_regen_summary |>
  filter(treatment == "31C") |>
  left_join(phys_growth_heat_penalty, by = "thicket") |>
  arrange(phys_growth_heat_penalty)

write_csv(heat_regen,
          file.path(TBL_DIR, "34_heat_tolerance_regeneration_screen.csv"))

# ---- Plot -----------------------------------------------------------------
p_growth_regen <- ggplot(wounded,
                         aes(new_corallites_label, pct_growth,
                             colour = thicket, shape = thicket)) +
  geom_point(position = position_jitter(width = 0.12, height = 0),
             size = 2.5, alpha = 0.9) +
  stat_summary(aes(group = new_corallites_label),
               fun = mean, geom = "point", shape = 95,
               size = 8, colour = "black", show.legend = FALSE) +
  facet_wrap(~ treatment_label, nrow = 1) +
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Source patch",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_discrete(drop = FALSE,
                   labels = c("No new\ncorallites", "New corallites\nby Day 15")) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(x = NULL,
       y = "Skeletal growth\n(% mass gain)",
       title = "A. Regrowth does not come with lower growth") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

p_heat_regen <- ggplot(heat_regen,
                       aes(phys_growth_heat_penalty, new_corallites_pct,
                           colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0.55, ymin = 50, ymax = Inf,
           fill = "#009E73", alpha = 0.035) +
  annotate("text", x = 0.30, y = 91,
           label = "lower heat penalty,\nmore regrowth",
           hjust = 0, size = 2.7, colour = "grey30") +
  annotate("text", x = 0.99, y = 8,
           label = "higher heat penalty,\nless regrowth",
           hjust = 1, size = 2.7, colour = "grey30") +
  geom_point(size = 3.7, alpha = 0.95) +
  geom_text(aes(label = source_label),
            nudge_x = 0.025, nudge_y = 3,
            hjust = 0, size = 3.3, fontface = "bold") +
  scale_colour_manual(values = PAL_GENO, guide = "none") +
  scale_shape_manual(values = SHP_GENO, guide = "none") +
  scale_x_continuous(
    limits = c(0.22, 1.05),
    breaks = c(0.25, 0.5, 0.75, 1),
    labels = c("lower\npenalty", "0.5", "0.75", "higher\npenalty")
  ) +
  scale_y_continuous(limits = c(-4, 104),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = label_number(suffix = "%")) +
  labs(x = "Whole-coral heat penalty\n(condition + growth only)",
       y = "Heated wounded fragments with\nnew corallites by Day 15",
       title = "B. No heat-tolerance/regrowth trade-off") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9))

p <- p_growth_regen + p_heat_regen +
  plot_layout(widths = c(1.18, 1)) +
  plot_annotation(
    title = "Regeneration trade-off screens",
    subtitle = str_wrap(
      "The current data do not show a growth-regeneration or heat-tolerance-regeneration trade-off. Regenerated heated fragments tended to have at least as much growth, and source patch C was both less heat-penalized and more likely to bud new corallites under heat.",
      width = 112
    ),
    caption = str_wrap(
      "Panel A uses wounded fragments with same-fragment growth and Day 15 healing/regrowth scores. Black ticks mark group means. Panel B uses source-patch means; the heat-penalty axis excludes regeneration traits. A, C, and D are source-patch labels, not confirmed genetic individuals.",
      width = 112
    )
  ) &
  theme_pub(9) &
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"),
        legend.position = "bottom")

save_fig(p, "34_regeneration_tradeoff_screens", width = 200, height = 112)

cat("\n=== Growth-regeneration summary ===\n")
print(growth_regen_summary |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))

cat("\n=== Heat-tolerance-regeneration screen ===\n")
print(heat_regen |>
        select(source_label, phys_growth_heat_penalty, new_corallites_pct,
               tip_extension_pct, growth_mean) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))
