# =============================================================================
# Purpose: Four-panel exploratory trade-off summary for the team HTML.
#
#          This figure keeps the trade-off claim modest. It shows the one
#          pattern that could be a cost of heat tolerance (source C grows less
#          at 28C but loses less growth under heat), then checks whether
#          regeneration appears to come at the cost of growth or heat tolerance.
#
# Input:   output/tables/13_genet_emmeans.csv
#          output/tables/12_genet_treatment_effects.csv
#          data/processed/buoyant_weight_clean.rds
#          data/processed/physio_clean.rds
# Output:  output/tables/35_tradeoff_join_audit.csv
#          output/tables/35_tradeoff_summary.csv
#          figures/35_tradeoff_summary.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Growth means and growth heat loss by source --------------------------
growth_means <- read_csv(file.path(TBL_DIR, "13_genet_emmeans.csv"),
                         show_col_types = FALSE) |>
  filter(response == "Growth (% mass change)") |>
  mutate(
    source_label = str_to_upper(thicket),
    treatment_label = factor(dplyr::recode(treatment,
                                           `28C` = "28 °C",
                                           `31C` = "31 °C"),
                             levels = c("28 °C", "31 °C")),
    temp_x = if_else(treatment == "28C", 1, 2)
  )

growth_tradeoff <- growth_means |>
  select(thicket, source_label, treatment, mean, se, n) |>
  pivot_wider(names_from = treatment,
              values_from = c(mean, se, n),
              names_sep = "_") |>
  mutate(
    ambient_growth_pct = mean_28C,
    heated_growth_pct = mean_31C,
    growth_heat_drop_pct = mean_28C - mean_31C,
    growth_heat_drop_se = sqrt(se_28C^2 + se_31C^2)
  )

# ---- Same-fragment growth and Day-15 regeneration -------------------------
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

write_csv(join_audit, file.path(TBL_DIR, "35_tradeoff_join_audit.csv"))

if (any(!growth_morph$treatment_match) || any(!growth_morph$thicket_match)) {
  stop("Growth and morphology disagree on treatment or source-patch labels.")
}

wounded_heat <- growth_morph |>
  filter(wound_morph == "yes", treatment_morph == "31C") |>
  mutate(
    thicket = thicket_morph,
    source_label = str_to_upper(thicket),
    new_corallites_label = factor(
      if_else(new_corallites_on_tip == 1,
              "New corallites by Day 15",
              "No new corallites"),
      levels = c("No new corallites", "New corallites by Day 15")
    )
  )

growth_regen_summary <- wounded_heat |>
  group_by(new_corallites_label) |>
  summarise(
    n = n(),
    growth_mean = mean(pct_growth, na.rm = TRUE),
    growth_se = sd(pct_growth, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ---- Independent source-level heat-tolerance axis -------------------------
# This heat-penalty score excludes regeneration so panel D does not compare a
# response with itself.
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

source_regen_heat <- growth_morph |>
  filter(wound_morph == "yes", treatment_morph == "31C") |>
  mutate(thicket = thicket_morph,
         source_label = str_to_upper(thicket)) |>
  group_by(thicket, source_label) |>
  summarise(
    n_heated_wounded = n(),
    growth_mean_31C_wounded = mean(pct_growth, na.rm = TRUE),
    new_corallites_pct = 100 * mean(new_corallites_on_tip, na.rm = TRUE),
    tip_extension_pct = 100 * mean(tip_extension, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(phys_growth_heat_penalty, by = "thicket")

tradeoff_summary <- growth_tradeoff |>
  select(thicket, source_label, ambient_growth_pct, heated_growth_pct,
         growth_heat_drop_pct) |>
  left_join(source_regen_heat, by = c("thicket", "source_label")) |>
  arrange(source_label)

write_csv(tradeoff_summary, file.path(TBL_DIR, "35_tradeoff_summary.csv"))

# ---- Plot panels -----------------------------------------------------------
label_31 <- growth_means |>
  filter(treatment == "31C") |>
  mutate(label_x = temp_x + 0.08)

p_growth_slope <- ggplot(growth_means,
                         aes(temp_x, mean,
                             colour = thicket, group = thicket)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se,
                      shape = thicket),
                  size = 0.4) +
  geom_text(data = label_31,
            aes(label_x, mean, label = source_label, colour = thicket),
            inherit.aes = FALSE, hjust = 0, size = 3.0, fontface = "bold") +
  scale_colour_manual(values = PAL_GENO, guide = "none") +
  scale_shape_manual(values = SHP_GENO, guide = "none") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("28 °C", "31 °C"),
                     limits = c(0.9, 2.22)) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL,
       y = "Skeletal growth\n(% mass gain)",
       title = "A. Source patch C loses less growth under heat") +
  theme_pub(8.5) +
  theme(plot.title = element_text(size = 8.5),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(8, 18, 8, 8, "pt"))

p_growth_tradeoff <- ggplot(growth_tradeoff,
                            aes(ambient_growth_pct, growth_heat_drop_pct,
                                colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.4,
           fill = "#009E73", alpha = 0.035) +
  annotate("text", x = 5.03, y = 0.18, label = "steadier\nunder heat",
           size = 2.5, colour = "grey30", hjust = 0) +
  annotate("text", x = 7.00, y = 3.33,
           label = "faster at 28 °C,\nbigger heat drop",
           size = 2.5, colour = "grey30", hjust = 1, vjust = 1) +
  geom_errorbar(aes(ymin = growth_heat_drop_pct - growth_heat_drop_se,
                    ymax = growth_heat_drop_pct + growth_heat_drop_se),
                width = 0, linewidth = 0.3, alpha = 0.65) +
  geom_errorbar(aes(xmin = ambient_growth_pct - se_28C,
                    xmax = ambient_growth_pct + se_28C),
                orientation = "y", width = 0, linewidth = 0.3, alpha = 0.65) +
  geom_point(size = 3.1, alpha = 0.95) +
  geom_text(aes(label = source_label),
            nudge_x = 0.07, nudge_y = 0.03,
            hjust = 0, size = 3.0, fontface = "bold") +
  scale_colour_manual(values = PAL_GENO, guide = "none") +
  scale_shape_manual(values = SHP_GENO, guide = "none") +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  scale_y_continuous(labels = label_number(suffix = " pts")) +
  coord_cartesian(xlim = c(5.0, 7.05), ylim = c(0, 3.45),
                  clip = "off") +
  labs(x = "Growth at 28 °C",
       y = "Growth lost under heat\n(28 °C - 31 °C)",
       title = "B. Possible growth/tolerance trade-off") +
  theme_pub(8.5) +
  theme(plot.title = element_text(size = 8.5))

p_growth_regen <- ggplot(wounded_heat,
                         aes(new_corallites_label, pct_growth,
                             colour = thicket, shape = thicket)) +
  geom_point(position = position_jitter(width = 0.12, height = 0),
             size = 2.6, alpha = 0.9) +
  stat_summary(aes(group = new_corallites_label),
               fun = mean, geom = "point", shape = 95,
               size = 8, colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Source patch",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_discrete(labels = c("No new\ncorallites",
                              "New corallites\nby Day 15")) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  labs(x = NULL,
       y = "Skeletal growth\nat 31 °C",
       title = "C. Regrowth does not mean lower growth") +
  theme_pub(8.5) +
  theme(plot.title = element_text(size = 8.5),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

p_heat_regen <- ggplot(source_regen_heat,
                       aes(phys_growth_heat_penalty, new_corallites_pct,
                           colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0.55, ymin = 50, ymax = Inf,
           fill = "#009E73", alpha = 0.035) +
  annotate("text", x = 0.29, y = 94,
           label = "lower heat penalty,\nmore regrowth",
           hjust = 0, vjust = 1, size = 2.5, colour = "grey30") +
  annotate("text", x = 1.02, y = 6,
           label = "higher heat penalty,\nless regrowth",
           hjust = 1, size = 2.5, colour = "grey30") +
  geom_point(size = 3.3, alpha = 0.95) +
  geom_text(aes(label = source_label),
            nudge_x = 0.025, nudge_y = 3,
            hjust = 0, size = 3.0, fontface = "bold") +
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
       title = "D. No tolerance/regrowth trade-off") +
  theme_pub(8.5) +
  theme(plot.title = element_text(size = 8.5))

p_tradeoff <- (p_growth_slope + p_growth_tradeoff) /
  (p_growth_regen + p_heat_regen) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Trade-off checks: one possible growth cost, no regeneration cost",
    subtitle = str_wrap(
      "Source patch C grew less at 28 °C but lost less growth under heat. The regeneration checks point the other way: under heat, fragments and source patches with more regrowth did not show lower growth or higher heat penalty.",
      width = 115
    ),
    caption = str_wrap(
      "All panels are exploratory. Points are coral fragments in panel C and source-patch means in panels A, B, and D. Vertical bars in A and horizontal/vertical bars in B are standard errors (SE). Black ticks in C are group means. A, C, and D are source-patch labels, not confirmed genetic individuals.",
      width = 115
    )
  ) &
  theme_pub(8.5) &
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"),
        legend.position = "bottom")

save_fig(p_tradeoff, "35_tradeoff_summary", width = 210, height = 175)

cat("\n=== Trade-off summary ===\n")
print(tradeoff_summary |>
        select(source_label, ambient_growth_pct, growth_heat_drop_pct,
               phys_growth_heat_penalty, new_corallites_pct) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))
