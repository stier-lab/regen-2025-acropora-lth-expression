# =============================================================================
# Purpose: Exploratory growth trade-off screen.
#
#          Ask whether the source that is steadiest under heat also grows more
#          slowly at 28C. This is a descriptive source-level screen, not a formal
#          test of a trade-off, because there are only three source patches.
#
# Input:   output/tables/13_genet_emmeans.csv
# Output:  output/tables/33_growth_heat_tradeoff.csv
#          figures/33_growth_heat_tradeoff.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

growth_means <- read_csv(file.path(TBL_DIR, "13_genet_emmeans.csv"),
                         show_col_types = FALSE) |>
  filter(response == "Growth (% mass change)") |>
  mutate(
    source_label = str_to_upper(thicket),
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
    growth_heat_drop_se = sqrt(se_28C^2 + se_31C^2),
    interpretation = case_when(
      source_label == "C" ~ "lower growth at 28 °C, smaller heat drop",
      TRUE ~ "higher growth at 28 °C, larger heat drop"
    )
  ) |>
  arrange(source_label)

write_csv(growth_tradeoff,
          file.path(TBL_DIR, "33_growth_heat_tradeoff.csv"))

label_positions <- growth_means |>
  filter(treatment == "31C") |>
  mutate(label_x = temp_x + 0.08)

p_slope <- ggplot(growth_means,
                  aes(temp_x, mean, colour = thicket, group = thicket)) +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se,
                      shape = thicket),
                  size = 0.45) +
  geom_text(data = label_positions,
            aes(label_x, mean, label = source_label, colour = thicket),
            inherit.aes = FALSE, hjust = 0, size = 3.2, fontface = "bold") +
  scale_colour_manual(values = PAL_GENO, guide = "none") +
  scale_shape_manual(values = SHP_GENO, guide = "none") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("28 °C", "31 °C"),
                     limits = c(0.9, 2.24)) +
  scale_y_continuous(labels = label_number(suffix = "%")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL,
       y = "Skeletal growth\n(% mass gain)",
       title = "A. Growth drops less in source patch C") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(8, 20, 8, 8, "pt"))

p_trade <- ggplot(growth_tradeoff,
                  aes(ambient_growth_pct, growth_heat_drop_pct,
                      colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.4,
           fill = "#009E73", alpha = 0.035) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.4, ymax = Inf,
           fill = "#D55E00", alpha = 0.04) +
  geom_errorbar(aes(ymin = growth_heat_drop_pct - growth_heat_drop_se,
                    ymax = growth_heat_drop_pct + growth_heat_drop_se),
                width = 0, linewidth = 0.35, alpha = 0.65) +
  geom_errorbar(aes(xmin = ambient_growth_pct - se_28C,
                    xmax = ambient_growth_pct + se_28C),
                orientation = "y", width = 0, linewidth = 0.35, alpha = 0.65) +
  geom_point(size = 3.5, alpha = 0.95) +
  geom_text(aes(label = source_label),
            nudge_x = 0.08, nudge_y = 0.04,
            hjust = 0, size = 3.2, fontface = "bold") +
  annotate("text", x = 5.04, y = 0.20, label = "steadier\nunder heat",
           size = 2.7, colour = "grey30", hjust = 0) +
  annotate("text", x = 6.98, y = 3.35, label = "faster at 28 °C,\nbigger heat drop",
           size = 2.7, colour = "grey30", hjust = 1, vjust = 1) +
  scale_colour_manual(values = PAL_GENO, guide = "none") +
  scale_shape_manual(values = SHP_GENO, guide = "none") +
  scale_x_continuous(labels = label_number(suffix = "%")) +
  scale_y_continuous(labels = label_number(suffix = " pts")) +
  coord_cartesian(xlim = c(5.0, 7.05), ylim = c(0, 3.45),
                  clip = "off") +
  labs(x = "Growth at 28 °C",
       y = "Growth lost under heat\n(28 °C - 31 °C)",
       title = "B. Possible growth-tolerance trade-off") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9))

p_tradeoff <- p_slope + p_trade +
  plot_layout(widths = c(0.95, 1.25)) +
  plot_annotation(
    title = "A possible growth trade-off, not yet a strong claim",
    subtitle = str_wrap(
      "Source patch C grew less at 28 °C than A or D, but it lost much less growth under heat. With only three source patches, this is a pattern to flag, not proof of an evolutionary or energetic trade-off.",
      width = 105
    ),
    caption = str_wrap(
      "Points are source-patch means from the end-of-experiment growth table; error bars show +/- 1 SE. A, C, and D are source-patch labels, not confirmed genetic individuals.",
      width = 105
    )
  ) &
  theme_pub(9) &
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"))

save_fig(p_tradeoff, "33_growth_heat_tradeoff", width = 190, height = 105)

cat("\n=== Growth heat trade-off screen ===\n")
print(growth_tradeoff |>
        select(source_label, ambient_growth_pct, heated_growth_pct,
               growth_heat_drop_pct, interpretation) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))
