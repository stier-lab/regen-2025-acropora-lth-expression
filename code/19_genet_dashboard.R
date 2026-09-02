# =============================================================================
# Purpose: Integrative cross-response genet dashboard. For each of the 3
#          genets (a, c, d), compute a standardized "heat sensitivity" effect
#          across each usable response variable (PAM, color, growth, symbionts,
#          and morphology time-to-onset traits). Build a forest plot of all
#          effect sizes to identify which genet is most/least thermally
#          sensitive across every dimension.
#
#          Also produces a composite heat-sensitivity ranking per genet: the
#          mean standardized heat effect across finite response-level effects.
#
# What & why: synthesis script. Other analyses test whether heat affects a
#   given response; this one addresses the cross-cutting question of which
#   genet best tolerates heat overall. To put effects measured in different
#   units on one axis (Fv/Fm, a colour score, growth,
#   symbiont counts, and wound-healing hazard ratios), each effect is rescaled
#   WITHIN its own response so the largest-magnitude effect = 1 (the "z"
#   column — a row-max standardization, not a statistical z-score). Positive
#   means the genet's phenotype is worse at 31 °C than 28 °C, i.e. more heat-
#   sensitive. Averaging these standardized values per genet yields a single
#   resilience ranking; in this dataset genet C is the most resilient. Nothing
#   here fits a new model — every number is a re-summary of upstream tables.
# Input:   output/tables/12_genet_treatment_effects.csv   (continuous responses)
#          output/tables/14_cox_hazard_ratios.csv         (per-genet HR for KM)
#          output/tables/15_genet_pca_displacement.csv    (multivariate)
# Output:  figures/19_genet_dashboard.{pdf,png}            — forest plot
#          figures/19b_genet_resilience_ranking.{pdf,png}  — composite score
#          figures/19d_wound_healing_heat_penalties.{pdf,png}
#                                                               — morphology-only
#                                                                 genet penalty
#          output/tables/19_genet_resilience_summary.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (TBL_DIR, FIG_DIR), theme_pub(),
# save_fig(), and the genet colour palette PAL_GENO used below.
source(here::here("code", "00_setup.R"))

# ---- Load upstream effect tables ------------------------------------------
# Three sources, one per analysis domain: per-genet continuous LMM contrasts
# (script 12), wound-healing Cox hazard ratios (script 14), and the
# multivariate per-genet PCA displacement (script 15).
cont <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                 show_col_types = FALSE)
cox  <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                 show_col_types = FALSE)
pca  <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"),
                 show_col_types = FALSE)

# ---- Standardize continuous-response effects ------------------------------
# Continuous responses (from script 12) give estimate = mean(28C) - mean(31C)
# at the response-specific endpoint per genet x wound, with SE:
# PAM/color at experimental Day 14, symbionts at the final biopsy, and growth
# as a start-to-end response. We collapse over wound (both wounded and
# unwounded corals contribute), then row-max scale within each response.
cont_eff <- cont |>
  filter(response %in% c("pam_fvfm", "color_dscale", "growth_pct",
                         "log_zoox_density"),
         is.finite(estimate)) |>
  group_by(response, thicket) |>
  summarise(estimate = mean(estimate, na.rm = TRUE),
            se = mean(SE, na.rm = TRUE),
            .groups = "drop") |>
  # Re-scale per-response so values are comparable
  group_by(response) |>
  mutate(z = estimate / max(abs(estimate), na.rm = TRUE)) |>
  ungroup() |>
  mutate(metric = "Delta phenotype (28C - 31C, response endpoint)")

# ---- Standardize Cox hazard ratios (per genet) ----------------------------
# HR < 1 means 31C delays/prevents trait expression — i.e., heat sensitivity.
# Convert to log(HR) so negative = heat-impaired, positive = heat-promoted.
cox_per_genet <- cox |>
  filter(grepl("^genet=", scope), is.finite(HR_31_vs28),
         HR_31_vs28 > 0, HR_31_vs28 < Inf) |>
  mutate(thicket = sub("\\s+first-observed approximation$", "",
                       sub("^genet=", "", scope)),
         logHR = log(HR_31_vs28),
         response = paste0("morph_", trait))

cox_eff <- cox_per_genet |>
  group_by(response) |>
  mutate(hr_scale = max(abs(logHR), na.rm = TRUE),
         z = if_else(is.finite(hr_scale) & hr_scale > 0,
                     -logHR / hr_scale, NA_real_)) |>
  # Sign flipped: more negative HR (delayed healing) → more positive z (more impaired)
  ungroup() |>
  filter(is.finite(z)) |>
  select(response, thicket, estimate = logHR, z) |>
  mutate(metric = "-log(HR 31C/28C) for healing trait onset",
         se = NA_real_)

# ---- Combine and forest plot ---------------------------------------------
all_eff <- bind_rows(
  cont_eff |> select(response, thicket, estimate, se, z, metric),
  cox_eff
) |>
  mutate(
    response_label = case_when(
      response == "pam_fvfm"          ~ "PAM Fv/Fm\nDay-14 contrast",
      response == "color_dscale"      ~ "Color (D-scale)\nDay-14 contrast",
      response == "growth_pct"        ~ "Growth\nstart-to-end contrast",
      response == "log_zoox_density"  ~ "log symbionts per cm2\nfinal biopsy contrast",
      grepl("^morph_", response)      ~ paste0(str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))), "\ntime-to-onset HR"),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_pct","log_zoox_density")
        ~ "Physiology",
      grepl("hole|polyp|smoothed|pigment", response_label, ignore.case = TRUE)
        ~ "Wound closure",
      grepl("tip|corallite", response_label, ignore.case = TRUE)
        ~ "Regeneration",
      TRUE
        ~ "Other"
    )
  )

# Sort by domain then by response_label for the y-axis (so rows group by domain)
all_eff <- all_eff |>
  mutate(response_label = factor(response_label,
                                  levels = unique(response_label[order(domain, response_label)])))

# Forest plot: each point is one genet's standardized heat effect on one
# response, faceted by domain. Dashed line at x = 0 marks "no heat effect".
dash_mean <- all_eff |>
  group_by(thicket) |>
  summarise(mean_penalty = mean(z, na.rm = TRUE),
            n_responses = n(),
            .groups = "drop") |>
  mutate(genet_label = factor(str_to_upper(thicket), levels = c("A", "D", "C")),
         label_x = mean_penalty + 0.035,
         label_hjust = 0)

p_dash_mean <- ggplot(dash_mean,
                      aes(mean_penalty, genet_label, fill = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.035) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.045) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.25,
           alpha = 0.88) +
  geom_text(aes(x = label_x,
                label = sprintf("%.2f", mean_penalty),
                hjust = label_hjust),
            size = 3, colour = "grey20") +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.3, 0.6),
                     labels = c("little/no\npenalty", "", "larger\npenalty")) +
  coord_cartesian(xlim = c(-0.08, 0.68), clip = "off") +
  labs(x = "Average heat penalty",
       y = NULL,
       title = "A. Average across responses") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        legend.position = "none")

p_dash_detail <- ggplot(all_eff,
                        aes(z, response_label,
                            colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.025) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.035) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_point(size = 3.0, alpha = 0.95) +
  scale_colour_manual(values = PAL_GENO, name = "Genet",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Genet",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("better/faster\nat 31C", "",
                                "no\npenalty", "",
                                "worse/slower\nat 31C")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Heat effect within each response",
       y = NULL,
       title = "B. Response-level heat effects") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"))

p_dash <- p_dash_mean + p_dash_detail +
  plot_layout(widths = c(0.82, 2.35), guides = "collect") +
  plot_annotation(
    title = "Genet C is least harmed by heat",
    subtitle = str_wrap(
      "Each response is first converted to a 31C vs 28C heat effect and put on a common within-response scale. C stays close to little/no penalty, while A and D show larger losses or delays.",
      width = 105
    ),
    caption = str_wrap(
      "Positive values mean worse physiology or slower milestone onset at 31C; negative values mean higher physiology or faster onset at 31C. Physiology uses endpoint contrasts averaged across wound state; morphology uses wounded-only Cox time-to-onset effects.",
      width = 105
    )
  ) &
  theme_pub(9) &
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

save_fig(p_dash, "19_genet_dashboard", width = 200, height = 185)

# ---- Composite thermal resilience score ----------------------------------
# Per genet: mean standardized heat sensitivity across finite response effects.
# Lower composite = smaller heat penalty and therefore greater resilience.
resilience <- all_eff |>
  group_by(thicket) |>
  summarise(
    mean_sensitivity   = mean(z, na.rm = TRUE),
    median_sensitivity = median(z, na.rm = TRUE),
    n_responses        = n(),
    .groups = "drop"
  ) |>
  left_join(pca |> select(thicket, pca_displacement = displacement),
            by = "thicket") |>
  mutate(rank_overall = rank(mean_sensitivity))

composite_n_note <- resilience |>
  arrange(thicket) |>
  transmute(label = paste0(str_to_upper(thicket), " n=", n_responses)) |>
  pull(label) |>
  paste(collapse = ", ")

# This summary table is the primary output: script 20 reads it (and its
# rank_overall column) back in as the "Genet resilience" rows.
write_csv(resilience, file.path(TBL_DIR, "19_genet_resilience_summary.csv"))

# Bar chart of the composite ranking; each bar is annotated with that genet's
# multivariate PCA displacement as an independent cross-check on the ordering.
p_rank <- ggplot(resilience,
                  aes(reorder(thicket, -mean_sensitivity),
                      mean_sensitivity, fill = thicket)) +
  geom_col(width = 0.55, alpha = 0.85, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("n=%d\nPCA = %.2f",
                                n_responses, pca_displacement)),
            vjust = -0.25, lineheight = 0.95, size = 3, colour = "grey20") +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
  labs(x = "Genet",
       y = "Mean relative heat penalty\n(row-scaled, unitless)",
       title = "Composite of standardized heat effects",
       subtitle = str_wrap(
         "Lower bars mean smaller heat penalties. Inputs are modeled contrasts or time-to-onset effects, not same-day raw values.",
         width = 55
       ),
       caption = str_wrap(
         paste0("Finite effects averaged per genet: ", composite_n_note,
                ". Missing or non-finite hazard ratios are omitted. PCA displacement is a separate descriptive check using centered/scaled final physiology values."),
         width = 70
       )) +
  coord_cartesian(ylim = c(0, max(resilience$mean_sensitivity) * 1.35)) +
  theme_pub(10)

save_fig(p_rank, "19b_genet_resilience_ranking", width = 130, height = 110)

# ---- Decomposed dashboard: heat-only vs heat-while-wounded ---------------
# The composite above pools across wound state. To answer "is genet C
# resilient to heat per se, or only to heat-while-wounded?" we split the
# continuous-response standardized effect into two scopes. Morphology Cox
# effects remain wounded-only by design.
cont_by_wound <- cont |>
  # Drop morphology rows — they have wound=NA (wounded-only by construction)
  # and are captured separately in cox_by_wound.
  filter(!is.na(wound),
         response %in% c("pam_fvfm", "color_dscale", "growth_pct",
                          "log_zoox_density")) |>
  group_by(response, thicket, wound) |>
  summarise(estimate = mean(estimate, na.rm = TRUE),
            se = mean(SE, na.rm = TRUE),
            .groups = "drop") |>
  group_by(response, wound) |>
  mutate(z = estimate / max(abs(estimate), na.rm = TRUE)) |>
  ungroup() |>
  mutate(scope = if_else(wound == "yes",
                         "heat while wounded",
                         "heat only (unwounded)"))

# Cox results are wounded-only by design — flag them as the third scope
cox_by_wound <- cox_eff |>
  select(response, thicket, z) |>
  mutate(scope = "heat while wounded",
         wound = "yes")

decomp <- bind_rows(
  cont_by_wound |> select(response, thicket, wound, scope, z),
  cox_by_wound  |> select(response, thicket, wound, scope, z)
) |>
  mutate(
    response_label = case_when(
      response == "pam_fvfm"         ~ "PAM Fv/Fm\nDay-14 contrast",
      response == "color_dscale"     ~ "Color (D)\nDay-14 contrast",
      response == "growth_pct"       ~ "Growth\nstart-to-end contrast",
      response == "log_zoox_density" ~ "log symbionts per cm2\nfinal biopsy contrast",
      grepl("^morph_", response)     ~ paste0(str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))), "\ntime-to-onset HR"),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_pct","log_zoox_density")
        ~ "Physiology",
      grepl("hole|polyp|smoothed|pigment", response_label, ignore.case = TRUE)
        ~ "Wound closure",
      grepl("tip|corallite", response_label, ignore.case = TRUE)
        ~ "Regeneration",
      TRUE ~ "Other"
    )
  )

# Same forest plot as p_dash, but now faceted domain (rows) × scope (columns)
# so heat-only and heat-while-wounded sensitivity can be compared side by side.
p_decomp <- ggplot(decomp,
                    aes(z, response_label,
                        colour = thicket, shape = thicket)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_grid(domain ~ scope, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = PAL_GENO, name = "Genet") +
  scale_shape_manual(values = c(a = 16, c = 17, d = 15), name = "Genet") +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1", "-0.5", "0", "0.5", "1\nrow max")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Relative heat penalty within response and scope (unitless; max abs(effect) = 1)",
       y = NULL,
       title = "Decomposed resilience: heat-only vs heat-while-wounded",
       subtitle = str_wrap(
         "Physiology panels use endpoint 31 °C vs 28 °C contrasts split by wound state; morphology panels use wounded-only time-to-onset hazard ratios.",
         width = 115
       ),
       caption = str_wrap(
         "Values are row-max scaled within each response and scope. Positive values mean a larger heat penalty; negative values mean the 31 °C group reached a milestone faster or had a higher endpoint value.",
         width = 115
       )) +
  theme_pub(9) +
  theme(panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold"))

save_fig(p_decomp, "19c_decomposed_resilience", width = 200, height = 175)

# ---- Focused wound-healing heat penalties --------------------------------
# Same row-max-scaled Cox effects as the dashboard, filtered to the wounded-only
# morphology milestones so the genet pattern is visible without the physiology
# rows. Positive values mean the milestone was delayed under 31C.
morph_penalty <- cox_eff |>
  mutate(
    trait = sub("^morph_", "", response),
    trait_label = recode(
      trait,
      axial_polyp_formation = "Axial polyp formed",
      wound_smoothed = "Wound smoothed",
      pigment_over_wound = "Pigment over wound",
      tip_exist = "Tip visible",
      tip_extension = "Tip extension",
      new_corallites_on_tip = "New corallites on tip"
    ),
    trait_label = factor(
      trait_label,
      levels = rev(c(
        "Wound smoothed",
        "Pigment over wound",
        "Axial polyp formed",
        "Tip visible",
        "Tip extension",
        "New corallites on tip"
      ))
    )
  )

morph_mean <- morph_penalty |>
  group_by(thicket) |>
  summarise(mean_penalty = mean(z, na.rm = TRUE),
            n_traits = n(),
            .groups = "drop") |>
  mutate(genet_label = factor(str_to_upper(thicket), levels = c("A", "D", "C")),
         label_x = if_else(mean_penalty >= 0,
                           mean_penalty + 0.025,
                           mean_penalty / 2),
         label_hjust = if_else(mean_penalty >= 0, 0, 0.5))

p_morph_mean <- ggplot(morph_mean,
                       aes(mean_penalty, genet_label, fill = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.035) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.045) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_col(width = 0.55, colour = "black", linewidth = 0.25,
           alpha = 0.88) +
  geom_text(aes(x = label_x,
                label = sprintf("%.2f", mean_penalty),
                hjust = label_hjust),
            size = 3, colour = "grey20") +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2),
                     labels = c("earlier", "no\ndelay", "delayed")) +
  coord_cartesian(xlim = c(-0.3, 0.3), clip = "off") +
  labs(x = "Mean heat effect",
       y = NULL,
       title = "A. Average across milestones") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        legend.position = "none")

p_morph_traits <- ggplot(morph_penalty,
                         aes(z, trait_label,
                             colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.025) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.035) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_point(size = 3.1, alpha = 0.95,
             position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = PAL_GENO, name = "Genet",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Genet",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("earlier\nat 31C", "",
                                "no\ndelay", "",
                                "delayed\nat 31C")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Heat effect on timing within each milestone",
       y = NULL,
       title = "B. Which steps were delayed?") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2))

p_morph <- p_morph_mean + p_morph_traits +
  plot_layout(widths = c(0.8, 2.2), guides = "collect") +
  plot_annotation(
    title = "Heat delays wound healing less in genet C",
    subtitle = str_wrap(
      "Dots left of the dashed line reached a milestone earlier at 31C; dots right of the line were delayed. C clusters near no delay or earlier onset, while A and D show the largest delays for some steps.",
      width = 105
    ),
    caption = str_wrap(
      "Scores compare genets within each milestone after putting Cox timing effects on a common scale. n = 8 wounded fragments per genet and milestone (4 at 28C, 4 at 31C); pigment over wound for genet A is omitted because the hazard ratio was non-finite.",
      width = 105
    )
  ) &
  theme_pub(9) &
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

save_fig(p_morph, "19d_wound_healing_heat_penalties", width = 185,
         height = 115)

# Per-genet × scope mean sensitivity
resilience_decomp <- decomp |>
  group_by(thicket, scope) |>
  summarise(mean_sensitivity = mean(z, na.rm = TRUE),
            n_responses      = n(),
            .groups = "drop")

write_csv(resilience_decomp,
          file.path(TBL_DIR, "19c_resilience_decomp_by_scope.csv"))

cat("\n=== Genet resilience summary ===\n")
print(resilience |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\n=== Decomposed resilience by scope ===\n")
print(resilience_decomp |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\nWrote 19_genet_dashboard.{pdf,png}, 19b_genet_resilience_ranking.{pdf,png},",
    "19c_decomposed_resilience.{pdf,png},",
    "19d_wound_healing_heat_penalties.{pdf,png},",
    "19_genet_resilience_summary.csv,",
    "19c_resilience_decomp_by_scope.csv\n")
