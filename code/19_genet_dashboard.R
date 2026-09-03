# =============================================================================
# Purpose: Integrative cross-response source-thicket dashboard. For each of the
#          3 source thickets (a, c, d), compute a standardized "heat sensitivity" effect
#          across each usable response variable (PAM, color, growth, symbionts,
#          and morphology time-to-onset traits). Build a forest plot of all
#          effect sizes to identify which source thicket is most/least thermally
#          sensitive across every dimension.
#
#          Also produces a composite heat-sensitivity ranking per source thicket: the
#          mean standardized heat effect across finite response-level effects.
#
# What & why: synthesis script. Other analyses test whether heat affects a
#   given response; this one addresses the cross-cutting question of which
#   source thicket best tolerates heat overall. To put effects measured in different
#   units on one axis (Fv/Fm, a colour score, growth,
#   symbiont counts, and wound-healing hazard ratios), each effect is rescaled
#   WITHIN its own response so the largest-magnitude effect = 1 (the "z"
#   column — a row-max standardization, not a statistical z-score). Positive
#   means the source thicket's phenotype is worse at 31 °C than 28 °C, i.e. more heat-
#   sensitive. Averaging these standardized values per source thicket yields a single
#   resilience ranking; in this dataset source C is the least heat-sensitive group. Nothing
#   here fits a new model — every number is a re-summary of upstream tables.
# Input:   output/tables/12_genet_treatment_effects.csv   (continuous responses)
#          output/tables/14_cox_hazard_ratios.csv         (per-source HR for KM)
#          output/tables/15_genet_pca_displacement.csv    (multivariate)
# Output:  figures/19_genet_dashboard.{pdf,png}            — forest plot
#          figures/19b_genet_resilience_ranking.{pdf,png}  — composite score
#          figures/19e_source_heat_penalty_summary.{pdf,png}
#                                                               — read-first
#                                                                 domain summary
#          figures/19f_source_physiology_heat_penalties.{pdf,png}
#                                                               — physiology-only
#                                                                 detail
#          figures/19d_wound_healing_heat_penalties.{pdf,png}
#                                                               — morphology-only
#                                                                 source penalty
#          output/tables/19_genet_resilience_summary.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (TBL_DIR, FIG_DIR), theme_pub(),
# save_fig(), and the source-thicket colour palette PAL_GENO used below.
source(here::here("code", "00_setup.R"))

# ---- Load upstream effect tables ------------------------------------------
# Three sources, one per analysis domain: per-source continuous LMM contrasts
# (script 12), wound-healing Cox hazard ratios (script 14), and the
# multivariate per-source PCA displacement (script 15).
cont <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                 show_col_types = FALSE)
cox  <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                 show_col_types = FALSE)
pca  <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"),
                 show_col_types = FALSE)

# ---- Standardize continuous-response effects ------------------------------
# Continuous responses (from script 12) give estimate = mean(28C) - mean(31C)
# at the response-specific endpoint per source thicket x wound, with SE:
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

# ---- Standardize Cox hazard ratios (per source thicket) -------------------
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
      response == "pam_fvfm"          ~ "Photosynthesis score\nDay 14",
      response == "color_dscale"      ~ "Color score\nDay 14",
      response == "growth_pct"        ~ "Growth\nstart to end",
      response == "log_zoox_density"  ~ "Symbiont density\nfinal tissue sample",
      grepl("^morph_", response)      ~ paste0(str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))), "\ntime-to-step"),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_pct","log_zoox_density")
        ~ "Condition + growth",
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

# ---- Simple read-first summaries -----------------------------------------
# These are intentionally aggregated views for the team-summary HTML. The
# detailed dashboard below keeps the response-level audit trail.
source_order <- c("A", "D", "C")

domain_summary <- all_eff |>
  mutate(summary_group = case_when(
    domain == "Condition + growth" ~ "Whole-coral condition + growth",
    domain == "Wound closure" ~ "Wound closure",
    domain == "Regeneration" ~ "Regeneration",
    TRUE                     ~ domain
  )) |>
  group_by(summary_group, thicket) |>
  summarise(mean_penalty = mean(z, na.rm = TRUE),
            n_responses = n(),
            .groups = "drop")

overall_summary <- all_eff |>
  group_by(thicket) |>
  summarise(mean_penalty = mean(z, na.rm = TRUE),
            n_responses = n(),
            .groups = "drop") |>
  mutate(summary_group = "All responses")

quick_summary <- bind_rows(overall_summary, domain_summary) |>
  mutate(
    summary_group = factor(summary_group,
                           levels = c("All responses", "Whole-coral condition + growth",
                                      "Wound closure", "Regeneration")),
    source_label = factor(str_to_upper(thicket), levels = source_order),
    label_y = case_when(
      mean_penalty >= 0       ~ mean_penalty + 0.06,
      abs(mean_penalty) < 0.05 ~ -0.08,
      TRUE                    ~ mean_penalty / 2
    ),
    label_vjust = if_else(mean_penalty >= 0, 0, 0.5)
  )

p_quick <- ggplot(quick_summary,
                  aes(source_label, mean_penalty, fill = thicket)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#009E73", alpha = 0.03) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#D55E00", alpha = 0.035) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25,
           alpha = 0.88) +
  geom_text(aes(y = label_y,
                label = sprintf("%.2f", mean_penalty),
                vjust = label_vjust),
            size = 3.0, colour = "grey20") +
  facet_wrap(~ summary_group, nrow = 1) +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
  scale_y_continuous(
    breaks = c(-0.4, 0, 0.5, 1),
    labels = c("earlier/\nbetter", "little/no\npenalty", "0.5", "largest\npenalty")
  ) +
  coord_cartesian(ylim = c(-0.48, 1.12), clip = "off") +
  labs(x = NULL,
       y = "Average heat penalty",
       title = "Source patch C has the smallest heat penalty",
       subtitle = str_wrap(
         "Bars average 31 °C vs 28 °C heat penalties after each measurement is put on its own 0-to-1 scale. The detailed measurement-by-measurement check is shown below.",
         width = 105
       )) +
  theme_pub(9) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

save_fig(p_quick, "19e_source_heat_penalty_summary", width = 190, height = 110)

physiology_detail <- cont_eff |>
  mutate(
    response_label = case_when(
      response == "pam_fvfm"          ~ "Photosynthesis score",
      response == "color_dscale"      ~ "Color",
      response == "growth_pct"        ~ "Growth",
      response == "log_zoox_density"  ~ "Symbionts",
      TRUE                            ~ response
    ),
    response_label = factor(
      response_label,
      levels = rev(c("Photosynthesis score", "Color", "Growth", "Symbionts"))
    )
  )

p_physiology <- ggplot(physiology_detail,
                       aes(z, response_label,
                           colour = thicket, shape = thicket)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.025) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.035) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_point(size = 3.4, alpha = 0.95,
             position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Source patch",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("better\nat 31 °C", "",
                                "little/no\npenalty", "",
                                "worse\nat 31 °C")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Heat penalty within each whole-coral measurement",
       y = NULL,
       title = "Whole-coral condition: source patch C loses less under heat",
       subtitle = str_wrap(
         "Each point compares 28 °C and 31 °C for one source patch and one measurement. C sits closer to little/no penalty for all four whole-coral measurements.",
         width = 100
       )) +
  theme_pub(9) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        legend.position = "bottom")

save_fig(p_physiology, "19f_source_physiology_heat_penalties", width = 160,
         height = 82)

# Forest plot: each point is one source thicket's standardized heat effect on one
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
       title = "A. Average across measurements") +
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
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Source patch",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("better/faster\nat 31 °C", "",
                                "no\npenalty", "",
                                "worse/slower\nat 31 °C")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Heat penalty within each measurement",
       y = NULL,
       title = "B. Measurement-level heat penalties") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"))

p_dash <- p_dash_mean + p_dash_detail +
  plot_layout(widths = c(0.82, 2.35), guides = "collect") +
  plot_annotation(
    title = "Source patch C is least harmed by heat",
    subtitle = str_wrap(
      "A, C, and D are source-patch labels, not confirmed genetic individuals. Each measurement is first converted to a 31 °C vs 28 °C heat penalty and put on a common 0-to-1 scale. C stays close to little/no penalty, while A and D show larger losses or delays.",
      width = 105
    ),
    caption = str_wrap(
      "Positive values mean worse whole-coral condition or slower healing/regrowth at 31 °C; negative values mean higher condition or faster healing/regrowth at 31 °C. Whole-coral measurements use final or full-study values averaged across wound state. Healing/regrowth measurements use wounded fragments and the time when each step appeared.",
      width = 105
    )
  ) &
  theme_pub(9) &
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

save_fig(p_dash, "19_genet_dashboard", width = 200, height = 185)

# ---- Composite thermal resilience score ----------------------------------
# Per source thicket: mean standardized heat sensitivity across finite response effects.
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
# rank_overall column) back in as the "source resilience" rows.
write_csv(resilience, file.path(TBL_DIR, "19_genet_resilience_summary.csv"))

# Bar chart of the composite ranking; each bar is annotated with that source's
# multivariate PCA displacement as an independent cross-check on the ordering.
p_rank <- ggplot(resilience,
                  aes(reorder(thicket, -mean_sensitivity),
                      mean_sensitivity, fill = thicket)) +
  geom_col(width = 0.55, alpha = 0.85, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("n=%d\nPCA = %.2f",
                                n_responses, pca_displacement)),
            vjust = -0.25, lineheight = 0.95, size = 3, colour = "grey20") +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
	  labs(x = "Source patch",
	       y = "Mean heat penalty\n(0-to-1 scale)",
	       title = "Composite of standardized heat penalties",
	       subtitle = str_wrap(
	         "Lower bars mean smaller heat penalties. Inputs are 31 °C vs 28 °C comparisons or healing/regrowth timing effects, not same-day raw values.",
      width = 55
    ),
    caption = str_wrap(
      paste0("Finite effects averaged per source patch: ", composite_n_note,
             ". Missing timing values are omitted. The multivariate shift is a separate descriptive check using final whole-coral condition values put on the same scale."),
      width = 70
    )) +
  coord_cartesian(ylim = c(0, max(resilience$mean_sensitivity) * 1.35)) +
  theme_pub(10)

save_fig(p_rank, "19b_genet_resilience_ranking", width = 130, height = 110)

# ---- Decomposed dashboard: heat-only vs heat-while-wounded ---------------
# The composite above pools across wound state. To answer "is source C
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
      response == "pam_fvfm"         ~ "Photosynthesis score\nDay 14",
      response == "color_dscale"     ~ "Color score\nDay 14",
      response == "growth_pct"       ~ "Growth\nstart to end",
      response == "log_zoox_density" ~ "Symbiont density\nfinal tissue sample",
      grepl("^morph_", response)     ~ paste0(str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))), "\ntime to step"),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_pct","log_zoox_density")
        ~ "Condition + growth",
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
  scale_colour_manual(values = PAL_GENO, name = "Source patch") +
  scale_shape_manual(values = c(a = 16, c = 17, d = 15), name = "Source patch") +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1", "-0.5", "0", "0.5", "largest\npenalty")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Heat penalty within each measurement and wound-state group",
       y = NULL,
       title = "Heat penalty by wound state",
       subtitle = str_wrap(
         "Whole-coral panels compare 31 °C with 28 °C separately for wounded and unwounded fragments. Healing/regrowth panels use wounded fragments and the time when each step appeared.",
         width = 115
       ),
    caption = str_wrap(
      "Values are scaled within each measurement and wound-state group. Positive values mean a larger heat penalty; negative values mean the 31 °C group reached a healing/regrowth step faster or had a higher final value.",
      width = 115
    )) +
  theme_pub(9) +
  theme(panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold"))

save_fig(p_decomp, "19c_decomposed_resilience", width = 200, height = 175)

# ---- Focused wound-healing heat penalties --------------------------------
# Same row-max-scaled Cox effects as the dashboard, filtered to the wounded-only
# morphology milestones so the source-thicket pattern is visible without the physiology
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
  labs(x = "Mean heat penalty",
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
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
  scale_shape_manual(values = SHP_GENO, name = "Source patch",
                     labels = c(a = "A", c = "C", d = "D")) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("earlier\nat 31 °C", "",
                                "no\ndelay", "",
                                "delayed\nat 31 °C")) +
  coord_cartesian(xlim = c(-1.05, 1.05)) +
  labs(x = "Heat penalty on timing within each milestone",
       y = NULL,
       title = "B. Which steps were delayed?") +
  theme_pub(9) +
  theme(plot.title = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2))

p_morph <- p_morph_mean + p_morph_traits +
  plot_layout(widths = c(0.8, 2.2), guides = "collect") +
  plot_annotation(
    title = "Heat delayed healing and regrowth less in source patch C",
    subtitle = str_wrap(
      "Dots left of the dashed line reached a milestone earlier at 31 °C; dots right of the line were delayed. C clusters near no delay or earlier onset, while A and D show the largest delays for some steps.",
      width = 105
    ),
    caption = str_wrap(
      "Scores compare source patches within each healing/regrowth step after putting timing effects on a common scale. n = 8 wounded fragments per source patch and step (4 at 28 °C, 4 at 31 °C); pigment over wound for source A is omitted because the timing value could not be estimated.",
      width = 105
    )
  ) &
  theme_pub(9) &
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

save_fig(p_morph, "19d_wound_healing_heat_penalties", width = 185,
         height = 115)

# Per-source-thicket × scope mean sensitivity
resilience_decomp <- decomp |>
  group_by(thicket, scope) |>
  summarise(mean_sensitivity = mean(z, na.rm = TRUE),
            n_responses      = n(),
            .groups = "drop")

write_csv(resilience_decomp,
          file.path(TBL_DIR, "19c_resilience_decomp_by_scope.csv"))

cat("\n=== Source-thicket resilience summary ===\n")
print(resilience |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\n=== Decomposed resilience by scope ===\n")
print(resilience_decomp |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\nWrote 19_genet_dashboard.{pdf,png}, 19b_genet_resilience_ranking.{pdf,png},",
    "19c_decomposed_resilience.{pdf,png},",
    "19e_source_heat_penalty_summary.{pdf,png},",
    "19f_source_physiology_heat_penalties.{pdf,png},",
    "19d_wound_healing_heat_penalties.{pdf,png},",
    "19_genet_resilience_summary.csv,",
    "19c_resilience_decomp_by_scope.csv\n")
