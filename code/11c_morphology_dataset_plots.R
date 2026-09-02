# =============================================================================
# Purpose: Comprehensive comparison plots for the two morphology datasets.
#
# What & why: 04_physio_morphology.R and 11_microscope_physio.R each make
#   dataset-specific trajectory plots. This script adds the cross-dataset view:
#   it keeps the main gross morphology cohort and the microscope/photo cohort
#   separate, but puts their observed trait trajectories, endpoint states, and
#   first-observed timing summaries into one reproducible plotting branch.
#
# Input:   data/processed/physio_clean.rds
#          data/processed/microscope_physio_clean.rds
# Output:  output/tables/11c_morphology_dataset_plot_design.csv
#          output/tables/11c_morphology_dataset_trait_trajectory.csv
#          output/tables/11e_morphology_dataset_event_summary.csv
#          output/tables/11e_morphology_dataset_treatment_deltas.csv
#          output/tables/11f_morphology_dataset_endpoint_tests.csv
#          output/tables/11f_morphology_dataset_logrank_tests.csv
#          output/tables/11f_morphology_dataset_analysis_manifest.csv
#          output/diagnostics/11f_morphology_dataset_diagnostic_checks.csv
#          output/diagnostics/11f_morphology_dataset_analysis_report.md
#          figures/11c_morphology_dataset_all_trait_trajectories.{pdf,png}
#          figures/11d_morphology_dataset_shared_trait_comparison.{pdf,png}
#          figures/11e_morphology_dataset_event_summary.{pdf,png}
#          figures/diagnostics/11f_morphology_dataset_scoring_coverage.{pdf,png}
#          figures/diagnostics/11f_morphology_dataset_endpoint_balance.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(survival))

DIAG_OUT <- file.path(OUT_DIR, "diagnostics")
DIAG_FIG <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_FIG, recursive = TRUE, showWarnings = FALSE)

physio_path <- file.path(DATA_PROC, "physio_clean.rds")
micro_path <- file.path(DATA_PROC, "microscope_physio_clean.rds")
missing_inputs <- c(physio_path, micro_path)[!file.exists(c(physio_path, micro_path))]
if (length(missing_inputs)) {
  stop("Missing cleaned morphology input(s): ",
       paste(missing_inputs, collapse = ", "), call. = FALSE)
}

ph <- readRDS(physio_path)
micro <- readRDS(micro_path)

trait_meta <- tibble::tribble(
  ~trait,                  ~short_label,       ~label,                    ~phase,              ~main, ~micro,
  "tissue_over_wound",     "Tissue cover",     "Tissue over wound",       "Wound covering",    FALSE, TRUE,
  "polyps_out",            "Polyps out",       "Polyps out",             "General condition", TRUE,  FALSE,
  "axial_polyp_formation", "Axial polyp",      "Axial polyp formation",  "Wound closure",     TRUE,  TRUE,
  "wound_smoothed",        "Wound smoothed",   "Wound smoothed",         "Wound closure",     TRUE,  TRUE,
  "pigment_over_wound",    "Pigment",          "Pigment over wound",     "Wound covering",    TRUE,  TRUE,
  "tip_exist",             "Tip exists",       "Tip exists",             "Regeneration",      TRUE,  TRUE,
  "tip_extension",         "Tip extension",    "Tip extension",          "Regeneration",      TRUE,  TRUE,
  "new_corallites_on_tip", "New corallites",   "New corallites on tip",  "Regeneration",      TRUE,  TRUE,
  "algae_on_wound",        "Algae",            "Algae on wound",         "Fouling",           TRUE,  FALSE
) |>
  mutate(
    trait_order = row_number(),
    trait = factor(trait, levels = trait),
    short_label = factor(short_label, levels = short_label),
    phase = factor(phase, levels = c("General condition", "Wound covering",
                                     "Wound closure", "Regeneration", "Fouling")),
    phase_short = dplyr::recode(
      as.character(phase),
      `General condition` = "Condition",
      `Wound covering` = "Covering",
      `Wound closure` = "Closure",
      Regeneration = "Regeneration",
      Fouling = "Fouling"
    )
  )

prep_long <- function(dat, dataset, dataset_label, trait_cols) {
  dat |>
    filter(as.character(wound) == "yes") |>
    select(id, date, day, treatment, tank, thicket, all_of(trait_cols)) |>
    pivot_longer(all_of(trait_cols), names_to = "trait",
                 values_to = "expressed") |>
    mutate(
      dataset = dataset,
      dataset_label = dataset_label,
      dataset_short = if_else(dataset == "physio_morphology", "Main", "Microscope"),
      dataset_id = paste(dataset, id, sep = "_"),
      treatment = factor(as.character(treatment), levels = c("28C", "31C")),
      expressed = as.integer(expressed)
    ) |>
    left_join(trait_meta, by = "trait") |>
    filter(!is.na(expressed))
}

main_traits <- trait_meta |> filter(main) |> pull(trait) |> as.character()
micro_traits <- trait_meta |> filter(micro) |> pull(trait) |> as.character()

long <- bind_rows(
  prep_long(ph, "physio_morphology", "Main morphology", main_traits),
  prep_long(micro, "microscope_physio", "Microscope photos", micro_traits)
) |>
  mutate(
    dataset_label = factor(dataset_label,
                           levels = c("Main morphology", "Microscope photos")),
    facet_trait = paste(dataset_short, short_label, sep = ": "),
    facet_event = paste(dataset_short, phase_short, sep = "\n")
  )

dataset_design <- long |>
  distinct(dataset, dataset_label, dataset_id, id, treatment, thicket, tank) |>
  count(dataset, dataset_label, treatment, thicket, tank,
        name = "n_corals") |>
  arrange(dataset, treatment, thicket, tank)
write_csv(dataset_design,
          file.path(TBL_DIR, "11c_morphology_dataset_plot_design.csv"))

trajectory <- long |>
  group_by(dataset, dataset_label, day, treatment, trait, short_label, label,
           trait_order, phase, phase_short, facet_trait) |>
  summarise(
    n_scored = n(),
    n_yes = sum(expressed == 1),
    prop = n_yes / n_scored,
    .groups = "drop"
  ) |>
  arrange(dataset, trait_order, treatment, day)
write_csv(trajectory,
          file.path(TBL_DIR, "11c_morphology_dataset_trait_trajectory.csv"))

n_text <- long |>
  distinct(dataset_label, treatment, dataset_id) |>
  count(dataset_label, treatment, name = "n_corals") |>
  mutate(part = paste0(dataset_label, " ", treatment, " n=", n_corals)) |>
  pull(part) |>
  paste(collapse = "; ")

facet_trait_levels <- trajectory |>
  distinct(dataset_label, trait_order, short_label, facet_trait) |>
  arrange(dataset_label, trait_order) |>
  pull(facet_trait)

trajectory_all <- trajectory |>
  mutate(facet_trait = factor(facet_trait, levels = facet_trait_levels))

p_all_traits <- ggplot(trajectory_all,
                       aes(day, prop, colour = treatment,
                           shape = treatment, group = treatment)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.25) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_wrap(~ facet_trait, ncol = 4) +
  labs(
    x = "Days post-wounding",
    y = "Corals expressing trait",
    title = "All scored wound-healing traits across both morphology datasets",
    subtitle = paste("Wounded corals only; cohorts are plotted separately and not pooled.",
                     n_text),
    caption = "Microscope pigment was scored through D7 only; blank D8-D15 pigment cells are omitted."
  ) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(4, "pt"),
    strip.text = element_text(size = 7.1),
    legend.position = "bottom"
  )
save_fig(p_all_traits, "11c_morphology_dataset_all_trait_trajectories",
         width = 220, height = 190)

shared_traits <- trait_meta |> filter(main, micro) |> pull(trait) |> as.character()
shared_levels <- trait_meta |>
  filter(as.character(trait) %in% shared_traits) |>
  pull(short_label)

trajectory_shared <- trajectory |>
  filter(as.character(trait) %in% shared_traits) |>
  mutate(short_label = factor(short_label, levels = shared_levels))

p_shared <- ggplot(trajectory_shared,
                   aes(day, prop, colour = treatment,
                       shape = treatment, group = treatment)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_grid(dataset_label ~ short_label) +
  labs(
    x = "Days post-wounding",
    y = "Corals expressing trait",
    title = "Shared morphology traits show the same closure-versus-regeneration split",
    subtitle = "Main morphology and microscope/photo cohorts are separate datasets; shared axes make the visual comparison explicit.",
    caption = "Points are observed treatment-level proportions. Microscope pigment was scored through D7 only."
  ) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(3.5, "pt"),
    strip.text = element_text(size = 7.1),
    legend.position = "bottom"
  )
save_fig(p_shared, "11d_morphology_dataset_shared_trait_comparison",
         width = 245, height = 125)

# ---- Event and endpoint summaries -----------------------------------------
events_by_coral <- long |>
  arrange(dataset, trait, treatment, dataset_id, day) |>
  group_by(dataset, dataset_label, dataset_id, id, treatment, thicket, tank,
           trait, short_label, label, trait_order, phase, phase_short,
           facet_event) |>
  summarise(
    n_scored_days = n(),
    first_scored_day = min(day),
    last_scored_day = max(day),
    event = any(expressed == 1),
    first_yes_day = if (any(expressed == 1)) min(day[expressed == 1]) else NA_real_,
    endpoint_day = max(day),
    endpoint_yes = expressed[which.max(day)],
    .groups = "drop"
  )

event_summary <- events_by_coral |>
  group_by(dataset, dataset_label, trait, short_label, label, trait_order,
           phase, phase_short, facet_event, treatment) |>
  summarise(
    n_corals = n(),
    n_events = sum(event),
    prop_ever = n_events / n_corals,
    median_first_yes_day = if (n_events > 0) median(first_yes_day[event]) else NA_real_,
    first_yes_q25 = if (n_events > 0) unname(quantile(first_yes_day[event], 0.25)) else NA_real_,
    first_yes_q75 = if (n_events > 0) unname(quantile(first_yes_day[event], 0.75)) else NA_real_,
    endpoint_day = max(endpoint_day),
    n_endpoint_scored = sum(!is.na(endpoint_yes)),
    n_endpoint_yes = sum(endpoint_yes == 1, na.rm = TRUE),
    prop_endpoint = n_endpoint_yes / n_endpoint_scored,
    .groups = "drop"
  ) |>
  arrange(dataset, trait_order, treatment)
write_csv(event_summary,
          file.path(TBL_DIR, "11e_morphology_dataset_event_summary.csv"))

treatment_deltas <- event_summary |>
  select(dataset, dataset_label, trait, short_label, label, phase, facet_event,
         trait_order, phase_short, treatment, n_corals, prop_ever,
         median_first_yes_day, endpoint_day, n_endpoint_scored,
         n_endpoint_yes, prop_endpoint) |>
  pivot_wider(
    names_from = treatment,
    values_from = c(n_corals, prop_ever, median_first_yes_day, endpoint_day,
                    n_endpoint_scored, n_endpoint_yes, prop_endpoint),
    names_sep = "_"
  ) |>
  mutate(
    endpoint_delta_31C_minus_28C = prop_endpoint_31C - prop_endpoint_28C,
    ever_delta_31C_minus_28C = prop_ever_31C - prop_ever_28C,
    median_first_day_delta_31C_minus_28C =
      median_first_yes_day_31C - median_first_yes_day_28C
  ) |>
  arrange(dataset, trait_order)
write_csv(treatment_deltas,
          file.path(TBL_DIR, "11e_morphology_dataset_treatment_deltas.csv"))

event_levels <- event_summary |>
  distinct(dataset_label, phase, facet_event) |>
  arrange(dataset_label, phase) |>
  pull(facet_event)

event_plot_data <- event_summary |>
  mutate(
    short_label = forcats::fct_rev(short_label),
    facet_event = factor(facet_event, levels = event_levels)
  )

endpoint_segments <- treatment_deltas |>
  mutate(
    short_label = forcats::fct_rev(short_label),
    facet_event = factor(facet_event, levels = event_levels)
  ) |>
  filter(!is.na(prop_endpoint_28C), !is.na(prop_endpoint_31C))

timing_points <- event_plot_data |>
  filter(!is.na(median_first_yes_day))

timing_segments <- treatment_deltas |>
  mutate(
    short_label = forcats::fct_rev(short_label),
    facet_event = factor(facet_event, levels = event_levels)
  ) |>
  filter(!is.na(median_first_yes_day_28C),
         !is.na(median_first_yes_day_31C))

p_endpoint <- ggplot() +
  geom_segment(
    data = endpoint_segments,
    aes(x = prop_endpoint_28C, xend = prop_endpoint_31C,
        y = short_label, yend = short_label),
    colour = "grey72", linewidth = 0.55
  ) +
  geom_point(
    data = event_plot_data,
    aes(x = prop_endpoint, y = short_label, colour = treatment,
        shape = treatment),
    size = 1.9
  ) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_y_discrete(drop = TRUE) +
  facet_wrap(~ facet_event, ncol = 4, scales = "free_y") +
  labs(x = "Corals expressing trait at last scored day", y = NULL) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(4, "pt"),
    strip.text = element_text(size = 7.1, lineheight = 0.95),
    legend.position = "bottom"
  )

p_timing <- ggplot() +
  geom_segment(
    data = timing_segments,
    aes(x = median_first_yes_day_28C, xend = median_first_yes_day_31C,
        y = short_label, yend = short_label),
    colour = "grey72", linewidth = 0.55
  ) +
  geom_point(
    data = timing_points,
    aes(x = median_first_yes_day, y = short_label, colour = treatment,
        shape = treatment),
    size = 1.9
  ) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_discrete(drop = TRUE) +
  facet_wrap(~ facet_event, ncol = 4, scales = "free_y") +
  labs(x = "Median first day observed, among corals with the trait", y = NULL) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(4, "pt"),
    strip.text = element_text(size = 7.1, lineheight = 0.95),
    legend.position = "bottom"
  )

p_events <- (p_endpoint / p_timing) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Endpoint state and first-observed timing across both morphology datasets",
    subtitle = "Every scored trait is summarized within its own cohort; grey connectors link 28C and 31C for the same trait.",
    caption = "Panel A: state at the last scored day for that trait. Panel B: first-observed timing only among corals where the trait occurred."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold", hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, hjust = 0,
                                 margin = margin(b = 6), colour = "grey30"),
    plot.caption = element_text(size = 7.5, hjust = 0, colour = "grey40")
  )

save_fig(p_events, "11e_morphology_dataset_event_summary",
         width = 245, height = 190)

# ---- Statistical tests paired with the cross-dataset plots -----------------
# These are intentionally simple, plot-matched tests: endpoint Fisher exact
# tests and first-observed-day log-rank tests within each dataset. They do not
# replace the main-cohort primary survival models in script 14, which use
# interval censoring and thicket adjustment. Their job is to provide a consistent
# table next to the cross-dataset descriptive figures.
endpoint_tests <- events_by_coral |>
  group_by(dataset, dataset_label, trait, short_label, label, trait_order,
           phase, phase_short) |>
  group_modify(\(.x, .y) {
    d <- .x |> filter(!is.na(endpoint_yes),
                      treatment %in% c("28C", "31C"))
    n_28 <- sum(d$treatment == "28C")
    n_31 <- sum(d$treatment == "31C")
    yes_28 <- sum(d$endpoint_yes[d$treatment == "28C"] == 1)
    yes_31 <- sum(d$endpoint_yes[d$treatment == "31C"] == 1)
    no_28 <- n_28 - yes_28
    no_31 <- n_31 - yes_31
    endpoint_day <- if (nrow(d) > 0) max(d$endpoint_day, na.rm = TRUE) else NA_real_
    if (n_28 == 0 || n_31 == 0) {
      return(tibble(
        endpoint_day = endpoint_day,
        n_28 = n_28, yes_28 = yes_28, prop_28 = NA_real_,
        n_31 = n_31, yes_31 = yes_31, prop_31 = NA_real_,
        risk_diff_31_minus_28 = NA_real_,
        odds_ratio_31_vs_28 = NA_real_,
        fisher_p = NA_real_,
        test_status = "not_estimable"
      ))
    }
    mat <- matrix(c(yes_28, no_28, yes_31, no_31),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(treatment = c("28C", "31C"),
                                  endpoint = c("yes", "no")))
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    odds_31 <- (yes_31 + 0.5) / (no_31 + 0.5)
    odds_28 <- (yes_28 + 0.5) / (no_28 + 0.5)
    tibble(
      endpoint_day = endpoint_day,
      n_28 = n_28,
      yes_28 = yes_28,
      prop_28 = yes_28 / n_28,
      n_31 = n_31,
      yes_31 = yes_31,
      prop_31 = yes_31 / n_31,
      risk_diff_31_minus_28 = yes_31 / n_31 - yes_28 / n_28,
      odds_ratio_31_vs_28 = odds_31 / odds_28,
      fisher_p = if (is.null(ft)) NA_real_ else ft$p.value,
      test_status = if (is.null(ft)) "not_estimable" else "computed"
    )
  }) |>
  ungroup() |>
  mutate(
    test = "Fisher exact test at last scored day",
    test_scope = "paired with figures 11c-11e; cohorts tested separately"
  ) |>
  arrange(dataset, trait_order)
write_csv(endpoint_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"))

logrank_tests <- events_by_coral |>
  mutate(event_time = if_else(event, first_yes_day, last_scored_day)) |>
  group_by(dataset, dataset_label, trait, short_label, label, trait_order,
           phase, phase_short) |>
  group_modify(\(.x, .y) {
    d <- .x |>
      filter(!is.na(event_time), treatment %in% c("28C", "31C")) |>
      mutate(treatment = factor(as.character(treatment),
                                levels = c("28C", "31C")),
             event_int = as.integer(event))
    n_28 <- sum(d$treatment == "28C")
    n_31 <- sum(d$treatment == "31C")
    n_event_28 <- sum(d$event_int[d$treatment == "28C"])
    n_event_31 <- sum(d$event_int[d$treatment == "31C"])
    if (n_28 == 0 || n_31 == 0 || sum(d$event_int) == 0) {
      return(tibble(
        n_28 = n_28, n_31 = n_31,
        n_event_28 = n_event_28, n_event_31 = n_event_31,
        logrank_chisq = NA_real_, logrank_df = NA_real_,
        logrank_p = NA_real_,
        test_status = "not_estimable"
      ))
    }
    fit <- tryCatch(
      survdiff(Surv(event_time, event_int) ~ treatment, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble(
        n_28 = n_28, n_31 = n_31,
        n_event_28 = n_event_28, n_event_31 = n_event_31,
        logrank_chisq = NA_real_, logrank_df = NA_real_,
        logrank_p = NA_real_,
        test_status = "not_estimable"
      ))
    }
    df <- length(fit$n) - 1
    tibble(
      n_28 = n_28,
      n_31 = n_31,
      n_event_28 = n_event_28,
      n_event_31 = n_event_31,
      logrank_chisq = unname(fit$chisq),
      logrank_df = df,
      logrank_p = pchisq(fit$chisq, df = df, lower.tail = FALSE),
      test_status = "computed"
    )
  }) |>
  ungroup() |>
  mutate(
    test = "Log-rank test on first-observed day",
    test_scope = paste(
      "paired with figures 11c-11e; cohorts tested separately;",
      "main-cohort primary inference remains script 14 interval-censored AFT"
    )
  ) |>
  arrange(dataset, trait_order)
write_csv(logrank_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"))

# ---- Diagnostic/QC plots for the plot-matched analyses ---------------------
coverage_grid <- long |>
  count(dataset, dataset_label, treatment, trait, short_label, label,
        trait_order, phase, phase_short, day, name = "n_scored") |>
  complete(
    nesting(dataset, dataset_label, treatment, trait, short_label, label,
            trait_order, phase, phase_short),
    day = 0:15,
    fill = list(n_scored = 0)
  ) |>
  mutate(
    dataset_label = factor(dataset_label,
                           levels = c("Main morphology", "Microscope photos")),
    short_label = forcats::fct_reorder(short_label, trait_order)
  )

coverage_checks <- coverage_grid |>
  group_by(dataset, dataset_label, treatment, trait, short_label, label,
           trait_order) |>
  summarise(
    check = "scoring_coverage",
    n_missing_days = sum(n_scored == 0),
    missing_days = paste(day[n_scored == 0], collapse = ","),
    expected_n_scored = max(n_scored),
    n_low_count_days = sum(n_scored > 0 & n_scored < max(n_scored)),
    low_count_days = paste(day[n_scored > 0 & n_scored < max(n_scored)],
                           collapse = ","),
    all_missing_after_d7 = all(day[n_scored == 0] %in% 8:15),
    all_d0_d7_scored = all(n_scored[day %in% 0:7] > 0),
    min_n_scored = min(n_scored),
    max_n_scored = max(n_scored),
    .groups = "drop"
  ) |>
  mutate(
    missing_days = if_else(missing_days == "", "none", missing_days),
    status = case_when(
      dataset == "microscope_physio" &
        as.character(trait) == "pigment_over_wound" &
        all_missing_after_d7 & all_d0_d7_scored ~ "HANDLED",
      dataset == "physio_morphology" &
        as.character(trait) == "pigment_over_wound" &
        (missing_days == "9" | low_count_days == "9") ~ "HANDLED",
      dataset == "microscope_physio" &
        as.character(trait) == "tissue_over_wound" &
        as.character(treatment) == "31C" &
        low_count_days == "1" ~ "HANDLED",
      n_missing_days == 0 & n_low_count_days == 0 ~ "PASS",
      TRUE ~ "WARN"
    ),
    notes = case_when(
      dataset == "physio_morphology" &
        as.character(trait) == "pigment_over_wound" &
        (missing_days == "9" | low_count_days == "9") ~ "source sheet has sparse/blank main-morphology pigment scoring at D9; adjacent D8/D10 and endpoint values are scored",
      dataset == "microscope_physio" &
        as.character(trait) == "tissue_over_wound" &
        as.character(treatment) == "31C" &
        low_count_days == "1" ~ "source sheet has two heated microscope tissue-cover blanks at D1; adjacent D0/D2 values are scored",
      status == "HANDLED" ~ "expected microscope pigment blanks after D7",
      status == "PASS" ~ "all D0-D15 cells represented",
      TRUE ~ "unexpected missing scored day(s)"
    )
  )

endpoint_diagnostic <- endpoint_tests |>
  transmute(
    dataset, dataset_label, treatment = NA_character_, trait, short_label,
    label, trait_order,
    check = "endpoint_test_estimability",
    n_missing_days = NA_integer_,
    missing_days = NA_character_,
    min_n_scored = pmin(n_28, n_31),
    max_n_scored = pmax(n_28, n_31),
    status = if_else(test_status == "computed", "PASS", "FAIL"),
    notes = if_else(
      test_status == "computed",
      "endpoint Fisher exact test computed",
      "endpoint test could not be estimated"
    )
  )

logrank_diagnostic <- logrank_tests |>
  transmute(
    dataset, dataset_label, treatment = NA_character_, trait, short_label,
    label, trait_order,
    check = "first_observed_timing_estimability",
    n_missing_days = NA_integer_,
    missing_days = NA_character_,
    min_n_scored = pmin(n_event_28, n_event_31),
    max_n_scored = pmax(n_event_28, n_event_31),
    status = case_when(
      test_status == "computed" ~ "PASS",
      as.character(phase) == "General condition" ~ "HANDLED",
      n_event_28 == 0 & n_event_31 == 0 ~ "HANDLED",
      TRUE ~ "WARN"
    ),
    notes = case_when(
      test_status == "computed" ~ "log-rank test computed",
      as.character(phase) == "General condition" ~ "baseline condition trait; first-observed timing test is not biologically meaningful",
      n_event_28 == 0 & n_event_31 == 0 ~ "no events in either treatment",
      TRUE ~ "timing test could not be estimated"
    )
  )

diagnostic_checks <- bind_rows(coverage_checks, endpoint_diagnostic,
                               logrank_diagnostic) |>
  arrange(dataset, trait_order, check, treatment)
write_csv(diagnostic_checks,
          file.path(DIAG_OUT, "11f_morphology_dataset_diagnostic_checks.csv"))

p_coverage <- ggplot(coverage_grid,
                     aes(day, forcats::fct_rev(short_label),
                         fill = n_scored)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "grey96", high = PAL_OKABE[[5]],
                      name = "Scored\ncorals") +
  scale_x_continuous(breaks = seq(0, 15, by = 3),
                     expand = expansion(mult = 0)) +
  facet_grid(dataset_label ~ treatment, scales = "free_y",
             space = "free_y") +
  labs(
    x = "Days post-wounding",
    y = NULL,
    title = "Scoring coverage for morphology plot-matched analyses",
    subtitle = "Tiles show the number of corals contributing to each plotted day, trait, cohort, and temperature.",
    caption = "The expected microscope pigment gap after D7 is treated as handled, not as a zero response."
  ) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(5, "pt"),
    strip.text = element_text(size = 7.5),
    legend.position = "right"
  )
save_fig(p_coverage, "11f_morphology_dataset_scoring_coverage",
         width = 180, height = 160, dir = DIAG_FIG)

endpoint_balance <- events_by_coral |>
  mutate(
    endpoint_state = if_else(endpoint_yes == 1, "Yes", "No"),
    dataset_label = factor(dataset_label,
                           levels = c("Main morphology", "Microscope photos")),
    short_label = forcats::fct_reorder(short_label, trait_order)
  ) |>
  count(dataset_label, treatment, trait, short_label, trait_order,
        endpoint_state, name = "n")

p_endpoint_balance <- ggplot(endpoint_balance,
                             aes(n, forcats::fct_rev(short_label),
                                 fill = endpoint_state)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(No = "grey78", Yes = "black"),
                    name = "Endpoint") +
  scale_x_continuous(breaks = scales::breaks_width(4),
                     expand = expansion(mult = c(0, 0.06))) +
  facet_grid(dataset_label ~ treatment, scales = "free_y",
             space = "free_y") +
  labs(
    x = "Corals at last scored day",
    y = NULL,
    title = "Endpoint balance for morphology endpoint tests",
    subtitle = "Stacked bars show yes/no counts entering each Fisher exact endpoint test.",
    caption = "Saturated endpoints are still valid for Fisher exact tests, but they are a warning against endpoint logistic models."
  ) +
  theme_pub(8) +
  theme(
    panel.spacing = unit(5, "pt"),
    strip.text = element_text(size = 7.5),
    legend.position = "bottom"
  )
save_fig(p_endpoint_balance, "11f_morphology_dataset_endpoint_balance",
         width = 180, height = 160, dir = DIAG_FIG)

analysis_manifest <- tibble::tribble(
  ~display_item, ~paired_analysis, ~diagnostic_artifact, ~scope, ~caveat,
  "figures/11c_morphology_dataset_all_trait_trajectories.pdf",
  "output/tables/11c_morphology_dataset_trait_trajectory.csv",
  "figures/diagnostics/11f_morphology_dataset_scoring_coverage.pdf",
  "all scored traits in both morphology cohorts",
  "descriptive trajectory figure; inferential tests are in 11f endpoint/log-rank tables",
  "figures/11d_morphology_dataset_shared_trait_comparison.pdf",
  "output/tables/11f_morphology_dataset_endpoint_tests.csv; output/tables/11f_morphology_dataset_logrank_tests.csv",
  "figures/diagnostics/11f_morphology_dataset_scoring_coverage.pdf",
  "shared traits only, cohorts kept separate",
  "plot-matched unadjusted tests; main morphology primary inference remains script 14",
  "figures/11e_morphology_dataset_event_summary.pdf",
  "output/tables/11e_morphology_dataset_event_summary.csv; output/tables/11e_morphology_dataset_treatment_deltas.csv; output/tables/11f_morphology_dataset_endpoint_tests.csv; output/tables/11f_morphology_dataset_logrank_tests.csv",
  "figures/diagnostics/11f_morphology_dataset_endpoint_balance.pdf",
  "endpoint state and first-observed timing for every scored trait",
  "Fisher exact handles saturated endpoint counts; log-rank is skipped/handled when no trait events occur",
  "figures/04_morphology_trajectories.pdf; figures/14_morphology_KM.pdf",
  "output/tables/04_morphology_trait_glmm_summaries.csv; output/tables/14_interval_survreg.csv; output/tables/14_cox_hazard_ratios.csv",
  "figures/diagnostics/B_<trait>.png; figures/diagnostics/14_cox_ph_<trait>.png; output/tables/25_model_diagnostic_coverage.csv",
  "formal main-cohort morphology models",
  "primary regeneration timing inference is interval-censored AFT from script 14",
  "figures/11_microscope_trait_trajectories.pdf",
  "output/tables/11_microscope_event_tests.csv; output/tables/11_microscope_event_summary.csv",
  "figures/diagnostics/11f_morphology_dataset_scoring_coverage.pdf; figures/diagnostics/11f_morphology_dataset_endpoint_balance.pdf",
  "microscope/photo validation cohort",
  "small wounded-only photo cohort; not pooled with main morphology"
)
write_csv(analysis_manifest,
          file.path(TBL_DIR, "11f_morphology_dataset_analysis_manifest.csv"))

diag_counts <- diagnostic_checks |>
  count(check, status, name = "n") |>
  arrange(check, status)
report <- c(
  "# Cross-dataset morphology analyses and diagnostics",
  "",
  sprintf("Generated: %s", format(Sys.time())),
  "",
  "## Analysis tables",
  "",
  "- `output/tables/11f_morphology_dataset_endpoint_tests.csv`: Fisher exact endpoint tests for every scored trait within each cohort.",
  "- `output/tables/11f_morphology_dataset_logrank_tests.csv`: first-observed-day log-rank tests for every scored trait within each cohort.",
  "- `output/tables/11f_morphology_dataset_analysis_manifest.csv`: display item to analysis/diagnostic crosswalk.",
  "",
  "## Diagnostic plots",
  "",
  "- `figures/diagnostics/11f_morphology_dataset_scoring_coverage.pdf`: scoring coverage by cohort, temperature, trait, and day.",
  "- `figures/diagnostics/11f_morphology_dataset_endpoint_balance.pdf`: endpoint yes/no balance entering Fisher tests.",
  "",
  "## Diagnostic check counts",
  "",
  "| Check | Status | n |",
  "|---|---|---|",
  sprintf("| %s | %s | %d |", diag_counts$check, diag_counts$status, diag_counts$n),
  "",
  "## Interpretation",
  "",
  "These plot-matched analyses keep the two morphology datasets separate. For the main morphology cohort, the formal primary inference remains the interval-censored survival analysis in `code/14_morphology_kaplan.R`; the 11f tables provide the same simple endpoint/timing summaries used for the microscope cohort so the cross-dataset figures have directly paired statistics."
)
writeLines(report,
           file.path(DIAG_OUT, "11f_morphology_dataset_analysis_report.md"))

cat("\n=== Cross-dataset morphology plotting branch ===\n")
cat("Design summary:", n_text, "\n")
cat("Wrote 11c all-trait trajectories, 11d shared-trait comparison, ",
    "11e endpoint/timing summary, 11f paired tests, diagnostic plots, ",
    "and companion CSV tables.\n", sep = "")
