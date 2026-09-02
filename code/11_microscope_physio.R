# =============================================================================
# Purpose: Separate analysis of the 16-coral microscope/photo wound-healing
#          cohort. This is NOT pooled with the main physio_morphology data.
#
# What & why: the main morphology analysis uses the larger physio/gene scoring
#   table. The microscope sheet is a distinct photo-only cohort with daily
#   scoring from D0-D15, wounded corals only, and only genets A/C. It is useful
#   as visual/validation evidence for the core biological contrast: tissue
#   coverage and wound smoothing happen readily, while new corallite formation
#   is reduced under heat. Because this dataset has only 16 corals and no
#   unwounded control, the analysis is deliberately descriptive plus
#   exploratory Fisher/log-rank checks.
#
# Input:   data/raw/microscope_physio/data.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/microscope_physio_clean.rds
#          output/tables/11_microscope_design_summary.csv
#          output/tables/11_microscope_design_crosswalk.csv
#          output/tables/11_microscope_trait_trajectory.csv
#          output/tables/11_microscope_trait_trajectory_by_genet.csv
#          output/tables/11_microscope_event_summary.csv
#          output/tables/11_microscope_event_tests.csv
#          output/tables/11_microscope_use_photo_series.csv
#          figures/11_microscope_trait_trajectories.{pdf,png}
#          figures/11b_microscope_trait_trajectories_by_genet.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(survival))

dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_PROC, recursive = TRUE, showWarnings = FALSE)

MICROSCOPE_SOURCE_URL <- "https://docs.google.com/spreadsheets/d/1rB5VFJqQov0ZXsPjRKnjwVDOtiPxgbfMMDm9l_9Y-Q4/edit"
MICROSCOPE_SOURCE_MODIFIED_UTC <- "2025-10-07T18:34:33.762Z"

raw_path <- file.path(DATA_RAW, "microscope_physio", "data.csv")
if (!file.exists(raw_path)) {
  stop("Missing microscope raw export: ", raw_path, call. = FALSE)
}

raw <- readr::read_csv(raw_path, show_col_types = FALSE, guess_max = 1000,
                       na = c("", "NA")) |>
  janitor::clean_names()

expected_cols <- c(
  "thicket", "wound", "tank", "sample", "treatment", "id", "date", "day",
  "tissue_over_wound", "hole_in_center", "polyp_in_hole", "wound_smoothed",
  "pigment_over_wound", "tip_exist", "tip_extension",
  "new_corallites_on_tip", "use_photo", "notes"
)
missing_cols <- setdiff(expected_cols, names(raw))
if (length(missing_cols)) {
  stop("Microscope raw export missing columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

yn_to_int <- function(x, column) {
  y <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
  y[y %in% c("", "na", "nan")] <- NA_character_
  unexpected <- setdiff(stats::na.omit(unique(y)), c("yes", "no"))
  if (length(unexpected)) {
    stop(sprintf("Unexpected yes/no code in %s: %s",
                 column, paste(unexpected, collapse = ", ")), call. = FALSE)
  }
  ifelse(is.na(y), NA_integer_, as.integer(y == "yes"))
}

trait_cols <- c(
  "tissue_over_wound", "hole_in_center", "polyp_in_hole", "wound_smoothed",
  "pigment_over_wound", "tip_exist", "tip_extension",
  "new_corallites_on_tip"
)

micro <- raw |>
  mutate(
    date = lubridate::mdy(date),
    day = as.integer(day),
    treatment = factor(paste0(as.integer(treatment), "C"),
                       levels = c("28C", "31C")),
    tank = as.integer(tank),
    thicket = stringr::str_to_lower(stringr::str_squish(thicket)),
    id = as.integer(id),
    wound = factor(stringr::str_to_lower(stringr::str_squish(wound)),
                   levels = c("no", "yes")),
    sample = stringr::str_to_lower(stringr::str_squish(sample)),
    notes = dplyr::na_if(as.character(notes), "NA"),
    across(all_of(c(trait_cols, "use_photo")),
           ~ yn_to_int(.x, dplyr::cur_column()))
  )

# The microscope sheet uses these as a single observable, as in the main
# morphology script. Stop if a future re-score makes them diverge.
if (!identical(micro$hole_in_center, micro$polyp_in_hole)) {
  stop("hole_in_center and polyp_in_hole diverge in microscope data; ",
       "do not combine without reviewing the scoring convention.", call. = FALSE)
}
micro <- micro |> mutate(axial_polyp_formation = hole_in_center)

saveRDS(micro, file.path(DATA_PROC, "microscope_physio_clean.rds"))

# ---- Design / provenance checks -------------------------------------------
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds")) |>
  filter(sample == "photo") |>
  transmute(
    id,
    meta_tank = as.integer(tank),
    meta_treatment = as.character(treatment),
    meta_thicket = as.character(thicket),
    meta_wound = as.character(wound),
    meta_sample = as.character(sample)
  )

design_crosswalk <- micro |>
  distinct(id, tank, treatment, thicket, wound, sample) |>
  mutate(treatment = as.character(treatment),
         wound = as.character(wound)) |>
  left_join(meta, by = "id") |>
  mutate(
    matches_metadata = !is.na(meta_tank) &
      tank == meta_tank &
      treatment == meta_treatment &
      thicket == meta_thicket &
      wound == meta_wound &
      sample == meta_sample
  ) |>
  arrange(treatment, thicket, tank, id)
write_csv(design_crosswalk, file.path(TBL_DIR, "11_microscope_design_crosswalk.csv"))

if (any(!design_crosswalk$matches_metadata)) {
  warning("Microscope IDs do not all match coral_metadata.rds; inspect ",
          "output/tables/11_microscope_design_crosswalk.csv", call. = FALSE)
}

design_summary <- tibble(
  dataset = "microscope_physio",
  source_title = "Microscope physio characterization",
  source_url = MICROSCOPE_SOURCE_URL,
  source_modified_utc = MICROSCOPE_SOURCE_MODIFIED_UTC,
  n_observations = nrow(micro),
  n_corals = n_distinct(micro$id),
  n_days = n_distinct(micro$day),
  day_min = min(micro$day, na.rm = TRUE),
  day_max = max(micro$day, na.rm = TRUE),
  n_treatments = n_distinct(micro$treatment),
  treatments = paste(sort(unique(as.character(micro$treatment))), collapse = ","),
  n_thickets = n_distinct(micro$thicket),
  thickets = paste(sort(unique(micro$thicket)), collapse = ","),
  n_tanks = n_distinct(micro$tank),
  tanks = paste(sort(unique(micro$tank)), collapse = ","),
  all_wounded = all(micro$wound == "yes"),
  all_photo_sample = all(micro$sample == "photo"),
  pigment_over_wound_scoring = "scored D0-D7; blank/not-scored D8-D15",
  inferential_scope = "separate photo-only validation/descriptive analysis"
)
write_csv(design_summary, file.path(TBL_DIR, "11_microscope_design_summary.csv"))

# ---- Trait summaries -------------------------------------------------------
trait_meta <- tibble::tribble(
  ~trait,                  ~label,                    ~event_interpretation,
  "tissue_over_wound",     "Tissue over wound",        "first observed tissue/coenosarc over wound",
  "axial_polyp_formation", "Axial polyp formation",   "first observed central axial polyp/hole",
  "wound_smoothed",        "Wound smoothed",          "first observed smoothed wound surface",
  "pigment_over_wound",    "Pigment over wound",      "first observed pigment over wound; scored D0-D7 only",
  "tip_exist",             "Tip exists",              "first observed visible axial tip",
  "tip_extension",         "Tip extension",           "first observed vertical tip extension",
  "new_corallites_on_tip", "New corallites on tip",   "first observed new corallites on the regenerating tip"
)
traits <- trait_meta$trait

long <- micro |>
  select(day, date, treatment, tank, thicket, id, use_photo, notes,
         all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  left_join(trait_meta, by = "trait") |>
  filter(!is.na(expressed))

prop_df <- long |>
  group_by(day, treatment, trait, label) |>
  summarise(
    n = n(),
    n_yes = sum(expressed == 1),
    prop = mean(expressed),
    .groups = "drop"
  )
write_csv(prop_df, file.path(TBL_DIR, "11_microscope_trait_trajectory.csv"))

prop_genet_df <- long |>
  group_by(day, treatment, thicket, trait, label) |>
  summarise(
    n = n(),
    n_yes = sum(expressed == 1),
    prop = mean(expressed),
    .groups = "drop"
  )
write_csv(prop_genet_df,
          file.path(TBL_DIR, "11_microscope_trait_trajectory_by_genet.csv"))

# ---- First-observed event timing ------------------------------------------
compute_events <- function(d, trait_name) {
  d |>
    select(id, treatment, tank, thicket, day, y = all_of(trait_name)) |>
    filter(!is.na(day), !is.na(y)) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      n_visits = n(),
      first_scored_day = min(day),
      last_scored_day = max(day),
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event_lower = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else if (first1 == 1) min(day) else day[first1 - 1]
      },
      event_upper = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) Inf else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(trait = trait_name)
}

events <- map_dfr(traits, ~ compute_events(micro, .x)) |>
  mutate(treatment = factor(treatment, levels = c("28C", "31C"))) |>
  left_join(trait_meta, by = "trait")

event_summary <- events |>
  group_by(trait, label, treatment, thicket) |>
  summarise(
    n_corals = n(),
    n_events = sum(event),
    prop_event = n_events / n_corals,
    median_event_day = if (n_events > 0) median(event_day[event == 1]) else NA_real_,
    iqr_low = if (n_events > 0) quantile(event_day[event == 1], 0.25) else NA_real_,
    iqr_high = if (n_events > 0) quantile(event_day[event == 1], 0.75) else NA_real_,
    last_scored_day = max(last_scored_day),
    .groups = "drop"
  )
write_csv(event_summary, file.path(TBL_DIR, "11_microscope_event_summary.csv"))

endpoint_tests <- map_dfr(traits, function(tr) {
  d <- long |> filter(trait == tr)
  endpoint_day <- max(d$day, na.rm = TRUE)
  e <- d |> filter(day == endpoint_day)
  yes28 <- sum(e$expressed[e$treatment == "28C"] == 1)
  no28 <- sum(e$expressed[e$treatment == "28C"] == 0)
  yes31 <- sum(e$expressed[e$treatment == "31C"] == 1)
  no31 <- sum(e$expressed[e$treatment == "31C"] == 0)
  n28 <- yes28 + no28
  n31 <- yes31 + no31
  mat <- matrix(c(yes28, no28, yes31, no31), nrow = 2, byrow = TRUE,
                dimnames = list(treatment = c("28C", "31C"),
                                state = c("yes", "no")))
  ft <- fisher.test(mat)
  odds31 <- (yes31 + 0.5) / (no31 + 0.5)
  odds28 <- (yes28 + 0.5) / (no28 + 0.5)
  tibble(
    trait = tr,
    endpoint_day = endpoint_day,
    n_28 = n28,
    yes_28 = yes28,
    prop_28 = yes28 / n28,
    n_31 = n31,
    yes_31 = yes31,
    prop_31 = yes31 / n31,
    risk_diff_31_minus_28 = yes31 / n31 - yes28 / n28,
    odds_ratio_31_vs_28 = odds31 / odds28,
    fisher_p = ft$p.value
  )
})

logrank_tests <- map_dfr(traits, function(tr) {
  d <- events |> filter(trait == tr)
  ok <- n_distinct(d$treatment) == 2 && sum(d$event, na.rm = TRUE) > 0
  if (!ok) {
    return(tibble(trait = tr, logrank_chisq = NA_real_,
                  logrank_df = NA_real_, logrank_p = NA_real_))
  }
  fit <- tryCatch(
    survdiff(Surv(event_day, event) ~ treatment, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(trait = tr, logrank_chisq = NA_real_,
                  logrank_df = NA_real_, logrank_p = NA_real_))
  }
  df <- length(fit$n) - 1
  tibble(
    trait = tr,
    logrank_chisq = unname(fit$chisq),
    logrank_df = df,
    logrank_p = pchisq(fit$chisq, df = df, lower.tail = FALSE)
  )
})

event_tests <- endpoint_tests |>
  left_join(logrank_tests, by = "trait") |>
  left_join(trait_meta, by = "trait") |>
  relocate(label, event_interpretation, .after = trait) |>
  mutate(
    test_scope = "exploratory microscope-only validation; n=8 corals per temperature",
    endpoint_test = "Fisher exact test at last scored day for each trait",
    time_to_event_test = "log-rank test on first observed day"
  )
write_csv(event_tests, file.path(TBL_DIR, "11_microscope_event_tests.csv"))

photo_series <- micro |>
  group_by(id, treatment, thicket, tank) |>
  summarise(
    n_scored_days = n_distinct(day),
    n_use_photo_days = sum(use_photo == 1, na.rm = TRUE),
    prop_use_photo_days = n_use_photo_days / n_scored_days,
    first_use_photo_day = if (n_use_photo_days > 0) min(day[use_photo == 1], na.rm = TRUE) else NA_integer_,
    last_use_photo_day = if (n_use_photo_days > 0) max(day[use_photo == 1], na.rm = TRUE) else NA_integer_,
    notes_with_flags = paste(unique(stats::na.omit(notes)), collapse = " | "),
    .groups = "drop"
  ) |>
  arrange(desc(n_use_photo_days), treatment, thicket, id)
write_csv(photo_series,
          file.path(TBL_DIR, "11_microscope_use_photo_series.csv"))

# ---- Figures ---------------------------------------------------------------
p_micro <- ggplot(prop_df, aes(day, prop, colour = treatment,
                               group = treatment)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_wrap(~ label, ncol = 4) +
  labs(x = "Days post-wounding",
       y = "Photo corals expressing trait",
       title = "Microscope wound-healing trajectories",
       subtitle = "Separate photo-only cohort; wounded corals only (n = 8 per temperature)") +
  theme_pub(9)
save_fig(p_micro, "11_microscope_trait_trajectories",
         width = 210, height = 130)

p_micro_genet <- ggplot(prop_genet_df,
                        aes(day, prop, colour = treatment,
                            linetype = thicket,
                            group = interaction(treatment, thicket))) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  geom_point(size = 1.3, alpha = 0.9) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_linetype_manual(values = c(a = "solid", c = "22"),
                        name = "Genet") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_wrap(~ label, ncol = 4) +
  labs(x = "Days post-wounding",
       y = "Photo corals expressing trait",
       title = "Microscope trajectories by genet",
       subtitle = "Photo-only cohort; genets A and C only") +
  theme_pub(9)
save_fig(p_micro_genet, "11b_microscope_trait_trajectories_by_genet",
         width = 210, height = 130)

cat("\n=== Microscope-only wound-healing branch ===\n")
cat("Source: Microscope physio characterization (Drive modified ",
    MICROSCOPE_SOURCE_MODIFIED_UTC, ")\n", sep = "")
cat("Cleaned observations:", nrow(micro), "; corals:", n_distinct(micro$id),
    "; days:", paste(range(micro$day), collapse = "-"), "\n")
print(as.data.frame(event_tests |>
  filter(trait %in% c("tissue_over_wound", "wound_smoothed",
                      "new_corallites_on_tip")) |>
  transmute(trait, endpoint_day,
            `28C` = sprintf("%d/%d", yes_28, n_28),
            `31C` = sprintf("%d/%d", yes_31, n_31),
            risk_diff_31_minus_28 = round(risk_diff_31_minus_28, 3),
            fisher_p = signif(fisher_p, 3),
            logrank_p = signif(logrank_p, 3))))

cat("Wrote microscope_physio_clean.rds, microscope CSV summaries, and ",
    "11_microscope_trait_trajectories.{pdf,png}\n", sep = "")
