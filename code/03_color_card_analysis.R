# =============================================================================
# Purpose: Color-card pigmentation (Siebeck et al. 2006 D-scale, D1-D6)
#          - Convert D-scale text to ordinal numeric (D1=1, ..., D6=6)
#          - Trajectory plot + paling/hole binary tallies
#          - Ordinal cumulative-link mixed model (clmm) for color
#
# What & why: corals were scored against the Siebeck CoralWatch "Coral Health
#   Chart" — a laminated reference card with six brown/tan shades (D1 = palest,
#   D6 = darkest). Shade tracks symbiont density and chlorophyll, so a DROP on
#   the D-scale ("paling") is an early, visible sign of bleaching stress before
#   the coral goes white. The question here is whether the 31 °C heat treatment
#   (and wounding) caused corals to pale over the experiment. This script does
#   the data prep, an end-of-experiment summary, and the trajectory figure; the
#   ordinal model (clmm) that respects the discrete, ranked nature of the scale
#   is in code/12b_color_clmm_robustness.R. Pigmentation is the whole-colony
#   stress signal that complements the symbiont/chlorophyll lab assays (code/06)
#   and the Fv/Fm photochemistry (code/02).
# Input:   data/raw/color_card/data.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/color_clean.rds
#          figures/03_color_trajectory.{pdf,png}
#          output/tables/03_color_end_proportions.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_RAW, DATA_PROC, TBL_DIR), the
# theme_pub() theme, the PAL_WOUND palette, and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Helper: parse the D-scale text into a number --------------------------
# Convert a Siebeck D-scale text entry to numeric. Single scores ("D5") map to
# their integer; split scores ("D3/D4") are averaged to the midpoint (3.5),
# following Molly's original convention (code/archive/molly_original/). This
# affects 40 of ~962 scored observations (D1/D2, D2/D3, D3/D4); taking only the
# first digit would bias those downward by 0.5.
# Mechanics: str_extract_all pulls every run of digits from the string; vapply
# returns one numeric per entry (NA when no digit was found, e.g. blank cells).
convert_color <- function(x) {
  nums <- stringr::str_extract_all(x, "\\d+")
  vapply(nums, function(v) if (length(v) == 0) NA_real_ else mean(as.numeric(v)),
         numeric(1))
}

# ---- Load ------------------------------------------------------------------
# clean_names() standardises the spreadsheet headers to snake_case.
cc_raw <- read_csv(file.path(DATA_RAW, "color_card", "data.csv"),
                   show_col_types = FALSE) |>
  janitor::clean_names()

# ---- Clean -----------------------------------------------------------------
# Type-coerce every column, tidy the categorical fields, and turn the D-scale
# text into the numeric color_num used downstream.
cc <- cc_raw |>
  # The raw header is "wounded"; rename to the project-wide "wound". thicket may
  # have varied capitalisation, so match it by prefix.
  rename(thicket = matches("^thicket"),
         wound   = wounded) |>
  mutate(
    date      = as_date(date),
    day       = as.integer(day),                 # day relative to wounding (D0)
    # Treatment recorded as set-point (28/31 °C); 28C is the reference level.
    treatment = factor(as.integer(treatment), levels = c(28, 31),
                       labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),  # genet ID; tidy case/whitespace
    id        = as.integer(id),                  # unique coral fragment ID
    wound     = factor(wound, levels = c("no", "yes")),  # "no" = reference level
    health_status = str_to_lower(str_squish(health_status)),
    # Paling and hole-at-center are yes/no observations; encode as 2-level
    # factors with "no" first so "yes" is the event of interest.
    paling    = factor(str_to_lower(str_squish(paling)),
                       levels = c("no", "yes")),
    hole_at_center = factor(str_to_lower(str_squish(hole_at_center)),
                             levels = c("no", "yes")),
    # Color: "D5" -> 5, "D3/D4" -> 3.5 (split scores averaged); NA/"" -> missing
    color_num = convert_color(color)
  ) |>
  filter(!is.na(color_num))            # drop rows with no usable color score

# Save the cleaned, one-row-per-observation table for downstream scripts.
saveRDS(cc, file.path(DATA_PROC, "color_clean.rds"))

# ---- End-of-experiment proportions ----------------------------------------
# Snapshot of the final timepoint: how many corals in each treatment×wound cell
# were still alive / had paled, and the mean color score. This is the summary
# underlying the trajectory plot and the model (12b).
end_day <- max(cc$day, na.rm = TRUE)             # last day with color data
end_props <- cc |>
  filter(day == end_day) |>
  group_by(treatment, wound) |>
  summarise(
    n_total       = n(),
    n_alive       = sum(health_status == "alive", na.rm = TRUE),
    n_paled       = sum(paling == "yes", na.rm = TRUE),
    mean_color    = mean(color_num, na.rm = TRUE),
    sd_color      = sd(color_num,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(prop_paled = n_paled / n_total)

write_csv(end_props, file.path(TBL_DIR, "03_color_end_proportions.csv"))

# ---- Trajectory plot -------------------------------------------------------
# Collapse to one mean color (± standard error) per day × treatment × wound cell.
traj <- cc |>
  group_by(day, treatment, wound) |>
  summarise(
    mean_color = mean(color_num, na.rm = TRUE),
    se = sd(color_num, na.rm = TRUE) / sqrt(n()),   # SE = SD / sqrt(n)
    n = n(),
    .groups = "drop"
  )

# One panel per temperature; line per wound status with a ±1 SE ribbon. Dotted
# vertical line at day 0 marks the wounding event.
p_color <- ggplot(traj, aes(day, mean_color,
                            colour = wound, fill = wound, group = wound)) +
  geom_ribbon(aes(ymin = mean_color - se, ymax = mean_color + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  # Reverse the y-axis so the DARKEST score (D6, healthy) sits at the bottom and
  # Lower, paler scores appear higher on this reversed axis.
  # Limits padded beyond 1-6 for visual margin.
  scale_y_reverse(breaks = 1:6, limits = c(6.2, 0.8)) +  # darker = higher D, plot 1 on top
  scale_colour_manual(values = PAL_WOUND,
                      name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND,
                    guide = "none") +
  labs(x = "Day relative to wounding (D0)",
       y = "Color score (Siebeck D-scale; lower = paler)",
       title = "Pigmentation trajectory",
       subtitle = "Mean ± 1 SE; paler scores appear higher on the reversed axis") +
  theme_pub(10)

save_fig(p_color, "03_color_trajectory", width = 170, height = 90)

cat("Color score: end-of-experiment proportions written.\n")
print(end_props)
