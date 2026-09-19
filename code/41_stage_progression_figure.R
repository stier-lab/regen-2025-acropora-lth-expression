# =============================================================================
# Purpose: One honest display of healing/regeneration milestone timing.
#
# What & why: scripts 36, 39 and 40 offer five designs for the same table. Three
# of them mislead at this sample size, for reasons worth stating in code so they
# are not reintroduced:
#
#   * 40_stage_mean_review plots a MEAN FIRST-OBSERVATION DAY among fragments
#     that reached the stage. For "New radial corallites bud" at 31 C that mean
#     is computed on 4 of 12 fragments, and it lands ~1 day right of the 28 C
#     mean. The eye reads "heated corals were a day slower"; the finding is that
#     two thirds of them never got there at all. A central tendency of the
#     survivors is not an estimate of anything when censoring is that heavy.
#   * 39a_stage_timing_violin_option draws kernel densities from <=12 points —
#     and, again, only from the fragments that reached the stage, so the panel
#     carrying the result summarises the biased third that did.
#   * 39c_stage_timing_heatmap_option encodes "percent reached" as a blue ramp,
#     while blue means "28 C" in every other figure in this repo.
#
# The design here keeps every fragment visible and every fragment in the
# denominator: individual first-observation days on the left, an explicit
# "never reached" zone on the right of the same axis, and the proportion
# reached with a binomial interval in the margin. A median tick is drawn ONLY
# where most of the group reached the stage (see MEDIAN_MIN_REACHED), so the
# figure cannot repeat script 40's error.
#
# Stage 1 comes from the 16-fragment microscope cohort and stages 2-6 from the
# 24-fragment main cohort, so the two are separated by a gap and labelled rather
# than left to a footnote.
#
# NOTE ON ORDERING: the milestones are numbered but are NOT a strict nested
# progression. Three of 24 fragments record a later milestone on an earlier day
# than an earlier one, one reached a later milestone without an earlier one, and
# at 28 C stage 5 is 11/12 while stage 6 is 12/12. Do not redraw this as a
# funnel or attrition cascade: that would impose an ordering the data do not
# obey. They are milestones that usually, not always, occur in order.
#
# Input:   output/tables/36_molly_stage_timing.csv
# Output:  figures/41_stage_progression.{pdf,png}
#          output/tables/41_stage_reached_proportions.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

## A median first-observation day is only shown when at least this fraction of
## the group reached the milestone. Below it, the reachers are a selected subset
## and a central tendency would mislead.
MEDIAN_MIN_REACHED <- 0.8

FINAL_DAY <- 15
NR_X      <- FINAL_DAY + 2.4   # x position of the "never reached" zone
NR_GAP    <- FINAL_DAY + 1.2   # where the panel break is drawn

stage_events <- readr::read_csv(
  file.path(TBL_DIR, "36_molly_stage_timing.csv"), show_col_types = FALSE
) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    treatment_label = factor(
      dplyr::recode(as.character(treatment), `28C` = "28 °C", `31C` = "31 °C"),
      levels = c("28 °C", "31 °C")),
    stage_order = as.integer(stage_order),
    reached = event_status == "Reached",
    ## short, verb-first labels; the Healing/Regeneration split is carried by the
    ## cohort annotation instead of being repeated in five row labels
    stage_short = dplyr::recode(
      as.character(stage),
      `Coenosarc barrier re-established` = "Tissue covers wound",
      `Axial polyp forms`                = "Main polyp forms",
      `Wound surface remodels`           = "Wound surface remodels",
      `Axial corallite or tip forms`     = "Tip forms",
      `Axial corallite or tip extends`   = "Tip extends",
      `New radial corallites bud`        = "New skeletal cups bud")
  )

## row order: milestone 1 at the top, 6 at the bottom, so the eye travels down
## the sequence and lands on the result
lev <- stage_events |> distinct(stage_order, stage_short) |> arrange(desc(stage_order))
stage_events <- stage_events |>
  mutate(stage_short = factor(stage_short, levels = lev$stage_short))

## ---- proportion reached, with a binomial interval -------------------------
prop <- stage_events |>
  group_by(stage_order, stage_short, treatment, treatment_label) |>
  summarise(n = dplyr::n(), k = sum(reached), .groups = "drop") |>
  rowwise() |>
  mutate(p = k / n,
         lo = stats::binom.test(k, n)$conf.int[1],
         hi = stats::binom.test(k, n)$conf.int[2]) |>
  ungroup() |>
  mutate(label = paste0(k, "/", n))

readr::write_csv(prop |> select(stage_order, stage_short, treatment, n, k, p, lo, hi),
                 file.path(TBL_DIR, "41_stage_reached_proportions.csv"))

## ---- medians, only where the group mostly reached the milestone -----------
med <- stage_events |>
  filter(reached) |>
  group_by(stage_order, stage_short, treatment, treatment_label) |>
  summarise(med = median(first_yes_day), .groups = "drop") |>
  left_join(prop |> select(stage_order, treatment, p), by = c("stage_order", "treatment")) |>
  filter(p >= MEDIAN_MIN_REACHED)

## ---- point positions: two sub-rows per milestone, fragments stacked --------
## The stack is centred on its sub-row and SCALED to fit inside it. A fixed step
## overflowed into the neighbouring milestone wherever many fragments shared a
## day, which made rows unreadable at the busiest days (6, 7 and 12).
dodge    <- c(`28 °C` = 0.20, `31 °C` = -0.20)
SUB_HALF <- 0.15                      # half-height available to one sub-row
pts <- stage_events |>
  mutate(x = ifelse(reached, first_yes_day, NR_X)) |>
  group_by(stage_short, treatment_label, x) |>
  mutate(k = dplyr::n(),
         step = pmin(0.05, (2 * SUB_HALF) / pmax(k - 1, 1)),
         stack = (row_number() - 1 - (k - 1) / 2) * step) |>
  ungroup() |>
  mutate(yoff = unname(dodge[as.character(treatment_label)]),
         y = as.numeric(stage_short) + yoff + stack)

nlev  <- length(levels(stage_events$stage_short))
bands <- tibble(y = seq_len(nlev)) |> filter(y %% 2 == 0)

p_main <- ggplot() +
  ## alternating milestone bands. annotate() with explicit bounds, because
  ## geom_tile() with a width wider than the scale was silently clipped away.
  annotate("rect", xmin = -0.5, xmax = NR_X + 1.3,
           ymin = bands$y - 0.46, ymax = bands$y + 0.46,
           fill = "grey96", colour = NA) +
  ## the "never reached" zone is part of the same axis, separated by a rule
  annotate("rect", xmin = NR_GAP + 0.25, xmax = NR_X + 1.3,
           ymin = 0.4, ymax = nlev + 0.6, fill = "grey90", colour = NA) +
  annotate("segment", x = NR_GAP + 0.25, xend = NR_GAP + 0.25, y = 0.4,
           yend = nlev + 0.6, colour = "grey50", linetype = "dotted") +
  ## median drawn first so points sit on top of it
  geom_segment(data = med,
               aes(x = med, xend = med,
                   y = as.numeric(stage_short) + unname(dodge[as.character(treatment_label)]) - 0.24,
                   yend = as.numeric(stage_short) + unname(dodge[as.character(treatment_label)]) + 0.24,
                   colour = treatment_label),
               linewidth = 1.6, alpha = 0.40, lineend = "butt", show.legend = FALSE) +
  geom_point(data = filter(pts, reached),
             aes(x, y, colour = treatment_label), size = 1.4, alpha = 0.95) +
  geom_point(data = filter(pts, !reached),
             aes(x, y, colour = treatment_label), shape = 4, size = 1.8,
             stroke = 0.75, show.legend = FALSE) +
  scale_colour_manual(values = unname(PAL_TEMP), name = NULL,
                      guide = guide_legend(override.aes = list(size = 2.4))) +
  scale_y_continuous(breaks = seq_along(levels(stage_events$stage_short)),
                     labels = levels(stage_events$stage_short),
                     expand = expansion(add = c(0.45, 0.55))) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15),
                     limits = c(-0.5, NR_X + 1.3),
                     expand = expansion(add = c(0, 0)),
                     sec.axis = dup_axis(breaks = NR_X, labels = "never\nreached", name = NULL)) +
  labs(x = "Day the milestone was first seen", y = NULL) +
  theme_pub(10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(hjust = 0))

p_prop <- ggplot(prop, aes(p, as.numeric(stage_short) +
                             unname(dodge[as.character(treatment_label)]),
                           colour = treatment_label)) +
  annotate("rect", xmin = -0.02, xmax = 1.45,
           ymin = bands$y - 0.46, ymax = bands$y + 0.46,
           fill = "grey96", colour = NA) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.5, alpha = 0.75) +
  geom_point(size = 1.8) +
  geom_text(aes(x = 1.12, label = label), size = 2.6, hjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = unname(PAL_TEMP), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "50", "100%"),
                     limits = c(-0.02, 1.45), expand = expansion(add = c(0, 0))) +
  scale_y_continuous(limits = c(0.4, length(levels(stage_events$stage_short)) + 0.6),
                     expand = expansion(add = c(0, 0))) +
  labs(x = "Reached by day 15", y = NULL) +
  theme_pub(10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- (p_main | p_prop) +
  patchwork::plot_layout(widths = c(3.1, 1), guides = "collect") +
  patchwork::plot_annotation(
    title = "Heat blocked the last regeneration milestone, not the earlier ones",
    subtitle = paste0(
      "Each point is one wounded fragment. Crosses are fragments that never reached the milestone by day 15 and stay in the denominator.\n",
      "Vertical ticks mark the median day, drawn only where at least ",
      round(MEDIAN_MIN_REACHED * 100), "% of the group reached the milestone."),
    caption = paste0(
      "Top row is the 16-fragment microscope cohort (8 per temperature); the other five are the main cohort (12 per temperature).\n",
      "Bars are 95% binomial intervals. Milestones are numbered but are not strictly nested: 3 of 24 fragments recorded one out of order."),
    theme = theme_pub(10) +
      theme(plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 8.2, colour = "grey30"),
            plot.caption = element_text(size = 7.2, colour = "grey40", hjust = 0)))

save_fig(fig, "41_stage_progression", width = 190, height = 118)

cat("\n== Reached by day 15 ==\n")
print(as.data.frame(prop |> select(stage_order, stage_short, treatment, label, p)),
      row.names = FALSE, digits = 2)
cat("\nMedians suppressed (reached < ", round(MEDIAN_MIN_REACHED * 100), "%):\n", sep = "")
print(as.data.frame(anti_join(prop, med, by = c("stage_order", "treatment")) |>
                      select(stage_short, treatment, label)), row.names = FALSE)
cat("\nwrote figures/41_stage_progression.{pdf,png}\n")
