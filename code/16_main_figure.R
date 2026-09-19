# =============================================================================
# Purpose: Build the publication Figure 1 — a four-panel narrative figure
#          for the LTH manuscript.
#          Panel A: Apex tank temperature trace (treatment validation)
#          Panel B: PAM Fv/Fm trajectory (whole-colony stress)
#          Panel C: Kaplan-Meier curves for the two diagnostic morphology
#                   traits (wound smoothed vs new corallites — closure
#                   preserved, regeneration impaired)
#          Panel D: PCA biplot (multivariate summary)
#
# What & why: this is the manuscript's single multi-panel figure. Each panel is
#   built as its own self-contained ggplot object (pA, pB, pC, pD), then patchwork
#   assembles them into a 2×2 grid at the end. The narrative reads left-to-right,
#   top-to-bottom: (A) show the treatment worked — the tanks held 28 vs 31 °C;
#   (B) show whole-colony physiological stress diverge over time (Fv/Fm); (C) the
#   mechanism — heat lets wounds close (left KM) but blocks new-corallite
#   regeneration (right KM); (D) summarize in multivariate space. Assembly
#   conventions used throughout:
#   each panel carries a letter via labs(tag = ...) (NOT a subtitle — panel
#   descriptions live in the caption); legends are suppressed per-panel with
#   guide = "none" and groups are instead DIRECT-LABELLED inside the panel
#   (annotate()) or collected once at the bottom via guides = "collect"; every
#   panel shares the same Okabe-Ito blue/orange (28 °C = #56B4E9, 31 °C = #D55E00)
#   and theme_pub(9) so the grid looks like one figure, not four. The final
#   composition step controls relative panel widths/heights and the figure title.
# Input:   data/processed/{apex_temperature_daily,pam_clean,physio_clean,
#                          coral_physio_wide}.rds
# Output:  figures/16_manuscript_fig1.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(survival))

# One cleaned table per panel — A pulls from the Apex temperature loggers,
# B from PAM fluorometry, C from the morphology scores, D from the wide
# one-row-per-coral physiology summary.
apex   <- readRDS(file.path(DATA_PROC, "apex_temperature_daily.rds"))
pam    <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
physio <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
wide   <- readRDS(file.path(DATA_PROC, "coral_physio_wide.rds"))

# Design n for the figure summary line, computed from the data (the old hardcoded
# "192 fragments" matched neither the 48 corals nor the 768 coral-day records).
n_corals <- dplyr::n_distinct(physio$id)
n_wound  <- physio |> filter(wound == "yes") |> distinct(id) |> nrow()
n_genet  <- dplyr::n_distinct(physio$thicket)
n_tank   <- dplyr::n_distinct(physio$tank)

# ---- Panel A: tank temperature (treatment validation) ---------------------
# Reduce the logger probes to the eight experimental tanks, fix unit drift, and
# tag each tank with its treatment so the trace shows 28 vs 31 °C separation.
temp_pa <- apex |>
  mutate(probe = str_squish(probe)) |>
  filter(str_detect(probe, "^Temp\\d+$"),
         is.finite(value_mean), value_mean > 0) |>
  mutate(
    tank = as.integer(str_extract(probe, "\\d+")),
    # apex_temperature_daily.rds is already in °C (the °F->°C conversion happens at
    # the raw-reading level in code/08, before aggregation), so no conversion here.
    value_c = value_mean,
    treatment = tank_treatment(tank)   # shared tank->treatment map (00_setup.R)
  ) |>
  filter(!is.na(treatment), value_c > 15, value_c < 40,
         between(date, as_date("2025-05-28"), as_date("2025-06-19")))

# With only two groups we direct-label instead of using a legend (CLAUDE.md
# tiered legend rule). Take one anchor point per treatment — the last date of one
# tank — to position the "28 °C"/"31 °C" text annotations placed below.
last_28 <- temp_pa |> filter(treatment == "28C") |>
  group_by(tank) |> slice_max(date, n = 1) |> ungroup() |>
  slice(1)
last_31 <- temp_pa |> filter(treatment == "31C") |>
  group_by(tank) |> slice_max(date, n = 1) |> ungroup() |>
  slice(1)

# group = tank draws one line per tank; colour encodes treatment. Dashed grey
# guides at 28/31 mark the target temperatures.
pA <- ggplot(temp_pa, aes(date, value_c, group = tank, colour = treatment)) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(size = 0.7, alpha = 0.7) +
  annotate("text", x = last_28$date, y = 27.5, label = "28 °C",
           hjust = 1.1, vjust = 1, size = 3, colour = "#0072B2",
           fontface = "bold") +
  annotate("text", x = last_31$date, y = 32, label = "31 °C",
           hjust = 1.1, vjust = -0.1, size = 3, colour = "#D55E00",
           fontface = "bold") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = NULL, guide = "none") +  # direct-labeled
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  labs(x = NULL, y = "Tank T (°C)",
       tag = "A") +   # panel descriptions belong in the figure caption, not subtitles
  theme_pub(9)

# ---- Panel B: PAM trajectory (whole-colony stress) ------------------------
# Collapse to group means ± standard error per day × treatment × wound, so the
# panel shows four mean trajectories with shaded SE ribbons.
pam_df <- pam |>
  group_by(day, treatment, wound) |>
  summarise(m = mean(fv_fm, na.rm = TRUE),
            se = sd(fv_fm, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Ribbon (mean±SE) + line + points, grouped by treatment×wound. Treatment is
# direct-labelled at the right edge (annotate) with clip = "off" so the labels
# can sit just outside the plotting area; the wound shape legend is collected
# later from panel D so it appears only once.
pB <- ggplot(pam_df, aes(day, m, colour = treatment,
                          fill = treatment, shape = wound)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se,
                  group = interaction(treatment, wound)),
              alpha = 0.15, colour = NA) +
  geom_line(aes(group = interaction(treatment, wound)),
            linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  # Embed small treatment labels at right edge instead of using a legend
  annotate("text", x = max(pam_df$day) + 0.5, y = 0.575,
           label = "31 °C", colour = "#D55E00", size = 2.7,
           fontface = "bold", hjust = 0) +
  annotate("text", x = max(pam_df$day) + 0.5, y = 0.675,
           label = "28 °C", colour = "#0072B2", size = 2.7,
           fontface = "bold", hjust = 0) +
  coord_cartesian(xlim = c(min(pam_df$day), max(pam_df$day) + 2),
                  clip = "off") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      guide = "none") +
  scale_fill_manual(values   = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(no = 16, yes = 17), name = "Wound",
                     guide = "none") +  # collected from panel D
  labs(x = "Days post-wounding", y = expression(F[v]/F[m]),
       tag = "B") +
  theme_pub(9)

# ---- Panel C: KM curves for the two diagnostic traits --------------------
# The mechanism, reduced from the full code/14 analysis to the two
# contrasting milestones: closure (wound_smoothed) vs regeneration
# (new_corallites_on_tip). Recompute first-observed event days and Kaplan-
# Meier curves inline here (same logic as code/14) so this figure is
# self-contained. See code/14 for the censoring/KM explanation.
focus_traits <- c(wound_smoothed = "Wound smoothed (closure)",
                  new_corallites_on_tip = "New corallites on tip (regeneration)")

# Collapse each coral's 0/1 history to its first-observed event day (or last day
# seen if never reached = right-censored), one record per coral × trait.
ph_events <- map_dfr(names(focus_traits), function(tr) {
  physio |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[tr]]) |>
    group_by(id, treatment) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(trait = unname(focus_traits[tr]))
})

# Kaplan-Meier computed manually: cumprod of (1 - events/at-risk) over event days, plotted
# as cumulative % expressing (1 - survival), seeded at day 0 = 0% expressed.
km_curves <- ph_events |>
  group_by(trait, treatment) |>
  group_modify(\(d, k) {
    d <- arrange(d, event_day)
    times <- sort(unique(d$event_day))
    n_at_risk <- sapply(times, \(t) sum(d$event_day >= t))
    n_events  <- sapply(times, \(t) sum(d$event_day == t & d$event == 1))
    surv      <- cumprod(1 - n_events / pmax(n_at_risk, 1))
    tibble(day = c(0, times), cum_event = 1 - c(1, surv))
  }) |> ungroup() |>
  mutate(trait = factor(trait, levels = unname(focus_traits)))

# Two side-by-side KM staircases (facet per trait); the colour legend is
# suppressed here and collected once at the bottom of the assembled figure.
pC <- ggplot(km_curves, aes(day, cum_event,
                            colour = treatment, group = treatment)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ trait, ncol = 2) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = NULL, guide = "none") +  # collected below
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Days post-wounding",
       y = "Corals expressing trait",
       tag = "C") +
  theme_pub(9)

# ---- Panel D: PCA biplot (multivariate summary) ---------------------------
# Compress four end-of-experiment response variables into two principal axes.
# pca_in keeps only complete cases; the parallel `keep`/`groups` carry the
# matching treatment/wound labels for colouring points. center+scale. = TRUE
# standardises the variables so no single one dominates by virtue of its units.
pca_vars <- c("pam_end", "color_end", "growth_pct", "zoox_end")
pretty   <- c(pam_end      = "Fv/Fm",
              color_end    = "Color",
              growth_pct = "Growth",
              zoox_end     = "Symbionts")
pca_in   <- wide[, pca_vars, drop = FALSE] |> drop_na()
keep     <- complete.cases(wide[, pca_vars])
groups   <- wide[keep, c("treatment", "wound")]

pca   <- prcomp(pca_in, center = TRUE, scale. = TRUE)
var_e <- summary(pca)$importance[2, ] * 100   # % variance per PC (for axis labels)
scores <- as.data.frame(pca$x) |> bind_cols(groups)   # coral positions in PC space
# Loadings = variable arrows. Scaled ×2.7 only so the arrows are visible against
# the point cloud (a biplot's arrow length is arbitrary; this is cosmetic).
load_a <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(label = unname(pretty[variable]),
         across(c(PC1, PC2), \(x) x * 2.7))

# Biplot layers, drawn back-to-front: zero gridlines -> 68% treatment ellipses
# (one SD-equivalent cloud) -> coral points (colour = treatment, shape = wound)
# -> loading arrows -> ggrepel labels nudged off the arrow tips so text never
# overlaps. Wound shapes are direct-labelled with a small annotate() legend.
pD <- ggplot(scores, aes(PC1, PC2, colour = treatment, shape = wound)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  stat_ellipse(aes(group = treatment, fill = treatment),
               geom = "polygon", level = 0.68,
               alpha = 0.10, colour = NA, show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.9, stroke = 0.4) +
  geom_segment(data = load_a,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(1.5, "mm")),
               colour = "grey25", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = load_a,
            aes(PC1 * 1.13, PC2 * 1.13, label = label),
            inherit.aes = FALSE,
            size = 3, colour = "grey15", fontface = "bold",
            box.padding = 0.2, point.padding = 0.5,
            min.segment.length = Inf, seed = 42) +
  # Direct-label wound status inside panel
  annotate("text", x = min(scores$PC1) + 0.2, y = max(scores$PC2) - 0.3,
           label = "● unwounded   ▲ wounded", colour = "grey25",
           size = 2.6, hjust = 0) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      guide = "none") +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(no = 16, yes = 17), guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%)", var_e[1]),
       y = sprintf("PC2 (%.0f%%)", var_e[2]),
       tag = "D") +
  theme_pub(9)

# ---- Compose ---------------------------------------------------------------
# Assemble the four panels with patchwork. `+` places panels side by side, `/`
# stacks rows. Top row A|B is equal width; bottom row gives the two-facet KM
# panel C extra room (1.4 vs 1) so its facets are not compressed by the biplot.
# heights = c(0.9, 1.1) makes the bottom row a bit taller; guides = "collect"
# gathers every surviving legend (the wound shapes) into ONE shared legend.
top_row <- pA + pB + patchwork::plot_layout(widths = c(1, 1))
bot_row <- pC + pD + patchwork::plot_layout(widths = c(1.4, 1))
fig1    <- (top_row / bot_row) +
  patchwork::plot_layout(heights = c(0.9, 1.1), guides = "collect") +
  patchwork::plot_annotation(
    title = expression("Heat compromises photochemistry, growth, and regeneration in"~italic(A.~pulchra)),
    subtitle = sprintf("n = %d corals (%d wounded) · %d genets · %d tanks · 28 vs 31 °C",
                       n_corals, n_wound, n_genet, n_tank),
    theme = theme(
      plot.title    = element_text(size = 11, face = "bold",
                                   margin = margin(b = 4, t = 4)),
      plot.subtitle = element_text(size = 9, colour = "grey30",
                                   margin = margin(b = 10))
    )
  ) &
  # `&` applies this theme to EVERY panel at once (patchwork operator), so the
  # single collected legend sits at the bottom across the whole figure.
  theme(legend.position = "bottom",
        legend.box      = "horizontal",
        legend.margin   = margin(t = 4))

# save_fig() writes the vector PDF (for the manuscript) + a PNG preview at the
# given size in mm; double-column width since this is the main figure.
save_fig(fig1, "16_manuscript_fig1", width = 220, height = 175)

cat("Wrote 16_manuscript_fig1.{pdf,png}\n")
