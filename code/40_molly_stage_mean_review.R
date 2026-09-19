# Descriptive companion and scoring audit requested by Molly on 2026-09-15.
# Retain non-reachers in denominators, never assign them a completion day.
source(here::here("code", "00_setup.R"))
events <- read_csv(file.path(TBL_DIR, "36_molly_stage_timing.csv"),
                   show_col_types = FALSE)
stage_names <- c("Living tissue covers wound", "Main polyp forms",
                 "Wound surface remodels", "Tip forms", "Tip extends",
                 "New skeletal cups bud")
means <- events |>
  group_by(stage_order, treatment) |>
  summarise(n_total = n(), n_reached = sum(event),
            mean_day = if (any(event)) mean(first_yes_day[event]) else NA_real_,
            q25 = if (any(event)) quantile(first_yes_day[event], .25) else NA_real_,
            q75 = if (any(event)) quantile(first_yes_day[event], .75) else NA_real_,
            .groups = "drop") |>
  mutate(stage = stage_names[stage_order],
         coverage = if_else(n_reached == n_total, "All fragments", "Reachers only"),
         row = 7 - stage_order + if_else(treatment == "28C", .14, -.14),
         count_label = paste0(n_reached, "/", n_total))
stopifnot(nrow(means) == 12, all(means$n_reached <= means$n_total))
write_csv(means, file.path(TBL_DIR, "40_stage_mean_review.csv"))
p <- ggplot(means, aes(mean_day, row, colour = treatment)) +
  geom_segment(aes(x = q25, xend = q75, yend = row), linewidth = .8) +
  geom_point(aes(shape = coverage), size = 2.8, stroke = .9) +
  geom_text(aes(x = 16.5, label = count_label), size = 3, show.legend = FALSE) +
  annotate("text", x = 16.5, y = 6.65, label = "Reached", size = 3) +
  scale_colour_manual(values = c("28C" = "#56B4E9", "31C" = "#D55E00"),
                      labels = c("28C" = "28 C", "31C" = "31 C"), name = NULL) +
  scale_shape_manual(values = c("All fragments" = 16, "Reachers only" = 1), name = NULL) +
  scale_y_continuous(breaks = 6:1, labels = stage_names, limits = c(.5, 6.8)) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15), limits = c(.5, 17.4)) +
  labs(x = "Average first-observation day among fragments that reached the stage", y = NULL,
       title = "Average days, shown alongside the number reaching each stage",
       subtitle = "Open circles leave out fragments that had not reached the stage by Day 15.",
       caption = "Lines show the middle half of observed days, not uncertainty in the mean.\nTissue cover: separate microscope experiment (8 per temperature); other stages: main experiment (12 per temperature).") +
  theme_pub(9) + theme(legend.position = "bottom")
save_fig(p, "40_stage_mean_review", width = 210, height = 125)

sequence <- events |> filter(stage_order %in% c(5, 6)) |>
  select(id, treatment, tank, thicket, stage_order, event, first_yes_day) |>
  pivot_wider(names_from = stage_order, values_from = c(event, first_yes_day)) |>
  filter(event_6 & (!event_5 | first_yes_day_6 < first_yes_day_5))
write_csv(sequence, file.path(TBL_DIR, "40_stage_sequence_review.csv"))
scores <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  semi_join(sequence, by = "id") |>
  select(id, treatment, tank, thicket, day, tip_exist, tip_extension,
         new_corallites_on_tip, notes) |> arrange(id, day)
write_csv(scores, file.path(TBL_DIR, "40_stage_sequence_original_scores.csv"))
print(sequence)
