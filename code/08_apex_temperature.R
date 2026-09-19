# Rebuild temperatures from original and recovered Apex exports.
# Clock/source conventions: docs/provenance/molly_followup_2026-09-11.md.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(xml2))
recovered <- file.path(DATA_RAW, "apex", "recovered_2026-09-10")
xml_files <- list.files(file.path(DATA_RAW, "apex"), pattern = "\\.xml$",
                        recursive = TRUE, full.names = TRUE)
stopifnot(length(xml_files) > 0)
parse_failures <- list()

# Preserve displayed logger times to match Molly's date-based assignments.
# UTC is a storage convention here, not a claim about the actual clock zone.
parse_apex <- function(path) {
  doc <- tryCatch(read_xml(path), error = function(e) {
    parse_failures[[path]] <<- tibble(source_file = basename(path),
                                    reason = conditionMessage(e))
    message("Excluded malformed XML: ", basename(path))
    NULL
  })
  if (is.null(doc)) return(tibble())
  probes <- xml_find_all(doc, "//record/probe[type='Temp']")
  tibble(
    source_path = path,
    serial = xml_text(xml_find_first(doc, "//serial")),
    logger_timezone = xml_text(xml_find_first(doc, "//timezone")),
    datetime = mdy_hms(xml_text(xml_find_first(probes, "../date")), tz = "UTC"),
    probe = str_squish(xml_text(xml_find_first(probes, "./name"))),
    raw_value = as.numeric(xml_text(xml_find_first(probes, "./value")))
  ) |> filter(str_detect(probe, "^Temp[0-9]+$")) |>
    mutate(value = if_else(raw_value > 60, (raw_value - 32) * 5 / 9, raw_value))
}
raw <- map_dfr(xml_files, parse_apex)
write_csv(bind_rows(parse_failures), file.path(TBL_DIR, "08_apex_parse_failures.csv"))
stopifnot(!anyNA(raw$datetime), !anyNA(raw$value))
manifest <- raw |> group_by(source_path, serial, logger_timezone) |>
  summarise(first_record = min(datetime), last_record = max(datetime),
            n_readings = n(), .groups = "drop") |>
  mutate(md5 = unname(tools::md5sum(source_path)),
         source_path = sub(paste0(DATA_RAW, "/"), "", source_path, fixed = TRUE))
write_csv(manifest, file.path(TBL_DIR, "08_apex_source_manifest.csv"))

# Count repeated exports once; never resolve conflicting values by row order.
duplicates <- raw |> group_by(serial, datetime, probe) |>
  summarise(n_copies = n(), n_values = n_distinct(value), .groups = "drop") |>
  filter(n_copies > 1)
write_csv(duplicates, file.path(TBL_DIR, "08_apex_duplicate_audit.csv"))
stopifnot(all(duplicates$n_values == 1), n_distinct(raw$serial) == 1)
readings <- raw |> distinct(serial, datetime, probe, value) |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+")),
         valid = is.finite(value) & value > 15 & value < 40)
write_csv(filter(readings, !valid), file.path(TBL_DIR, "08_apex_excluded_readings.csv"))
readings <- filter(readings, valid)

# Independently cross-check the compiled sheet without adding duplicate records.
compiled <- read_csv(file.path(recovered, "compiled_temperature.csv"),
                     show_col_types = FALSE) |>
  transmute(datetime = mdy_hms(.data[["/record/date"]], tz = "UTC"),
            probe = .data[["/record/probe/name"]],
            raw_value = .data[["/record/probe/value"]],
            compiled_c = if_else(raw_value > 60, (raw_value - 32) * 5 / 9, raw_value))
compiled_check <- compiled |>
  left_join(select(readings, datetime, probe, value),
            by = c("datetime", "probe"), relationship = "many-to-one") |>
  mutate(matches_xml = !is.na(value) & abs(compiled_c - value) < 1e-8)
write_csv(tibble(n_compiled = nrow(compiled_check),
                 n_duplicate_rows = sum(duplicated(compiled)),
                 n_matches_xml = sum(compiled_check$matches_xml),
                 n_unmatched = sum(!compiled_check$matches_xml)),
          file.path(TBL_DIR, "08_apex_compiled_crosscheck.csv"))
stopifnot(all(compiled_check$matches_xml))

# Check whether the combined sheet alone covers the assigned study period.
# Export record IDs are not part of the measurement key.
compiled_unique <- compiled |> distinct(datetime, probe, compiled_c)
study_xml <- readings |>
  filter(as.Date(datetime) >= as.Date("2025-05-28"),
         as.Date(datetime) <= as.Date("2025-06-19"),
         tank %in% c(TANK_28C, TANK_31C)) |>
  select(datetime, probe, value)
missing_compiled <- anti_join(study_xml, compiled_unique,
                              by = c("datetime", "probe", "value" = "compiled_c"))
write_csv(tibble(first_compiled = min(compiled_unique$datetime),
                 last_compiled = max(compiled_unique$datetime),
                 n_unique_compiled = nrow(compiled_unique),
                 n_study_xml = nrow(study_xml),
                 n_study_missing_from_compiled = nrow(missing_compiled)),
          file.path(TBL_DIR, "08_apex_compiled_standalone_check.csv"))
stopifnot(nrow(missing_compiled) == 0)

hourly <- readings |> mutate(datetime = floor_date(datetime, "hour")) |>
  group_by(datetime, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(), .groups = "drop")
daily <- readings |> mutate(date = as.Date(datetime)) |>
  group_by(date, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(),
            n_hours = n_distinct(floor_date(datetime, "hour")), .groups = "drop")
saveRDS(hourly, file.path(DATA_PROC, "apex_temperature.rds"))
saveRDS(daily, file.path(DATA_PROC, "apex_temperature_daily.rds"))

assignments <- read_csv(file.path(recovered, "tank_assignments.csv"),
                       show_col_types = FALSE) |>
  filter(project == "17. LTH expression by temperature") |>
  transmute(date = mdy(date), tank = as.integer(tank),
            treatment = paste0(treatment, "C"), phase = str_trim(.data[["acc.ramp.hold"]]))
stopifnot(!anyDuplicated(select(assignments, date, tank)),
          all(assignments$treatment == tank_treatment(assignments$tank)))
temp_daily <- daily |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+"))) |>
  inner_join(assignments, by = c("date", "tank"), relationship = "one-to-one") |>
  mutate(day = as.integer(date - as.Date("2025-06-04")))
coverage <- expand_grid(date = as.Date("2025-06-04") + 0:16,
                        tank = sort(c(TANK_28C, TANK_31C))) |>
  left_join(daily |> mutate(tank = as.integer(str_extract(probe, "[0-9]+"))),
            by = c("date", "tank")) |>
  mutate(day = as.integer(date - as.Date("2025-06-04")),
         treatment = tank_treatment(tank),
         across(c(n, n_hours), ~replace_na(.x, 0L)))
write_csv(coverage, file.path(TBL_DIR, "08_apex_daily_coverage.csv"))
summary <- temp_daily |> filter(day >= 0, day <= 15) |>
  group_by(treatment) |>
  summarise(n_tanks = n_distinct(tank), n_tank_days = n(),
            mean_c = mean(value_mean), min_daily_c = min(value_mean),
            max_daily_c = max(value_mean), min_hours = min(n_hours), .groups = "drop")
write_csv(summary, file.path(TBL_DIR, "08_apex_treatment_summary.csv"))
stopifnot(all(filter(coverage, day <= 15)$n == 144),
          all(filter(coverage, day <= 15)$n_hours == 24))
clock_sensitivity <- map_dfr(c(0, -3), function(offset) {
  readings |> mutate(date = as.Date(datetime + hours(offset)),
                     treatment = tank_treatment(tank)) |>
    filter(!is.na(treatment), date >= as.Date("2025-06-04"),
           date <= as.Date("2025-06-19")) |>
    group_by(date, tank, treatment) |> summarise(mean_c = mean(value), .groups = "drop") |>
    group_by(treatment) |> summarise(mean_c = mean(mean_c), .groups = "drop") |>
    mutate(clock_shift_hours = offset)
})
write_csv(clock_sensitivity, file.path(TBL_DIR, "08_apex_clock_sensitivity.csv"))

p_apex <- ggplot(temp_daily, aes(day, value_mean, group = tank, colour = treatment)) +
  annotate("rect", xmin = -7, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.45, alpha = 0.85) + geom_point(size = 1.1) +
  scale_colour_manual(values = c('28C' = "#56B4E9", '31C' = "#D55E00"),
                      labels = c('28C' = "28 °C", '31C' = "31 °C"), name = "Treatment") +
  scale_x_continuous(breaks = c(-7, -3, 0, 5, 10, 15)) +
  labs(x = "Day relative to tip clipping (Day 0 = 4 June 2025)",
       y = "Daily mean tank temperature (°C)",
       title = "Tank temperatures through the final sampling day",
       subtitle = "One line per tank; shaded area precedes clipping; dashed lines mark targets") +
  theme_pub(10)
save_fig(p_apex, "08_apex_temperature", width = 180, height = 105)

# Tanks were reused; do not label dates outside the main assignment as treatment.
p_full <- daily |> mutate(tank = as.integer(str_extract(probe, "[0-9]+"))) |>
  filter(tank %in% c(TANK_28C, TANK_31C)) |>
  ggplot(aes(date, value_mean)) +
  geom_line(linewidth = 0.3) + facet_wrap(~tank, ncol = 4) +
  geom_vline(xintercept = as.Date(c("2025-05-28", "2025-06-19")),
             linetype = "dotted", colour = "grey50") +
  labs(x = NULL, y = "Daily mean temperature (°C)",
       title = "Full recovered logger record, by tank",
       subtitle = "Dotted lines bound the main LTH assignment period; other dates include other uses") +
  theme_pub(9)
save_fig(p_full, "08b_apex_temperature_full", width = 200, height = 120)
print(summary)
print(count(coverage, day, n_hours), n = 40)
