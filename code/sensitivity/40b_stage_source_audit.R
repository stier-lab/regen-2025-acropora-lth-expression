# Verify disputed scores against archived and live-source snapshots by ID/day.
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
keys <- c("id", "day")
traits <- c("tip_exist", "tip_extension", "new_corallites_on_tip")
cols <- c(keys, traits)
raw <- read_csv("data/raw/physio_morphology/data.csv", show_col_types = FALSE)
book <- read_excel("data/raw/physio_morphology/physio_characterization_log_-_complete.xlsx",
                   sheet = "data")
live <- read_csv("data/external/molly_stage_review_live_2026-09-15.csv",
                 show_col_types = FALSE)
clean <- readRDS("data/processed/physio_clean.rds")
normalize <- function(x) x |> filter(id %in% c(115, 121)) |>
  select(all_of(cols)) |> arrange(id, day)
a <- normalize(raw)
b <- normalize(book)
c <- normalize(live)
stopifnot(nrow(a) == 32, !anyDuplicated(a[keys]),
          identical(as.data.frame(a), as.data.frame(b)),
          identical(as.data.frame(a), as.data.frame(c)),
          all(as.matrix(a[traits]) %in% c("yes", "no")))
encoded <- a |> mutate(across(all_of(traits), ~as.integer(.x == "yes")))
stopifnot(isTRUE(all.equal(as.data.frame(encoded),
                          as.data.frame(normalize(clean)), check.attributes = FALSE)))
write_csv(tibble(check = c("CSV vs workbook", "CSV vs live sheet", "CSV vs processed"),
                 n_rows = 32L, n_trait_values = 96L, agrees = TRUE),
          "output/tables/40b_stage_source_agreement.csv")
events <- read_csv("output/tables/36_molly_stage_timing.csv", show_col_types = FALSE) |>
  filter(stage_order == 6)
# Descriptive sensitivity only: omit disputed fragments, without recoding them.
counts <- bind_rows(
  events |> mutate(scenario = "Original observations"),
  events |> filter(!id %in% c(115,121)) |> mutate(scenario = "Omit both flagged fragments")
) |> group_by(scenario, treatment) |>
  summarise(n_total = n(), n_reached = sum(event), .groups = "drop")
write_csv(counts, "output/tables/40b_stage_exclusion_counts.csv")
print(counts)
