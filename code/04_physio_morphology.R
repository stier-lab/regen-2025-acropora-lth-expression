# =============================================================================
# Purpose: 9 binary morphological characterizations of wound healing over time.
#          - For wounded corals only, plot cumulative proportion expressing each
#            trait by day × treatment, with genet (thicket) as a third
#            dimension (linetype) so genotype-specific responses are visible.
#          - Logistic mixed models per trait with treatment × thicket
#            interaction so genet-level differences in heat response are
#            tested directly (rather than absorbed in a random intercept).
#
# What & why: each wounded coral was photographed/scored over time for nine
#   yes/no signs of how the cut surface was healing. These split into two
#   biological phases: WOUND CLOSURE (the tissue sealing the lesion — e.g.
#   "wound smoothed", "pigment over wound", "hole in center" disappearing) and
#   REGENERATION (rebuilding lost skeletal structure — e.g. "tip exists", "tip
#   extension", "new corallites on tip"). This distinction underlies the
#   paper's headline result: heat (31 °C) blocks REGENERATION but not wound
#   CLOSURE. Because the three genets (A/C/D) differ in heat resilience
#   (C > D > A), genet is modelled as a FIXED treatment × genet interaction —
#   we WANT to estimate each genotype's heat response, not average it away into
#   a random intercept. Each trait is a separate binary outcome, so each gets
#   its own logistic mixed model.
# Input:   data/raw/physio_morphology/data.csv
# Output:  data/processed/physio_clean.rds
#          figures/04_morphology_trajectories.{pdf,png}        — pooled-by-genet
#          figures/04b_morphology_trajectories_by_genet.{pdf,png} — split by genet
#          output/tables/04_morphology_trait_glmm_summaries.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_RAW, DATA_PROC, TBL_DIR), the
# theme_pub() theme, and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# guess_max = 2000 makes readr scan many rows before guessing column types
# (some trait columns are blank early on, which would otherwise mis-type them).
# suppressWarnings hides the expected parse messages from that coercion.
# clean_names() -> snake_case; drop rows missing the essentials (species, id).
raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "physio_morphology", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))

parse_lth_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(lubridate::as_date(x))
  lubridate::as_date(lubridate::parse_date_time(
    as.character(x),
    orders = c("mdy HMS", "mdy", "ymd HMS", "ymd"),
    quiet = TRUE
  ))
}

# ---- Clean -----------------------------------------------------------------
# Type-coerce the design columns and turn the nine yes/no trait columns into
# 0/1 integers so they can be modelled and averaged as proportions.
ph <- raw |>
  # Raw header "wounded" -> project-wide "wound"; match thicket by prefix.
  rename(thicket = matches("^thicket"),
         wound   = wounded) |>
  mutate(
    date      = parse_lth_date(date),
    day       = as.integer(day),                 # day relative to wounding (D0)
    # Treatment as set-point (28/31 °C); 28C is the reference level.
    treatment = factor(as.integer(treatment), levels = c(28, 31),
                       labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),  # genet ID; tidy case/whitespace
    id        = as.integer(id),                  # unique coral fragment ID
    wound     = factor(wound, levels = c("no", "yes")),  # "no" = reference level
    # Encode all nine trait columns at once: lower-case, trim, then test == "yes"
    # to get a clean 0/1 integer. Anything that isn't exactly "yes" -> 0 (NA
    # inputs stay NA via the comparison and are filtered later).
    across(c(polyps_out, hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension,
             new_corallites_on_tip, algae_on_wound),
           ~ as.integer(str_to_lower(str_squish(.x)) == "yes"))
  )

# ---- Combine hole_in_center + polyp_in_hole -> axial_polyp_formation --------
# In the raw spreadsheet these two columns are byte-identical (same 0/1 values
# AND the same NA pattern in every row). M. Brzezinski (pers. comm., 2026)
# confirmed this is not a data-entry slip but how the trait was scored: the
# central "hole" IS the axial polyp hole, which forms around the regenerating
# axial polyp, so the two co-occur (a hole without a polyp was seen once or
# twice). They are therefore one observable, combined here into a single trait,
# `axial_polyp_formation` (the axial corallite/calyx + polyp). We keep the two
# original columns in the saved table for provenance and assert they still match
# (if a future re-score makes them differ, this stops and prompts a rethink).
stopifnot(identical(ph$hole_in_center, ph$polyp_in_hole))
ph <- ph |> mutate(axial_polyp_formation = hole_in_center)

# Save the cleaned, one-row-per-observation table for downstream scripts (the two
# original columns retained alongside the combined axial_polyp_formation).
saveRDS(ph, file.path(DATA_PROC, "physio_clean.rds"))

# ---- Long-form for plotting ------------------------------------------------
# Trait columns paired with their display labels (single source of truth for both
# the modelled trait order and the facet labels). hole_in_center/polyp_in_hole are
# represented by the combined trait axial_polyp_formation (see above).
trait_labels <- c(
  polyps_out            = "Polyps out",
  axial_polyp_formation = "Axial polyp formation",
  wound_smoothed        = "Wound smoothed",
  pigment_over_wound    = "Pigment over wound",
  tip_exist             = "Tip exists",
  tip_extension         = "Tip extension",
  new_corallites_on_tip = "New corallites on tip",
  algae_on_wound        = "Algae on wound"
)
traits <- names(trait_labels)

# Restrict to WOUNDED corals (unwounded controls have no wound to heal), then
# stack the nine traits into one long trait/expressed column for faceting.
long <- ph |>
  filter(wound == "yes") |>
  select(day, treatment, tank, thicket, id, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  filter(!is.na(expressed))                      # drop unscored coral-day-traits

# Proportion of wounded corals expressing each trait, by day × treatment. Because
# "expressed" is 0/1, its mean IS the proportion showing the trait.
prop_df <- long |>
  group_by(day, treatment, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits, labels = unname(trait_labels)))

# Sample sizes for the figure subtitle, computed from the data (not hardcoded):
# wounded corals per treatment and per genet × treatment cell.
n_wound_per_trt  <- long |> distinct(id, treatment) |> count(treatment) |>
  pull(n) |> (\(x) if (length(unique(x)) == 1) unique(x) else paste(range(x), collapse = "–"))()
n_wound_per_cell <- long |> distinct(id, treatment, thicket) |>
  count(treatment, thicket) |>
  pull(n) |> (\(x) if (length(unique(x)) == 1) unique(x) else paste(range(x), collapse = "–"))()

# ---- Pooled trajectory figure (one panel per trait) ------------------------
# Lines coloured by temperature; one panel per healing trait. Blue = 28 °C,
# orange/red = 31 °C. Dotted line at day 0 = wounding.
p_traits <- ggplot(prop_df, aes(day, prop,
                                colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +          # y is a proportion, fixed 0-100%
  labs(x = "Day relative to wounding (D0)",
       y = "Wounded corals expressing trait",
       title = "Morphological wound-healing characteristics",
       subtitle = sprintf("Wounded corals only (n = %s per treatment)",
                          n_wound_per_trt)) +
  theme_pub(9)

save_fig(p_traits, "04_morphology_trajectories", width = 200, height = 130)

# ---- Per-genet trajectory plot --------------------------------------------
# Same proportions, but now also split by genet (thicket) so genotype-specific
# heat responses (the C > D > A resilience ordering) are visible by eye.
prop_genet_df <- long |>
  group_by(day, treatment, thicket, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits, labels = unname(trait_labels)))

# Temperature -> colour, genet -> linetype; group by their interaction so each
# temperature×genet combination is drawn as its own line.
p_traits_genet <- ggplot(prop_genet_df,
                          aes(day, prop, colour = treatment,
                              linetype = thicket,
                              group = interaction(treatment, thicket))) +
  geom_line(linewidth = 0.55, alpha = 0.85) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  # Genets A/C/D get distinct dash patterns ("22"/"44" are on/off dash specs).
  scale_linetype_manual(values = c(a = "solid", c = "22", d = "44"),
                        name = "Genet") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Day relative to wounding (D0)",
       y = "Wounded corals expressing trait",
       title = "Morphological wound-healing trajectories by genet",
       subtitle = sprintf("Wounded corals only (n = %s per genet × treatment cell)",
                          n_wound_per_cell)) +
  theme_pub(9)

save_fig(p_traits_genet, "04b_morphology_trajectories_by_genet",
         width = 210, height = 135)

# ---- Per-trait GLMM with genet × treatment ---------------------------------
# Question: does the temperature effect on trait expression depend on genet?
# Each trait is a binary outcome, so we fit one logistic mixed model per trait.
# Fixed effects: treatment × day × thicket (genet is FIXED so we can read off
# each genotype's heat response, not absorb it into a random intercept).
# Random intercepts (1|tank) + (1|id) handle shared-tank and repeated-measures
# non-independence.
fit_one <- function(trait_name) {
  d <- long |>
    filter(trait == trait_name, day >= 0) |>    # post-wounding days only (healing window)
    mutate(thicket = factor(thicket))
  # Skip traits that are all-0 or all-1: a binary model needs both outcomes to
  # estimate anything (avoids complete separation).
  if (length(unique(d$expressed)) < 2) return(NULL)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      # bobyqa optimizer + a high iteration cap (2e5) help these sparse binary
      # models converge; the whole fit is wrapped in tryCatch so one trait that
      # fails to converge returns NULL instead of stopping the loop.
      lme4::glmer(expressed ~ treatment * day * thicket + (1 | tank) + (1 | id),
                  family = binomial, data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5)))
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  # broom.mixed::tidy pulls the fixed-effect coefficient table into a tidy frame.
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed") |>
    mutate(trait = trait_name, n = nrow(d))
  # Type-II Wald ANOVA — reports significance of treatment, day, thicket and
  # their interactions. Type II respects marginality (tests each main effect
  # after the others but before its higher-order interactions).
  av <- as.data.frame(car::Anova(fit, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(trait = trait_name)
  list(tidy = tidy_fit, anova = av)
}
# Fit every trait, name the list, and drop the NULLs (skipped/failed traits).
trait_results <- map(traits, fit_one) |> setNames(traits) |> compact()
# Stack the per-trait coefficient and ANOVA tables into two long data frames.
trait_tidy <- map_dfr(trait_results, "tidy")
trait_anova <- map_dfr(trait_results, "anova")
write_csv(trait_tidy,  file.path(TBL_DIR, "04_morphology_trait_glmm_summaries.csv"))
write_csv(trait_anova, file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"))

cat("Wrote physio_clean.rds, 04_morphology_trajectories.{pdf,png},",
    "04b_morphology_trajectories_by_genet.{pdf,png},",
    "04_morphology_trait_glmm_summaries.csv,",
    "04_morphology_trait_anova_genet.csv\n")
