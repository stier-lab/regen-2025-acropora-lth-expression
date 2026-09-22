# =============================================================================
# Purpose: PAM Fv/Fm × treatment × wound × day
#          - Mixed model: Fv/Fm ~ treatment * wound * day + thicket + (1|tank) + (1|id)
#          - Genet (thicket) is a FIXED blocking term (only 3 genets — too few to
#            estimate a random-effect variance; same choice as the canonical
#            models in 12_models.R, Bolker 2008 / Gelman 2005)
#          - Within-coral repeated measures via random effect of id
#          - Pre-treatment baseline (day < 0) is shown in the figure but EXCLUDED
#            from the model, which fits a single linear day trend post-wounding
#          - Figure: trajectories with 95% CI, faceted by treatment, color by wound
#
# What & why: Fv/Fm is the dark-adapted maximum photochemical efficiency of
#   photosystem II, measured with a pulse-amplitude-modulation (PAM) fluorometer.
#   It is the standard, non-destructive readout of how stressed the coral's
#   symbiotic algae (Symbiodiniaceae) are: a healthy symbiont sits near ~0.6,
#   and the value drops as heat stress damages the photosynthetic machinery (the
#   first step toward bleaching). Here we ask whether the 31 °C heat treatment
#   depressed Fv/Fm over the experiment, and whether wounding and the genotype of
#   the host coral modified that response. This is the physiological evidence for
#   the headline result — that heat impairs the symbiosis (and thus the energy
#   budget needed for regeneration), even when the wound still closes.
# Input:   data/raw/pam/PAM_data.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/pam_clean.rds
#          figures/02_pam_fvfm_trajectory.{pdf,png}
#          output/tables/02_pam_treatment_contrasts.csv
#          output/models/02_pam_lmer.rds
# =============================================================================

# 00_setup.R loads packages, defines shared paths (DATA_RAW, DATA_PROC, TBL_DIR,
# MOD_DIR), the theme_pub() plot theme, the PAL_WOUND palette, and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# clean_names() converts the spreadsheet headers to predictable snake_case.
pam_raw <- read_csv(file.path(DATA_RAW, "pam", "PAM_data.csv"),
                    show_col_types = FALSE) |>
  janitor::clean_names()

# Coral metadata (genet/thicket assignment etc.) built in an earlier script.
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))

# ---- Clean -----------------------------------------------------------------
# Coerce every column to its type and recover Fv/Fm. The Fv/Fm column was a live
# spreadsheet formula in some rows ("=L2/100"), so it reads back as text; where
# that happens we recompute it from the raw "Y" reading using the PAM convention
# (Y / 1000) rather than trusting the broken cell.
pam <- pam_raw |>
  rename(thicket = matches("^thicket"),
         fv_fm   = matches("^fv_fm")) |>
  mutate(
    date       = as_date(date),
    day        = as.integer(day),                # day relative to wounding (D0)
    # Treatment recorded as the numeric set-point (28/31 °C); make it an ordered
    # 2-level factor with 28C first so model contrasts read as "31C vs 28C".
    treatment  = factor(as.integer(treatment), levels = c(28, 31),
                        labels = c("28C", "31C")),
    tank       = as.integer(tank),
    wound      = factor(wound, levels = c("no", "yes")),  # "no" = baseline level
    thicket    = str_to_lower(str_squish(thicket)),       # genet ID; tidy whitespace/case
    id         = as.integer(id),                 # unique coral fragment ID
    location   = str_to_lower(location),         # probe placement: top vs bottom
    f          = as.numeric(f),                  # raw PAM fluorescence channels
    m          = as.numeric(m),
    y          = as.numeric(y),
    e          = as.numeric(e),
    # PAM convention: Fv/Fm = Y / 1000 when reported as raw "Y". Only fall back to
    # this when the fv_fm cell failed to parse as a number (the broken-formula rows).
    fv_fm      = if_else(is.na(suppressWarnings(as.numeric(fv_fm))),
                         y / 1000,
                         suppressWarnings(as.numeric(fv_fm)))
  ) |>
  # Keep only biologically valid Fv/Fm: it is a ratio bounded in (0, 1), so 0,
  # negatives, and >1 are measurement/parsing errors and are dropped.
  filter(!is.na(fv_fm), fv_fm > 0, fv_fm < 1)

# Average the two within-coral probe readings (top + bottom) into one value per
# coral per day. These are spatial subsamples, not independent fragments or
# interchangeable technical replicates. This average describes the two sampled
# positions; code/42 retains their paired differences for the location question.
pam_avg <- pam |>
  group_by(date, day, treatment, tank, thicket, wound, id) |>
  summarise(fv_fm = mean(fv_fm, na.rm = TRUE), .groups = "drop")

# Save the cleaned, one-row-per-coral-per-day table for downstream scripts.
saveRDS(pam_avg, file.path(DATA_PROC, "pam_clean.rds"))

# ---- Location (top vs bottom) sensitivity check ----------------------------
# Molly's original analysis (code/archive/molly_original/LTH_PAM.R) noted that
# top vs bottom probe placement appeared to differ. The primary pipeline
# averages the two as technical replicates (above). Here we test that decision:
# fit the model on the UN-averaged data with location as a fixed factor and
# report the location terms. If location and its interactions are
# non-significant, that does not establish equivalence. This legacy screen is
# superseded for interpretation by the tank-aware paired analysis in code/42.
if ("location" %in% names(pam) &&
    length(unique(na.omit(pam$location))) > 1) {
  pam_loc <- pam |> mutate(location = factor(location))
  # lmerTest::lmer (not lme4::lmer) so anova() returns Satterthwaite p-values.
  # REML = FALSE (ML) here because this is a test about FIXED effects (the
  # location terms): likelihood-ratio / ML F-tests on fixed effects are only
  # valid under ML, whereas REML is reserved for the final variance estimates.
  # check.conv.singular = "ignore": this saturated 4-way model can drive a
  # variance component to 0; we tolerate that since we need the fixed-effect
  # tests, not the random-effect estimates.
  mod_loc <- lmerTest::lmer(
    fv_fm ~ treatment * wound * day * location + (1 | tank) + (1 | id),
    data = pam_loc, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
  )
  # Pull the ANOVA table and keep only rows that mention "location" — the terms
  # that decide whether probe placement matters.
  loc_anova <- as.data.frame(anova(mod_loc)) |>
    tibble::rownames_to_column("term") |>
    filter(grepl("location", term))
  write_csv(loc_anova, file.path(TBL_DIR, "02b_pam_location_sensitivity.csv"))
  cat("\n=== PAM location (top/bottom) sensitivity — terms involving location ===\n")
  print(loc_anova)
  # The decision that matters is whether location INTERACTS with the
  # experimental factors. A main-effect offset is averaged out; an interaction
  # would mean averaging distorts the treatment/wound/day effects.
  pcol <- intersect(c("Pr(>F)"), names(loc_anova))
  interaction_terms <- loc_anova[grepl(":", loc_anova$term), , drop = FALSE]
  inter_sig <- length(pcol) == 1 &&
    any(interaction_terms[[pcol]] < 0.05, na.rm = TRUE)
  main_sig <- length(pcol) == 1 &&
    any(loc_anova[loc_anova$term == "location", pcol] < 0.05, na.rm = TRUE)
  cat(sprintf(
    "  location main effect %s; location x (experimental factor) interactions %s.\n",
    if (main_sig) "SIGNIFICANT (real top/bottom offset)" else "n.s.",
    if (inter_sig) "SIGNIFICANT — averaging may distort effects; see 02b table"
    else "n.s. — not evidence of equivalence; use the paired analysis in code/42"))
}

# ---- Mixed model -----------------------------------------------------------
# The primary model. Day is continuous (a linear time trend), treatment and
# wound are factors, and they are fully crossed so we can read off whether heat
# and wounding interact and whether their effects change through time.
# Genet (thicket) enters as a FIXED blocking term: with only 3 genets a random
# "genet variance" is near-unidentified, so we estimate the genet means directly
# (matches 12_models.R; Bolker et al. 2008, Gelman 2005).
# Random intercepts account for the remaining non-independence in the design:
#   (1|tank)    — corals sharing a tank share its micro-environment
#   (1|id)      — repeated measures on the same fragment across days
# Pre-treatment baseline (day < 0) is dropped here so the single linear day term
# describes the POST-wounding trajectory only; pooling the baseline into one
# slope would dilute/bias the treatment × day heat signal. The figure below
# still plots the baseline point for context.
# REML = TRUE here (unlike the ML sensitivity model): this is the FINAL model, so
# we want the less-biased REML estimates of the variance components.
pam_mod <- pam_avg |> filter(day >= 0)
mod <- lme4::lmer(
  fv_fm ~ treatment * wound * day + thicket + (1 | tank) + (1 | id),
  data = pam_mod, REML = TRUE
)
saveRDS(mod, file.path(MOD_DIR, "02_pam_lmer.rds"))

# Estimated marginal means: model-predicted Fv/Fm for every treatment×wound cell
# at each observed (post-wounding) day, then pairwise contrasts. adjust = "tukey"
# controls the family-wise error rate across the pairwise comparisons within each
# day. Grid is the measured days, not a hardcoded set, so no day is interpolated.
emm <- emmeans::emmeans(mod, ~ treatment * wound | day,
                        at = list(day = sort(unique(pam_mod$day))))
contrasts_tbl <- as_tibble(pairs(emm, adjust = "tukey"))
write_csv(contrasts_tbl, file.path(TBL_DIR, "02_pam_treatment_contrasts.csv"))

# ---- Figure ----------------------------------------------------------------
# Collapse to one mean (± standard error) per day × treatment × wound cell for
# plotting. These are the raw cell means, not the model-adjusted means.
plot_df <- pam_avg |>
  group_by(day, treatment, wound) |>
  summarise(
    mean = mean(fv_fm, na.rm = TRUE),
    se   = sd(fv_fm, na.rm = TRUE) / sqrt(n()),  # SE = SD / sqrt(n)
    n    = n(),
    .groups = "drop"
  )

# One panel per temperature; within each, a line per wound status with a ±1 SE
# ribbon. The dotted vertical line at day 0 marks the wounding event.
p_pam <- ggplot(plot_df, aes(day, mean,
                             colour = wound, fill = wound, group = wound)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  scale_colour_manual(values = PAL_WOUND,
                      name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND,
                    guide = "none") +
  labs(x = "Day relative to wounding (D0)",
       y = expression(F[v]/F[m]),
       title = "Photochemical efficiency",
       subtitle = "Mean ± 1 SE across genets and tanks") +
  theme_pub(10)

save_fig(p_pam, "02_pam_fvfm_trajectory", width = 170, height = 90)

# Console summary
cat("\n=== PAM mixed model ===\n")
print(summary(mod)$coefficients)
cat("\nWrote: pam_clean.rds, 02_pam_lmer.rds, 02_pam_treatment_contrasts.csv,",
    "02_pam_fvfm_trajectory.{pdf,png}\n")
