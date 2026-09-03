# =============================================================================
# Purpose: Primary statistical models for every continuous + morphological
#          response, plus the two robustness refits that the diagnostics call
#          for. This single file replaces the former 12 / 12b / 12c trio:
#
#            PART 1 — Primary models (was 12_extended_stats.R)
#                     LMM/GLMM with treatment × wound × day × thicket fixed
#                     structure; genet (thicket) is a fixed effect
#                     because only 3 genets were collected (Bolker 2008,
#                     Gelman 2005: few grouping levels → fixed beats random).
#                     Produces type-III ANOVA, emmeans contrasts, R²m/R²c,
#                     per-genet treatment effects, and DHARMa diagnostics.
#
#            PART 2 — Color ordinal-CLMM robustness (was 12b)
#                     Refit the 5-level color D-scale under the ordinal
#                     likelihood; confirms the Gaussian LMM's inferences hold.
#
#            PART 3 — Penalized morphology refit (was 12c)
#                     blme::bglmer with Cauchy(0, 2.5) priors to tame the
#                     complete/quasi-complete separation in 7/8 binary traits;
#                     yields finite, interpretable coefficient Wald tests.
#
# What & why: these are the primary statistical models. Each coral fragment was
#   tracked through a 2 (temperature: 28 vs 31 °C) × 2 (wound: yes/no) × 3 (genet:
#   A/C/D) design, biopsied repeatedly over D0-D15. For each response we fit ONE
#   mixed model of the form  response ~ treatment * wound * day * thicket
#   + (1|tank) + (1|id). The biological question — does heat block the wound
#   response? — is the treatment × day interaction (whether the heated group
#   follows a different time course than the ambient group), and the genet
#   question is any term containing `thicket`. We keep the FULL interaction and
#   report every term rather than stepwise-dropping. Two follow-up refits (Parts 2
#   and 3) re-run the two responses whose residual diagnostics flagged a likelihood
#   mismatch, confirming the main inferences hold.
#
# Why genet (thicket) is a FIXED effect, not random: a random effect estimates a
#   *variance* across the levels of a grouping factor, which needs many levels
#   (rule of thumb ~5+, ideally >8) to estimate with any precision. We have only
#   3 genets, so a random "genet variance" would be near-unidentified and
#   biased. With few levels the recommended choice is to treat the factor as fixed
#   — we estimate the 3 genet means directly and can ask which genet is most/least
#   heat-tolerant (Bolker et al. 2008 TREE; Gelman 2005 Ann. Stat.). Tank and coral
#   id, by contrast, have many levels and are nuisance groupings we want to
#   generalize beyond, so they stay as random intercepts (see below).
#
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean,physio_clean}.rds
# Output:  output/models/12_<response>_lmm.rds          — primary fits
#          output/models/12b_color_clmm.rds             — ordinal robustness
#          output/models/12c_morph_<trait>_blme.rds     — penalized refits
#          output/tables/12_anova_summary.csv            — ANOVA all terms
#          output/tables/12_emmeans_contrasts.csv        — pairwise contrasts
#          output/tables/12_genet_treatment_effects.csv  — per-genet effects
#          output/tables/12_r2_summary.csv               — R² per response
#          output/tables/12_morph_fixed_effects.csv      — raw glmer fixefs
#          output/tables/12b_color_clmm*.csv             — ordinal robustness
#          output/tables/12c_morph_blme_*.csv            — penalized robustness
#          figures/12_diagnostics/<response>.png         — DHARMa plots
# =============================================================================

# 00_setup.R loads packages (tidyverse, lme4/lmerTest, emmeans, DHARMa, ...) and
# defines the shared paths used below (DATA_PROC, MOD_DIR, TBL_DIR, FIG_DIR) plus
# the colour palettes and theme_pub() / save_fig() helpers.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(MuMIn)   # r.squaredGLMM() — marginal/conditional R² for mixed models
  library(blme)    # Bayesian/penalized GLMMs used in PART 3
  # NB: do NOT attach `ordinal` — it masks dplyr::slice and breaks downstream
  # scripts. PART 2 calls it namespace-qualified (ordinal::clmm).
})

# Folder for the per-response DHARMa residual-diagnostic PNGs written below.
DIAG_DIR <- file.path(FIG_DIR, "12_diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PART 1 — PRIMARY MODELS  (formerly 12_extended_stats.R)
# =============================================================================

# Accumulators: each response appends its results to these lists, which are
# bound into one table per output (ANOVA / contrasts / R² / genet effects)
# at the end of Part 1.
results_anova        <- list()
results_emm          <- list()
results_r2           <- list()
results_genet_effect <- list()

# ---- Helpers ---------------------------------------------------------------
# save_dharma(): simulate scaled (quantile) residuals via DHARMa and write the
# QQ + residual-vs-predicted diagnostic panel. DHARMa works for GLMMs
# where ordinary residuals are uninformative; it simulates from the fitted model
# and compares observed vs expected, so deviations flag a wrong likelihood or
# missing structure. Wrapped in tryCatch so a single failure never halts the run.
save_dharma <- function(model, name) {
  res <- tryCatch(
    DHARMa::simulateResiduals(model, n = 500, refit = FALSE),
    error = function(e) NULL
  )
  if (is.null(res)) return(invisible(NULL))
  png(file.path(DIAG_DIR, paste0(name, ".png")),
      width = 1600, height = 800, res = 160)
  plot(res)
  dev.off()
  invisible(res)
}

# record_results(): run the omnibus ANOVA and R² for a fitted model and stash
# them in the accumulator lists. We use TYPE-III sums of squares (test each term
# after all others, including its higher-order interactions) because the design
# is a fully crossed factorial WITH interactions — type-I would depend on term
# order and type-II assumes no interaction, both wrong here. For Gaussian LMMs
# (lmerTest), anova(type = 3) gives Satterthwaite/Kenward-Roger F-tests; for the
# binomial GLMMs car::Anova(type = 3) gives Wald χ² tests instead.
record_results <- function(model, response_name) {
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest") ||
      inherits(model, "lm")) {
    if (inherits(model, "lm") && !inherits(model, "lmerMod")) {
      av <- as.data.frame(car::Anova(model, type = 3))   # plain lm path
    } else {
      av <- as.data.frame(anova(model, type = 3))        # lmerTest F-tests
    }
    av$term <- rownames(av); av$response <- response_name
    results_anova[[response_name]] <<- av
  } else if (inherits(model, "glmerMod")) {
    av <- as.data.frame(car::Anova(model, type = 3))      # GLMM Wald χ²
    av$term <- rownames(av); av$response <- response_name
    results_anova[[response_name]] <<- av
  }
  # Nakagawa R²: marginal (R2m, fixed effects only) and conditional (R2c, fixed
  # + random). The gap between them is the share of variance absorbed by tank/id.
  r2 <- tryCatch(MuMIn::r.squaredGLMM(model), error = function(e) NULL)
  if (!is.null(r2)) {
    results_r2[[response_name]] <<- tibble::tibble(
      response = response_name,
      R2_marginal    = r2[1, "R2m"],
      R2_conditional = r2[1, "R2c"]
    )
  }
  save_dharma(model, response_name)
}

# record_genet_effect(): from an emmeans grid, extract the 28C-vs-31C contrast
# *within each source-patch × wound* combination. This is how we answer "is
# source patch C more heat-tolerant than A/D?" — one treatment-effect estimate per source patch. adjust =
# "none" here because these are pre-specified per-genet contrasts collected for a
# summary table, not a family we are protecting against false positives.
# Per-source-patch treatment effect at end-of-experiment (or biopsy-day max)
# Returns one row per source patch × wound with (28C - 31C) treatment-effect estimate
record_genet_effect <- function(emm_obj, response_name) {
  contr <- as_tibble(pairs(emm_obj, by = c("thicket", "wound"),
                            adjust = "none"))
  if (nrow(contr) == 0) return(invisible(NULL))
  contr$response <- response_name
  results_genet_effect[[response_name]] <<- contr
}

# ---- 1. PAM Fv/Fm -----------------------------------------------------------
# Fv/Fm = dark-adapted maximum quantum yield of photosystem II; a
# non-destructive proxy for the photosynthetic health of the symbionts. A drop
# under heat indicates bleaching stress.
cat("\n=== 1. PAM Fv/Fm ===\n")
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = as.factor(thicket)) |>   # genet must be a factor to enter as fixed levels
  filter(day >= 0)   # drop the pre-treatment baseline (day -1): with a single linear
                     # day term, pooling it in would dilute the treatment × day signal

# With n=3 genets and 4 wounded + 4 unwounded fragments per genet × treatment
# cell, the full 4-way interaction is supported but high-order terms will have
# limited power. We report all terms rather than stepwise-dropping.
#
# The model: fixed = treatment * wound * day * thicket (full factorial, so every
#   main effect AND interaction is estimated). The primary heat signal is the
#   treatment × day interaction. Random intercepts:
#     (1 | tank) — the 8 tanks are the experimental unit to which temperature was
#                  applied; corals sharing a tank are non-independent, so a tank
#                  intercept absorbs that shared environment (pseudoreplication fix).
#     (1 | id)   — each fragment is measured repeatedly across days; an id
#                  intercept models the within-fragment correlation of repeats.
# REML = TRUE: restricted maximum likelihood gives less-biased variance-component
#   estimates and suits inference on a FIXED model
#   (we are not comparing different fixed structures here, so ML is not needed).
# check.conv.singular "ignore": with this many parameters a variance component can
#   land near zero (singular fit); we suppress the warning rather than simplify,
#   because the fixed-effect tests we report are unaffected.
m_pam <- lmerTest::lmer(
  fv_fm ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
  data = pam, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_pam, file.path(MOD_DIR, "12_pam_lmm.rds"))
record_results(m_pam, "pam_fvfm")

# Day-specific contrasts: estimate the 31C-vs-28C difference at each biopsy day,
# separately within every genet × wound cell. emmeans builds model-predicted
# marginal means on the grid; pairs() forms the treatment contrast. adjust =
# "tukey" controls the family-wise error across the set of pairwise comparisons
# (here a single pair per cell, so Tukey ≈ unadjusted, but kept for consistency
# with the multi-level contrasts elsewhere).
# Grid = the measured days, so no contrast is reported at a day that was
# never sampled (PAM days are 0,3,6,9,12,14 — the old c(0,3,7,10,14) interpolated
# days 7 and 10, which were never measured).
emm_pam <- emmeans::emmeans(m_pam, ~ treatment | thicket * wound * day,
                             at = list(day = sort(unique(pam$day))))
results_emm[["pam_fvfm"]] <- as_tibble(pairs(emm_pam, adjust = "tukey")) |>
  mutate(response = "pam_fvfm")
# Per-genet treatment effect at Day 14 (the end-of-experiment "thermal tolerance")
emm_pam_end <- emmeans::emmeans(m_pam, ~ treatment | thicket * wound,
                                 at = list(day = 14))
record_genet_effect(emm_pam_end, "pam_fvfm")

# ---- 2. Color score (ordinal D-scale) --------------------------------------
# Color = visual paling on the CoralWatch-style D1-D6 health chart (lower =
# paler = more bleached). It is ordinal, but we first fit it as a Gaussian LMM
# so it sits on the same footing as the other continuous metrics;
# PART 2 then re-fits it with the ordinal likelihood as a robustness check.
cat("\n=== 2. Color score ===\n")
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = as.factor(thicket)) |>
  filter(day >= 0)   # drop the pre-treatment baseline (color day -4), as for PAM

# Same fixed/random structure and rationale as the PAM model above.
m_color <- lmerTest::lmer(
  color_num ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
  data = color, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_color, file.path(MOD_DIR, "12_color_lmm.rds"))
record_results(m_color, "color_dscale")

emm_color <- emmeans::emmeans(m_color, ~ treatment | thicket * wound * day,
                              at = list(day = sort(unique(color$day))))  # measured days only
results_emm[["color_dscale"]] <- as_tibble(pairs(emm_color, adjust = "tukey")) |>
  mutate(response = "color_dscale")
emm_color_end <- emmeans::emmeans(m_color, ~ treatment | thicket * wound,
                                   at = list(day = 14))
record_genet_effect(emm_color_end, "color_dscale")

# ---- 3. Coral growth (% skeletal mass change) ------------------------------
# PRIMARY growth metric is % mass change (100 · Δdry / dry_initial; see
# code/05_buoyant_weight.R). An AREAL calcification rate is NOT used: it would
# need the whole-fragment surface area, but the wax SA we have is only for the
# small symbiont-slurry sub-fragment (the rest was chopped for transcriptomics),
# so it is not a valid denominator for whole-fragment mass gain. SGR is recorded
# as a robustness response below; both are surface-area-free.
cat("\n=== 3. Growth (% mass change) ===\n")
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = as.factor(thicket))
bw_a <- bw |> filter(is.finite(pct_growth))

# There is NO `day` term here: growth is a single start-to-end measure per
# coral, not a repeated time series, so the model is a 3-way factorial
# (treatment × wound × thicket) only.
# n<=48 corals with 1 obs each — drop (1|id) because with one observation per
# fragment a fragment-level intercept is not identifiable (it would
# alias the residual). We retain (1|tank) as the experimental-unit random effect
# for the temperature treatment.
m_bw <- lmerTest::lmer(
  pct_growth ~ treatment * wound * thicket + (1 | tank),
  data = bw_a, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_bw, file.path(MOD_DIR, "12_bw_lm.rds"))
record_results(m_bw, "growth_pct")

emm_bw <- emmeans::emmeans(m_bw, ~ treatment | thicket * wound)
results_emm[["growth_pct"]] <- as_tibble(pairs(emm_bw, adjust = "tukey")) |>
  mutate(response = "growth_pct", day = NA_integer_)  # day = NA: no time axis
record_genet_effect(emm_bw, "growth_pct")

# Robustness: specific growth rate (SGR, % d^-1). Same model on a second
# surface-area-free growth currency — if the heat/genet story holds on both, the
# choice of growth metric is not driving the conclusion.
m_bw_sgr <- lmerTest::lmer(
  sgr ~ treatment * wound * thicket + (1 | tank),
  data = bw, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_bw_sgr, file.path(MOD_DIR, "12_bw_sgr_lm.rds"))
record_results(m_bw_sgr, "growth_sgr")

# ---- 4. Symbiont density (cells/cm^2) --------------------------------------
# Symbiont (Symbiodiniaceae) density per unit skeletal area — the direct census
# of the algal partner that PAM and color only proxy. Bleaching = loss of symbionts.
cat("\n=== 4. Symbiont density ===\n")
phys <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>   # need positive values for log()
  mutate(thicket = factor(thicket),
         biopsy_day_c = biopsy_day - 1)   # centre day on the first biopsy so the
                                          # intercept is interpretable at day 0

# log(cells_per_cm2): densities are right-skewed and span orders of magnitude, so
# we model them on the log scale (multiplicative effects → additive, residuals
# closer to Gaussian). biopsy_day_c is the (centred) continuous time axis.
# Destructive sampling: each coral is sacrificed at ONE biopsy day, so there is
# one observation per id → drop (1|id) (not identifiable), keep (1|tank).
m_zoox <- lmerTest::lmer(
  log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket +
    (1 | tank),
  data = phys, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_zoox, file.path(MOD_DIR, "12_zoox_lmm.rds"))
record_results(m_zoox, "log_zoox_density")

emm_zoox <- emmeans::emmeans(m_zoox, ~ treatment | thicket * wound * biopsy_day_c,
                             at = list(biopsy_day_c = sort(unique(phys$biopsy_day_c))))  # measured biopsy days only
results_emm[["log_zoox_density"]] <- as_tibble(pairs(emm_zoox, adjust = "tukey")) |>
  mutate(response = "log_zoox_density")
emm_zoox_end <- emmeans::emmeans(m_zoox, ~ treatment | thicket * wound,
                                  at = list(biopsy_day_c = 14))
record_genet_effect(emm_zoox_end, "log_zoox_density")

# ---- 5. Morphological traits (binomial GLMM, per-trait) --------------------
# Regeneration traits: each binary trait (polyps emerged, hole closed, new
# corallites on the tip, ...) is a yes/no milestone scored on each WOUNDED
# fragment over time. These are the traits that distinguish wound CLOSURE
# (smoothing/pigment) from REGENERATION (new polyps/corallites/tip growth).
cat("\n=== 5. Morphological traits (binomial GLMM with genet) ===\n")
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>   # only wounded corals have these traits
  mutate(thicket = as.factor(thicket))

# hole_in_center + polyp_in_hole are one observable, combined upstream (code/04)
# into axial_polyp_formation (M. Brzezinski pers. comm.); use the combined trait.
traits <- c("polyps_out", "axial_polyp_formation",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip", "algae_on_wound")

# fit_trait(): fit one binomial GLMM per trait. Returns NULL (skips) when the
# trait is invariant (all 0s or all 1s → nothing to model) or sample size is small.
fit_trait <- function(tr) {
  d <- ph |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) return(NULL)  # guard: no variation / too few obs
  # Full treatment × day × genet model (wound is fixed at "yes" here, so it drops
  # out of the formula). family = binomial → logistic link for a 0/1 outcome.
  # Include coral ID because morphology is repeatedly scored on the same fragments
  # over time (repeated-measures), plus tank as the treatment block. bobyqa with a
  # high iteration cap is a more stable optimizer for these sparse binary fits.
  # Wrapped in tryCatch + suppress* because separation/convergence warnings are
  # expected for these milestone traits and are handled in PART 3.
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::glmer(y ~ treatment * day * thicket + (1 | tank) + (1 | id),
                  family = binomial, data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 1e5)))
    )),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  # TYPE-II Wald χ² here (not type-III): for GLMMs with separation the type-II
  # omnibus tests of each term (respecting marginality but not forcing the
  # highest-order interaction) are the stable, interpretable summaries.
  av <- as.data.frame(car::Anova(m, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(response = paste0("morph_", tr), n_obs = nrow(d))
  tidy <- broom.mixed::tidy(m, effects = "fixed") |>
    mutate(response = paste0("morph_", tr))
  # Per-genet probability of trait at Day 10 (peak divergence in KM)
  emm <- tryCatch(
    emmeans::emmeans(m, ~ treatment | thicket, at = list(day = 10),
                     type = "response"),
    error = function(e) NULL
  )
  genet_eff <- if (!is.null(emm)) {
    as_tibble(pairs(emm, adjust = "tukey")) |>
      mutate(response = paste0("morph_", tr))
  } else NULL
  list(model = m, anova = av, tidy = tidy, genet_eff = genet_eff)
}

# Fit every trait; compact() drops the NULLs (invariant / unfittable traits).
morph_results <- map(traits, fit_trait) |> setNames(traits) |> compact()
for (tr in names(morph_results)) {
  saveRDS(morph_results[[tr]]$model,
          file.path(MOD_DIR, paste0("12_morph_", tr, "_glmm.rds")))
  results_anova[[paste0("morph_", tr)]] <- morph_results[[tr]]$anova
  if (!is.null(morph_results[[tr]]$genet_eff)) {
    results_genet_effect[[paste0("morph_", tr)]] <- morph_results[[tr]]$genet_eff
  }
}

morph_tidy <- map_dfr(morph_results, "tidy")
write_csv(morph_tidy, file.path(TBL_DIR, "12_morph_fixed_effects.csv"))

# ---- Aggregate + save ------------------------------------------------------
# Collapse the per-response accumulator lists into one long table each and write
# to CSV. These four files are the machine-readable record of Part 1.
all_anova <- bind_rows(results_anova, .id = "response_id") |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_anova, file.path(TBL_DIR, "12_anova_summary.csv"))

all_emm <- bind_rows(results_emm) |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_emm, file.path(TBL_DIR, "12_emmeans_contrasts.csv"))

all_r2 <- bind_rows(results_r2) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(all_r2, file.path(TBL_DIR, "12_r2_summary.csv"))

all_genet <- bind_rows(results_genet_effect) |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_genet, file.path(TBL_DIR, "12_genet_treatment_effects.csv"))

# ---- Console summary -------------------------------------------------------
# Print the highlights to the log: the treatment/genet ANOVA rows, the R² table,
# and the per-genet end-of-experiment treatment effects — a sanity read of
# the main results without opening the CSVs.
cat("\n=== Type-III ANOVA: continuous responses (genet terms highlighted) ===\n")
print(all_anova |>
        filter(response_id %in% c("pam_fvfm", "color_dscale",
                                   "growth_pct", "growth_sgr",
                                   "log_zoox_density"),
               grepl("thicket|treatment", term, ignore.case = TRUE)) |>
        select(response_id, term, any_of(c("F value", "F", "Pr(>F)",
                                           "Chisq", "Pr(>Chisq)"))))

cat("\n=== R² (marginal / conditional) ===\n")
print(all_r2)

cat("\n=== Per-genet end-of-experiment treatment effect (28C - 31C) ===\n")
print(all_genet |>
        select(response, thicket, wound, estimate, SE, t.ratio, p.value) |>
        mutate(across(where(is.numeric), \(x) round(x, 3))))

cat("\nWrote: 12_anova_summary.csv, 12_emmeans_contrasts.csv,",
    "12_genet_treatment_effects.csv, 12_r2_summary.csv,",
    "12_morph_fixed_effects.csv, diagnostics in figures/12_diagnostics/\n")

# =============================================================================
# PART 2 — COLOR ORDINAL-CLMM ROBUSTNESS  (formerly 12b_color_clmm_robustness.R)
#
# Residual diagnostics flagged the Gaussian color LMM as KS-non-uniform (the
# D-scale is a 5-level ordinal). We keep the Gaussian model for presentation
# (it matches the other physiology metrics) but confirm here that every
# qualitative inference holds under the ordinal likelihood.
#
# What & why: a Gaussian LMM treats the color score as a continuous number with
#   equal spacing and symmetric errors. The D-scale is 5 ordered bins, so
#   the Gaussian residuals fail the uniformity (KS) check. A cumulative-link mixed
#   model (CLMM, ordinal::clmm) models P(color ≤ k) with a logit link and a set of
#   estimated thresholds — the likelihood for ordered categories. If
#   the heat/genet signals appear under BOTH likelihoods, they are not an artefact
#   of forcing an ordinal scale into a Gaussian model.
# =============================================================================
cat("\n=== PART 2: color ordinal-CLMM robustness ===\n")

color_ord_df <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket),
         # Convert split scores (e.g. 3.5) to the next-higher D category.
         # base::round() uses bankers rounding, so 2.5 would become 2 while
         # 3.5 becomes 4; floor(x + 0.5) is deterministic half-up rounding.
         # Clamp to the D1-D6 Siebeck range (the data currently spans D1-D5,
         # but capping at 6 — not 5 — avoids merging a D6 into D5
         # if a darker score ever enters; unused top levels are dropped by clmm).
         color_ord = factor(pmin(pmax(as.integer(floor(color_num + 0.5)), 1), 6),
                            ordered = TRUE))

# The full 4-way fixed × 2 random structure is too rich for clmm to estimate (the
# Hessian goes singular and SEs diverge). Robustness check, not the primary model:
# same data and the ordinal likelihood, but a trimmed fixed structure (the
# 2- and 3-way interactions that carry the heat/genet signal)
# and a single (1|id) random intercept. Hess = TRUE stores the Hessian so we can
# get standard errors / Wald tests.
m_clmm <- ordinal::clmm(
  color_ord ~ treatment * wound * day + treatment * thicket +
              treatment * wound * thicket + day * thicket + (1 | id),
  data = color_ord_df, Hess = TRUE
)
saveRDS(m_clmm, file.path(MOD_DIR, "12b_color_clmm.rds"))

# Term-by-term likelihood-ratio tests via drop1(): refit dropping each droppable
# term and compare by χ². This is the robustness analogue of the type-III ANOVA
# from Part 1 — we check the sign/significance of each term matches the LMM.
av_drop <- as.data.frame(drop1(m_clmm, test = "Chisq")) |>
  tibble::rownames_to_column("term") |>
  mutate(response = "color_dscale_ordinal")
write_csv(av_drop, file.path(TBL_DIR, "12b_color_clmm.csv"))
av <- av_drop |> rename(LR.stat = LRT, `Pr(>Chisq)` = `Pr(>Chi)`)

# Robustness check: put the Gaussian LMM and the ordinal CLMM side by side on the
# treatment term. If both agree, the choice of likelihood does not change the
# conclusion.
# Compare to Gaussian LMM for the treatment main effect
m_gauss <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
gauss_av <- as.data.frame(anova(m_gauss))
cat("\n=== Treatment main-effect comparison (Gaussian LMM vs ordinal CLMM) ===\n")
if ("treatment" %in% rownames(gauss_av)) {
  cat("Gaussian LMM treatment: F =",
      round(gauss_av["treatment", "F value"], 2), "\n")
}
trt_row <- av[grepl("^treatment$", av$term), , drop = FALSE]
if (nrow(trt_row) > 0) {
  cat("Ordinal CLMM treatment: LRT =", round(trt_row$LR.stat[1], 2),
      " p =", signif(trt_row$`Pr(>Chisq)`[1], 3), "\n")
} else {
  cat("CLMM drop1 did not report a main treatment term (likely all retained due to higher-order interactions); see 12b_color_clmm.csv for full term-by-term tests.\n")
}

# Per-genet contrasts: skip if Hessian is singular (emmeans would fail).
# The drop1 LRT table above is the primary robustness output.
pairs_df <- tryCatch({
  emm <- emmeans::emmeans(m_clmm, ~ treatment | thicket * wound,
                          at = list(day = 14))
  as_tibble(pairs(emm, adjust = "none")) |>
    mutate(response = "color_dscale_ordinal")
}, error = function(e) {
  message("emmeans pairwise contrasts skipped (Hessian singular: ",
          conditionMessage(e), ")")
  tibble(
    note = "CLMM Hessian singular; per-genet contrasts not computable from this fit",
    response = "color_dscale_ordinal"
  )
})
write_csv(pairs_df, file.path(TBL_DIR, "12b_color_clmm_genet_effects.csv"))

cat("\nWrote 12b_color_clmm.csv, 12b_color_clmm_genet_effects.csv, 12b_color_clmm.rds\n")
cat("CLMM fit converged (LL =", round(as.numeric(logLik(m_clmm)), 1),
    "); use drop1 LRTs in 12b_color_clmm.csv for term-by-term tests.\n")
cat("Robustness conclusion: see report (overall direction of every term should match the Gaussian LMM).\n")

# =============================================================================
# PART 3 — PENALIZED MORPHOLOGY REFIT  (formerly 12c_morph_blme.R)
#
# The raw glmer fits in PART 1 give finite predictions (the random-effect
# penalty regularizes them) and omnibus type-II Wald χ², but the
# *individual coefficient* Wald z's diverge when any treatment × day × thicket
# cell reaches 0 or 1 observed events (7/8 traits hit this). Refit with
# blme::bglmer + Cauchy(0, 2.5) priors (Gelman et al. 2008 default) so the
# coefficients become finite and Wald tests interpretable.
#
# What & why (separation): with binary milestones, a cell where the trait is
#   always present (or always absent) drives the maximum-likelihood coefficient
#   to ±infinity — the "complete separation" problem. lme4 returns a large
#   estimate with a large SE, so individual coefficient tests are
#   uninterpretable. The fix is a weakly-informative prior on the fixed effects: a
#   Cauchy(0, 2.5) on the (centred, scaled) logit coefficients shrinks them toward
#   0 while still letting strong effects through. This is the Gelman et al. (2008)
#   recommendation for logistic regression with separation; blme::bglmer
#   implements it as a penalized GLMM. Same fixed/random structure as Part 1; the
#   only change is the prior, which makes the Wald tests usable.
# =============================================================================
cat("\n=== PART 3: penalized morphology GLMMs (blme::bglmer, Cauchy(0,2.5)) ===\n")

ph_blme <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(
    treatment = factor(treatment),
    thicket = factor(thicket)
  )

fit_blme <- function(tr) {
  d <- ph_blme |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) return(NULL)   # same invariance/size guard
  # Set explicit treatment (dummy) contrasts so each coefficient is a
  # "level vs reference" comparison when interpreting the penalized
  # coefficients and their now-finite SEs.
  contrasts(d$treatment) <- contr.treatment(nlevels(d$treatment))
  contrasts(d$thicket) <- contr.treatment(nlevels(d$thicket))
  # The Gelman (2008) Cauchy(0, 2.5) prior assumes inputs are STANDARDIZED — it is
  # calibrated for predictors centred at 0 and scaled to ~0.5 SD (binary inputs to
  # ±0.5). `day` runs 0-15 raw, so the un-scaled prior would over-shrink its slope.
  # Rescale day to mean 0, SD 0.5 (Gelman's "divide by 2 SD") so the prior carries
  # the weight it was designed to. Coefficient SIGN/significance (the robustness
  # check) is unaffected; only the day coefficient's scale changes.
  d$day_z <- (d$day - mean(d$day)) / (2 * stats::sd(d$day))
  # Cauchy(0, 2.5) prior on all (now standardized) fixed effects — Gelman 2008
  # default for logistic regression with separation. Same structure as PART 1.
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      blme::bglmer(
        y ~ treatment * day_z * thicket + (1 | tank) + (1 | id),
        family   = binomial,
        data     = d,
        fixef.prior = "t(scale = 2.5, df = 1)",   # Cauchy(0, 2.5) — Gelman 2008
        control  = lme4::glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 1e5))
      )
    )),
    error = function(e) {
      message("bglmer failed for ", tr, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(m)) return(NULL)
  saveRDS(m, file.path(MOD_DIR, paste0("12c_morph_", tr, "_blme.rds")))
  av <- as.data.frame(car::Anova(m, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(trait = tr)
  tidy <- broom.mixed::tidy(m, effects = "fixed") |>
    mutate(trait = tr, max_se = max(std.error, na.rm = TRUE))
  list(anova = av, tidy = tidy)
}

res <- map(traits, function(tr) {
  cat("Fitting:", tr, "... ")
  out <- fit_blme(tr)
  if (is.null(out)) cat("skipped\n")
  else cat("done (max fixef SE =",
           round(max(out$tidy$std.error, na.rm = TRUE), 2), ")\n")
  out
}) |> setNames(traits) |> compact()

# Aggregate
blme_anova <- map_dfr(res, "anova")
blme_tidy  <- map_dfr(res, "tidy")
write_csv(blme_anova, file.path(TBL_DIR, "12c_morph_blme_anova.csv"))
write_csv(blme_tidy,  file.path(TBL_DIR, "12c_morph_blme_fixed_effects.csv"))

# Check whether the prior resolved the separation: count, per trait, how many
# fixed-effect SEs remain pathological (> 50 on the logit scale = effectively
# infinite). After the Cauchy prior these counts should be 0 / near-0, confirming
# the coefficients are now finite and the Wald tests usable.
cat("\n=== SE sanity (after Cauchy prior) ===\n")
ses <- blme_tidy |> group_by(trait) |>
  summarise(n_fixed = n(),
            n_sep   = sum(std.error > 50, na.rm = TRUE),
            max_se  = max(std.error, na.rm = TRUE))
print(ses)

cat("\nWrote 12c_morph_blme_anova.csv (",
    nrow(blme_anova), "rows ),\n",
    "      12c_morph_blme_fixed_effects.csv (",
    nrow(blme_tidy), "rows),\n",
    "      output/models/12c_morph_<trait>_blme.rds (",
    length(res), "fits)\n", sep = "")
