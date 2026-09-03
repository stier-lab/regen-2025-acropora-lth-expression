# =============================================================================
# Purpose: Aggregate every statistical result in the repo into one tidy
#          spreadsheet that the manuscript can pull from directly. One row
#          per (response, term, contrast, model) with a unified schema so
#          when an analysis updates, the paper text can be updated by
#          re-reading this single CSV.
#
#          Schema:
#            domain            (Physiology | Morphology | Survival | Multivariate)
#            response          PAM Fv/Fm, Color D-scale, Growth %, etc.
#            model_type        LMM | GLMM (binomial) | LM | Cox PH | PCA | LRT
#            term              Treatment, Wound, Day, Treatment:Wound, ...
#                                or contrast label (e.g. "28C - 31C | a | wounded")
#            test              ANOVA F | LRT chi-sq | Wald z | t | HR | LRT
#            statistic         numeric test stat
#            df1, df2          degrees of freedom (df1 = numerator/test df,
#                                df2 = denominator if Satterthwaite/F)
#            n                 sample size used for the test
#            estimate          model-scale estimate (effect, mean diff, log HR)
#            units             units / scale ("Fv/Fm units", "D-scale", "%",
#                                "ln cells cm-2", "log HR", etc.)
#            pct_change        signed % change of 31C relative to 28C
#                                (where computable on natural scale)
#            ci_low, ci_high   95% CI for the estimate
#            p_value           p-value
#            qualitative       short biological-language summary
#            source_script     code file that produced the row
#            source_artifact   CSV / model RDS it came from
#
# What & why: this is the single source of truth for the manuscript's numbers.
#   Every other script writes its own results CSV in its own shape; here we read
#   them all back in and recast each into the one shared row schema described
#   above, so that every effect size, test statistic, df, p-value and CI cited
#   in the paper can be traced to exactly one row — and to the script that made
#   it (the source_script / source_artifact columns). Each "Block" below is an
#   adapter with the same job: read one upstream table, rename/compute its
#   columns into the common schema, and tag its domain. The blocks are then
#   row-bound into `master` and written out twice — a tidy machine-readable CSV
#   and a paper-ready formatted CSV. Nothing here fits a model; it only
#   transcribes and harmonizes results that already exist. Blocks wrapped in
#   `if (file.exists(...))` are optional sensitivity analyses that may not have
#   been run; they contribute zero rows (an empty tibble) when their
#   source CSV is absent, so the table degrades gracefully.
# Input:   output/tables/*.csv, output/models/*.rds
# Output:  output/tables/20_master_results.csv             — tidy spreadsheet
#          output/tables/20_master_results_paper_ready.csv — manuscript-formatted
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Helpers --------------------------------------------------------------
# Formatters and lookups shared by the blocks below.

# Render a p-value in journal format: "<0.001" below that threshold, else 3 dp.
fmt_p   <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
# Fixed-decimal number formatter (NA-safe) for the paper-ready table.
fmt_num <- function(x, d = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
}

# Turn an estimate + p-value into a plain-language direction sentence. Defined
# here as a reusable qualitative-summary helper; the blocks below currently
# build their own bespoke sentences inline rather than calling this.
qual_dir <- function(estimate, p, response_label,
                     positive = "higher", negative = "lower",
                     ref = "31 °C vs 28 °C") {
  if (is.na(estimate) || is.na(p)) return(NA_character_)
  if (p >= 0.05) return(paste0("no detectable difference (", ref, ")"))
  direction <- if (estimate > 0) positive else negative
  paste0(response_label, " ", direction, " under ", ref)
}

# Lookup: internal response_id -> manuscript-ready response label.
response_label_map <- c(
  pam_fvfm         = "PAM Fv/Fm",
  color_dscale     = "Color (Siebeck D)",
  growth_pct       = "Growth (% mass change)",   # from 12 (LMM ANOVA, emmeans)
  pct_growth       = "Growth (% mass change)",   # from 05 (tank-permutation / descriptive lm)
  log_zoox_density = "ln symbionts cm-2"
)

# Lookup: internal response_id -> natural (back-transformed) units for display.
natural_units <- c(
  pam_fvfm         = "Fv/Fm",
  color_dscale     = "D-scale units",
  growth_pct       = "% over 15 d",
  pct_growth       = "% over 15 d",
  log_zoox_density = "ln cells cm-2"
)

baseline_means <- function() {
  # Compute mean of the response at 28C (end of experiment) so the per-genet
  # treatment effects can be expressed as % change. Used for pct_change column.
  # Returns a named list (one tibble of thicket x wound baselines per response).
  # Note: log_zoox_density gets an EMPTY baseline tibble on purpose — its %
  # change is derived analytically from the log estimate (see attach_pct), not
  # from a baseline mean, so no 28C mean is needed.
  pam  <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
  col  <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
  bw   <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
  zx   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))

  list(
    pam_fvfm     = pam |>
      filter(day == max(day, na.rm = TRUE), treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(fv_fm, na.rm = TRUE), .groups = "drop"),
    color_dscale = col |>
      filter(day == max(day, na.rm = TRUE), treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(color_num, na.rm = TRUE), .groups = "drop"),
    growth_pct = bw |>
      filter(treatment == "28C", is.finite(pct_growth)) |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(pct_growth, na.rm = TRUE), .groups = "drop"),
    log_zoox_density = tibble(
      thicket = character(),
      wound = factor(levels = c("no", "yes")),
      baseline = numeric()
    )
  )
}

# Compute the baselines once, up front, for reuse by Block 2's % change.
bases <- baseline_means()

# ===========================================================================
# Block 1 — Continuous LMM ANOVA (script 12)
# ===========================================================================
# Block 1 covers continuous responses only. Morph rows in 12_anova are
# duplicates of those in 04_morphology_trait_anova_genet.csv (same model,
# different script that writes the table) — drop them here and pick up
# the morphology ANOVAs in Block 4 below.
# Emits one row per omnibus fixed-effect test (Treatment, Wound, Day, and their
# interactions); test column is "ANOVA F" when an F value exists (lmerTest /
# car F-tests) and "Wald chi-sq" otherwise.
anova12 <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"),
                    show_col_types = FALSE) |>
  filter(!grepl("^morph_", response_id),
         !grepl("^morph_", response)) |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = if_else(response_id %in% c("growth_pct"),
                              "LM", "LMM"),
    term            = term,
    test            = if_else(is.na(`F value`), "Wald chi-sq", "ANOVA F"),
    statistic       = coalesce(`F value`, Chisq),
    df1             = coalesce(NumDF, Df),   # numerator/test df (lmerTest NumDF; car Df)
    df2             = DenDF,                 # denominator df (Satterthwaite) for F-tests
    n               = n_obs,
    estimate        = NA_real_,
    units           = NA_character_,
    pct_change      = NA_real_,
    ci_low          = NA_real_,
    ci_high         = NA_real_,
    p_value         = coalesce(`Pr(>F)`, `Pr(>Chisq)`),
    qualitative     = paste0(term, " effect on ", response,
                             if_else(coalesce(`Pr(>F)`, `Pr(>Chisq)`) < 0.05,
                                     " (significant)", " (n.s.)")),
    source_script   = "code/12_models.R",
    source_artifact = "output/tables/12_anova_summary.csv"
  )

# ===========================================================================
# Block 2 — Per-genet treatment contrasts at end of experiment (script 12)
# ===========================================================================
# One row per genet x wound contrast (28C - 31C) at Day 14, with estimate, SE,
# and the % change of 31C relative to 28C attached for natural-scale responses.
genet_eff <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                      show_col_types = FALSE)

# Attach baseline and compute pct change. The estimate is 28C - 31C, so the
# heat-driven change of 31C relative to 28C is -estimate (hence the minus
# signs). Three back-transform rules by response type:
#   growth_pct       -> estimate is already a % difference, so % change = -estimate
#   log_zoox_density -> estimate is on the log scale, so % change = exp(-estimate)-1
#   everything else  -> divide -estimate by the 28C baseline mean
attach_pct <- function(df, response_key) {
  base_df <- bases[[response_key]]
  if (is.null(base_df)) {
    df$pct_change <- NA_real_
    return(df)
  }
  df |>
    left_join(base_df, by = c("thicket", "wound")) |>
    mutate(pct_change = if (response_key == "growth_pct") {
      -estimate
    } else if (response_key == "log_zoox_density") {
      (exp(-estimate) - 1) * 100
    } else {
      -estimate / baseline * 100
    })
}

# Apply attach_pct to each response separately (each needs its own baseline),
# then stack the results back together.
genet_eff_pct <- bind_rows(
  attach_pct(filter(genet_eff, response == "pam_fvfm"),       "pam_fvfm"),
  attach_pct(filter(genet_eff, response == "color_dscale"),   "color_dscale"),
  attach_pct(filter(genet_eff, response == "growth_pct"),   "growth_pct"),
  attach_pct(filter(genet_eff, response == "log_zoox_density"),"log_zoox_density")
)

# Recast the per-genet contrasts into the common schema. test = Wald z when a
# z-ratio is present (e.g. GLMM/emmeans on z scale), else Satterthwaite t; the
# 95% CI is built from the estimate +/- 1.96*SE.
genet_rows <- genet_eff_pct |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = if_else(response %in% c("growth_pct"),
                              "LM", "LMM"),
    term            = sprintf("contrast: 28C - 31C | genet=%s | wound=%s",
                              thicket, wound),
    test            = if_else(!is.na(z.ratio), "Wald z", "Satterthwaite t"),
    statistic       = coalesce(t.ratio, z.ratio),
    df1             = df,
    df2             = NA_real_,
    n               = NA_real_,
    estimate        = estimate,
    units           = natural_units[response] |> coalesce(""),
    pct_change      = pct_change,
    ci_low          = estimate - 1.96 * SE,
    ci_high         = estimate + 1.96 * SE,
    p_value         = p.value,
    qualitative     = sprintf(
      "Genet %s, %s: %s under 31C %s",
      thicket, if_else(wound == "yes", "wounded", "unwounded"),
      if_else(p.value < 0.05,
              if_else(estimate > 0, "lower", "higher"),
              "no change"),
      if_else(!is.na(pct_change), sprintf("(%+.1f%%)", pct_change), "")
    ),
    source_script   = "code/12_models.R",
    source_artifact = "output/tables/12_genet_treatment_effects.csv"
  )

# ===========================================================================
# Block 3 — R² per response (script 12)
# ===========================================================================
# Marginal (fixed effects only) and conditional (fixed + random) R^2 from each
# LMM. pivot_longer makes one row per R^2 metric per response.
r2_rows <- read_csv(file.path(TBL_DIR, "12_r2_summary.csv"),
                    show_col_types = FALSE) |>
  pivot_longer(c(R2_marginal, R2_conditional),
               names_to = "metric", values_to = "value") |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = "LMM",
    term            = metric,
    test            = "R-squared",
    statistic       = value,
    df1             = NA_real_, df2 = NA_real_, n = NA_real_,
    estimate        = value, units = "proportion",
    pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value         = NA_real_,
    qualitative     = sprintf("%s = %.3f", metric, value),
    source_script   = "code/12_models.R",
    source_artifact = "output/tables/12_r2_summary.csv"
  )

# ===========================================================================
# Block 4 — Morphology GLMM ANOVA (penalized fits)
# ===========================================================================
# Morphology traits are binary (trait expressed yes/no) and prone to perfect
# separation, so they are fit with blme (a weakly-informative Cauchy(0,2.5)
# prior on the fixed effects) instead of plain glmer. This first table is the
# omnibus Wald chi-sq per fixed effect.
morph_anova <- read_csv(file.path(TBL_DIR, "12c_morph_blme_anova.csv"),
                        show_col_types = FALSE) |>
  transmute(
    domain          = "Morphology",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "GLMM (binomial, blme Cauchy(0,2.5))",
    term            = term,
    test            = "Wald chi-sq",
    statistic       = Chisq,
    df1             = Df,
    df2             = NA_real_, n = NA_real_,
    estimate        = NA_real_, units = NA_character_,
    pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value         = `Pr(>Chisq)`,
    qualitative     = paste0(term, " effect on ", trait,
                             if_else(`Pr(>Chisq)` < 0.05,
                                     " (significant)", " (n.s.)")),
    source_script   = "code/12_models.R",
    source_artifact = "output/tables/12c_morph_blme_anova.csv"
  )

# Primary morphology fixed-effects: pull from the blme (Cauchy-penalized)
# refit so that Wald tests are interpretable under separation. Each coefficient
# is an odds ratio after exp(); pct_change is (OR - 1) * 100.
# A `separated` flag catches coefficients where even the prior could not tame
# the estimate (giant SE or implausibly large effect). For those rows we blank
# out the unstable estimate/CI/p and keep only the omnibus chi-sq from above,
# so the table never reports an artefactual "significant" odds ratio.
morph_fixed <- read_csv(
  file.path(TBL_DIR, "12c_morph_blme_fixed_effects.csv"),
  show_col_types = FALSE
) |>
  filter(term != "(Intercept)") |>
  mutate(separated = std.error > 50 | abs(estimate) > 10) |>
  transmute(
    domain          = "Morphology",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "GLMM (binomial, blme Cauchy(0,2.5))",
    term            = paste0("fixed: ", term),
    test            = "Wald z",
    statistic       = statistic,
    df1             = NA_real_, df2 = NA_real_, n = NA_real_,
    estimate        = if_else(separated, NA_real_, estimate),
    units           = if_else(separated, "log-odds [residual separation]",
                              "log-odds"),
    pct_change      = if_else(separated, NA_real_, (exp(estimate) - 1) * 100),
    ci_low          = if_else(separated, NA_real_, estimate - 1.96 * std.error),
    ci_high         = if_else(separated, NA_real_, estimate + 1.96 * std.error),
    p_value         = if_else(separated, NA_real_,
                              2 * pnorm(-abs(statistic))),
    qualitative     = if_else(
      separated,
      sprintf("%s: residual separation after Cauchy(0,2.5) prior — report omnibus chi-sq only", term),
      sprintf("%s: OR=%.2f (%s, blme penalized)", term, exp(estimate),
              if_else(2 * pnorm(-abs(statistic)) < 0.05, "sig", "n.s."))
    ),
    source_script   = "code/12_models.R",
    source_artifact = "output/tables/12c_morph_blme_fixed_effects.csv"
  )

# ===========================================================================
# Block 5 — Cox PH hazard ratios (script 14)
# ===========================================================================
# Wound-healing milestones (e.g. hole closure, polyp/tip regeneration) are
# time-to-event outcomes, modelled with Cox proportional hazards. HR_31_vs28 is
# the hazard ratio for reaching the milestone under 31C vs 28C: HR > 1 = faster
# onset under heat, HR < 1 = delayed. Filter drops non-finite/degenerate HRs.
cox_rows <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                     show_col_types = FALSE) |>
  filter(is.finite(HR_31_vs28), HR_31_vs28 > 0, HR_31_vs28 < Inf) |>
  transmute(
    domain          = "Survival",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "Cox PH",
    term            = sprintf("HR 31C/28C [%s]", scope),
    test            = "Wald z",
    statistic       = z,
    df1             = 1, df2 = NA_real_, n = n,
    estimate        = HR_31_vs28,
    units           = "hazard ratio",
    pct_change      = (HR_31_vs28 - 1) * 100,
    ci_low          = HR_lo,
    ci_high         = HR_hi,
    p_value         = p,
    qualitative     = sprintf(
      "%s: %s under 31C (HR=%.2f, %s)",
      trait,
      case_when(HR_31_vs28 > 1.1 & p < 0.05 ~ "faster onset",
                HR_31_vs28 < 0.9 & p < 0.05 ~ "delayed onset",
                TRUE                         ~ "no detectable change"),
      HR_31_vs28,
      if_else(p < 0.05, "sig", "n.s.")
    ),
    source_script   = "code/14_morphology_kaplan.R",
    source_artifact = "output/tables/14_cox_hazard_ratios.csv"
  )

# Interval-censored Weibull AFT alternative to the Cox model: because traits
# were only observed on discrete survey days, the true event time lies in an
# interval. time_ratio_31_vs28 > 1 means heat delays first expression.
interval_rows <- if (file.exists(file.path(TBL_DIR, "14_interval_survreg.csv"))) {
  read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"),
           show_col_types = FALSE) |>
    filter(is.finite(time_ratio_31_vs28), time_ratio_31_vs28 > 0) |>
    transmute(
      domain          = "Survival",
      response        = str_to_sentence(gsub("_", " ", trait)),
      model_type      = "interval-censored Weibull AFT",
      term            = "time ratio 31C/28C",
      test            = "Wald z",
      statistic       = z,
      df1             = 1, df2 = NA_real_, n = n,
      estimate        = time_ratio_31_vs28,
      units           = "time ratio",
      pct_change      = (time_ratio_31_vs28 - 1) * 100,
      ci_low          = ratio_lo,
      ci_high         = ratio_hi,
      p_value         = p,
      qualitative     = sprintf(
        "%s: %s at 31C (time ratio=%.2f, %s). %s",
        trait,
        case_when(time_ratio_31_vs28 > 1.1 & p < 0.05 ~ "later first expression",
                  time_ratio_31_vs28 < 0.9 & p < 0.05 ~ "earlier first expression",
                  TRUE                                  ~ "no detectable timing change"),
        time_ratio_31_vs28,
        if_else(p < 0.05, "sig", "n.s."),
        event_interpretation
      ),
      source_script   = "code/14_morphology_kaplan.R",
      source_artifact = "output/tables/14_interval_survreg.csv"
    )
} else tibble()

# Time-varying-coefficient Cox refit for the one case (pigmentation, source patch C)
# where the standard Cox model violated the proportional-hazards assumption.
# tt(treatment) lets the treatment effect change with log(t+1); a non-sig
# result here means the earlier per-genet HR was inflated by the PH violation.
cox_tt <- if (file.exists(file.path(TBL_DIR, "14c_cox_tt_pigment_genetC.csv"))) {
  read_csv(file.path(TBL_DIR, "14c_cox_tt_pigment_genetC.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain="Survival",
      response=str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH (time-varying)",
      term=sprintf("tt(treatment) [%s]", scope),
      test="Wald z",
      statistic=z,
      df1=1, df2=NA_real_, n=n,
      estimate=coef,
      units="log-HR slope vs log(t+1)",
      pct_change=NA_real_,
      ci_low=coef - 1.96 * se, ci_high=coef + 1.96 * se,
      p_value=p,
      qualitative=sprintf("PH-corrected per-genet test (%s): %s. %s",
                          scope,
                          if_else(p < 0.05, "sig", "n.s. — original per-genet HR was inflated by PH violation"),
                          note),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14c_cox_tt_pigment_genetC.csv"
    )
} else tibble()

# Likelihood-ratio tests for whether genet matters in the Cox models: one row
# for the genet main effect, one for the genet x treatment interaction. The two
# transmute() calls below pull those two tests out of the same source table.
cox_genet_lrt <- read_csv(file.path(TBL_DIR, "14_cox_genet_LRT.csv"),
                          show_col_types = FALSE)
cox_genet_rows <- bind_rows(
  cox_genet_lrt |>
    transmute(
      domain="Survival",
      response = str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH",
      term="LRT: genet main effect",
      test="LRT chi-sq",
      statistic=chisq_genet_main,
      df1=df_genet_main, df2=NA_real_, n=NA_real_,
      estimate=NA_real_, units=NA_character_,
      pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=p_genet_main,
      qualitative=sprintf("Genet main effect on %s healing: %s", trait,
                          if_else(p_genet_main < 0.05, "sig", "n.s.")),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14_cox_genet_LRT.csv"
    ),
  cox_genet_lrt |>
    transmute(
      domain="Survival",
      response = str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH",
      term="LRT: genet x treatment",
      test="LRT chi-sq",
      statistic=chisq_genet_trt_int,
      df1=df_genet_trt_int, df2=NA_real_, n=NA_real_,
      estimate=NA_real_, units=NA_character_,
      pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=p_genet_trt_int,
      qualitative=sprintf("Genet x treatment on %s healing: %s", trait,
                          if_else(p_genet_trt_int < 0.05, "sig", "n.s.")),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14_cox_genet_LRT.csv"
    )
)

# ===========================================================================
# Block 6 — Genet LRT comparisons (script 13)
# ===========================================================================
# Does adding the genet x treatment (x wound x day) interaction improve fit?
# One LRT per continuous response, comparing the null model to the full model.
# estimate carries the delta-AIC; the term string notes the exact interaction.
lrt13 <- read_csv(file.path(TBL_DIR, "13_genet_anova.csv"),
                  show_col_types = FALSE) |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = "LRT (null vs genet x trt)",
    term            = if_else(response == "growth_pct",
                              "LRT: adding genet x treatment x wound",
                              "LRT: adding genet x treatment x wound x day"),
    test            = "LRT chi-sq",
    statistic       = lrt_chisq,
    df1             = lrt_df, df2 = NA_real_, n = n_obs,
    estimate        = delta_aic,
    units           = "delta AIC",
    pct_change      = NA_real_,
    ci_low          = NA_real_, ci_high = NA_real_,
    p_value         = lrt_p,
    qualitative     = sprintf(
      "Adding genet x trt improves fit by dAIC=%.1f (%s)", delta_aic,
      if_else(lrt_p < 0.05, "sig", "n.s.")
    ),
    source_script   = "code/13_genet_interaction.R",
    source_artifact = "output/tables/13_genet_anova.csv"
  )

# ===========================================================================
# Block 7 — PCA + genet displacement (script 15)
# ===========================================================================
# Two products of the end-of-experiment multivariate PCA. First: the variable
# loadings on each principal component (pivoted to one row per variable x PC).
pca_load <- read_csv(file.path(TBL_DIR, "15_pca_loadings.csv"),
                     show_col_types = FALSE) |>
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "loading") |>
  transmute(
    domain="Multivariate", response="End-of-experiment physiology PCA",
    model_type="PCA",
    term=sprintf("loading[%s] -> %s", variable, PC),
    test="PCA loading",
    statistic=loading,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=loading, units="loading",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("%s loads %.2f on %s", variable, loading, PC),
    source_script="code/15_multivariate.R",
    source_artifact="output/tables/15_pca_loadings.csv"
  )

# Second: how far each genet's centroid moves in PC space from 28C to 31C —
# a single multivariate "how much did heat shift this genet" distance.
pca_disp <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"),
                     show_col_types = FALSE) |>
  transmute(
    domain="Multivariate", response="Per-genet thermal displacement (PCA)",
    model_type="PCA",
    term=sprintf("genet=%s", thicket),
    test="Euclidean displacement",
    statistic=displacement,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=displacement, units="PC units",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("Genet %s shifts %.2f PC-units under heat", thicket,
                        displacement),
    source_script="code/15_multivariate.R",
    source_artifact="output/tables/15_genet_pca_displacement.csv"
  )

# ===========================================================================
# Block 8 — Buoyant weight LM (script 05 — single fixed-effect table)
# ===========================================================================
# Coral-level descriptive linear model of areal calcification; one row per
# fixed-effect coefficient (intercept dropped), with t-test and 95% CI.
bw_lm <- read_csv(file.path(TBL_DIR, "05_buoyant_weight_lm.csv"),
                  show_col_types = FALSE) |>
  filter(term != "(Intercept)") |>
  transmute(
    domain="Physiology", response="Growth (% mass change)",
    model_type="LM (coral-level descriptive)",
    term=paste0("fixed: ", term),
    test="t",
    statistic=statistic,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=estimate, units="% over 15 d",
    pct_change=NA_real_, ci_low=conf.low, ci_high=conf.high,
    p_value=p.value,
    qualitative=sprintf("%s: est=%.2f%% (%s)", term, estimate,
                        if_else(p.value < 0.05, "sig", "n.s.")),
    source_script="code/05_buoyant_weight.R",
    source_artifact="output/tables/05_buoyant_weight_lm.csv"
  )

# Tank-level randomization test: the experimental unit for temperature is the
# tank (n = 8), not the coral, so this is the design-correct treatment test that
# avoids pseudoreplication from many corals within a tank.
bw_tank_test <- if (file.exists(file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"))) {
  read_csv(file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain="Physiology", response="Growth (% mass change)",
      model_type="tank-level randomization",
      term="treatment: 28C - 31C",
      test=test,
      statistic=estimate_28_minus_31,
      df1=NA_real_, df2=NA_real_, n=n_tanks,
      estimate=estimate_28_minus_31, units="% over 15 d",
      pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=p_two_sided,
      qualitative=sprintf("Tank-level heat effect %.2f pct-points over 15 d (p=%.3g)",
                          estimate_28_minus_31, p_two_sided),
      source_script="code/05_buoyant_weight.R",
      source_artifact="output/tables/05_buoyant_weight_tank_test.csv"
    )
} else tibble()

# ===========================================================================
# Block 9 — Color CLMM ordinal robustness (script 12b)
# ===========================================================================
# Robustness check: colour was modelled as continuous in Block 1, but it is an
# ordered Siebeck D-scale. A cumulative-link mixed model (CLMM) re-tests each
# effect treating colour as ordinal; matching significance confirms the main
# analysis is not an artefact of treating an ordinal score as continuous.
clmm_rows <- if (file.exists(file.path(TBL_DIR, "12b_color_clmm.csv"))) {
  read_csv(file.path(TBL_DIR, "12b_color_clmm.csv"),
           show_col_types = FALSE) |>
    filter(!is.na(LRT)) |>
    transmute(
      domain          = "Physiology",
      response        = "Color (Siebeck D) — ordinal robustness",
      model_type      = "CLMM (cumulative-link mixed)",
      term            = term,
      test            = "LRT chi-sq",
      statistic       = LRT,
      df1             = Df, df2 = NA_real_, n = NA_real_,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = `Pr(>Chi)`,
      qualitative     = sprintf("%s (ordinal CLMM): %s", term,
                                if_else(`Pr(>Chi)` < 0.05, "sig", "n.s.")),
      source_script   = "code/12_models.R",
      source_artifact = "output/tables/12b_color_clmm.csv"
    )
} else {
  tibble()
}

# ===========================================================================
# Block 10 — Raw summary stats for the manuscript narrative
#            (cells / growth means, so the prose cites the table, not ad-hoc numbers)
# ===========================================================================
# Raw per-treatment buoyant-weight growth means (mean, SD, n) so the prose can
# cite this table instead of recomputing numbers ad hoc.
bw_raw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  group_by(treatment) |>
  summarise(mean_pct = mean(pct_growth, na.rm = TRUE),
            sd_pct   = sd(pct_growth, na.rm = TRUE),
            n        = sum(!is.na(pct_growth)),
            .groups  = "drop")

bw_means_rows <- bw_raw |>
  transmute(
    domain="Physiology", response="Buoyant weight growth (%)",
    model_type="raw summary",
    term=sprintf("mean(%s)", treatment),
    test="mean +/- SD",
    statistic=mean_pct,
    df1=NA_real_, df2=NA_real_, n=n,
    estimate=mean_pct, units="%",
    pct_change=NA_real_, ci_low=mean_pct - sd_pct, ci_high=mean_pct + sd_pct,
    p_value=NA_real_,
    qualitative=sprintf("%s mean growth = %.2f%% (SD %.2f, n=%d)",
                        treatment, mean_pct, sd_pct, n),
    source_script="code/05_buoyant_weight.R (raw)",
    source_artifact="data/processed/buoyant_weight_clean.rds"
  )

# Derived headline number: the % reduction in mass gain at 31C relative to 28C
# (a reported number), computed straight from the two means.
bw_pct_drop <- bw_raw |>
  summarise(
    pct_drop_31C_vs_28C = (mean_pct[treatment == "31C"] -
                           mean_pct[treatment == "28C"]) /
                           mean_pct[treatment == "28C"] * 100
  ) |>
  transmute(
    domain="Physiology", response="Buoyant weight growth (%)",
    model_type="derived",
    term="percent reduction 31C vs 28C", test="ratio",
    statistic=pct_drop_31C_vs_28C,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=pct_drop_31C_vs_28C, units="%",
    pct_change=pct_drop_31C_vs_28C, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("31C corals grew %.0f%% less mass than 28C corals",
                        abs(pct_drop_31C_vs_28C)),
    source_script="code/05_buoyant_weight.R (raw)",
    source_artifact="data/processed/buoyant_weight_clean.rds"
  )

# Raw symbiont-density means at the first and last biopsy day per treatment, so
# the narrative can cite starting vs ending zooxanthellae densities directly.
zoox_raw <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0)

zoox_means_rows <- zoox_raw |>
  group_by(treatment, biopsy_day) |>
  summarise(mean_cells = mean(cells_per_cm2, na.rm = TRUE),
            n = n(), .groups = "drop") |>
  filter(biopsy_day %in% c(min(biopsy_day), max(biopsy_day))) |>
  transmute(
    domain="Physiology", response="Symbiont density (cells cm-2)",
    model_type="raw summary",
    term=sprintf("mean(%s, day %s)", treatment, biopsy_day),
    test="mean", statistic=mean_cells,
    df1=NA_real_, df2=NA_real_, n=n,
    estimate=mean_cells, units="cells cm-2",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("%s, biopsy day %s: %.2e cells cm-2 (n=%d)",
                        treatment, biopsy_day, mean_cells, n),
    source_script="code/06_symbiont_chl.R (raw)",
    source_artifact="data/processed/symbiont_chl_clean.rds"
  )

# ===========================================================================
# Block 11 — Time-series diagnostics (script 23): autocorrelation,
#            random slope, nonlinearity-of-time for repeated-measures responses
# ===========================================================================
# Checks that the repeated-measures LMMs are well specified (residual
# autocorrelation, need for random slopes, nonlinear time). df is read as text
# upstream, so as.numeric() coerces it (suppressWarnings hides NA coercions).
ts_rows <- if (file.exists(file.path(TBL_DIR, "23_timeseries_diagnostics.csv"))) {
  read_csv(file.path(TBL_DIR, "23_timeseries_diagnostics.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain          = "Time-series diagnostic",
      response        = response,
      model_type      = "diagnostic",
      term            = check,
      test            = "LRT / Wald",
      statistic       = statistic,
      df1             = suppressWarnings(as.numeric(df)), df2 = NA_real_,
      n               = NA_real_,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = p_value,
      qualitative     = conclusion,
      source_script   = "code/sensitivity/23_timeseries_diagnostics.R",
      source_artifact = "output/tables/23_timeseries_diagnostics.csv"
    )
} else tibble()

# ===========================================================================
# Block 12 — Cox proportional-hazards tests (script 14, all overall models)
# ===========================================================================
# Schoenfeld-residual (cox.zph) test of the proportional-hazards assumption for
# every Cox model in Block 5. ph_ok flags whether the assumption held; failures
# point to the time-varying refit captured in cox_tt above.
coxph_rows <- if (file.exists(file.path(TBL_DIR, "14_cox_ph_tests.csv"))) {
  read_csv(file.path(TBL_DIR, "14_cox_ph_tests.csv"), show_col_types = FALSE) |>
    transmute(
      domain          = "Survival diagnostic",
      response        = str_to_sentence(gsub("_", " ", trait)),
      model_type      = "cox.zph (Schoenfeld)",
      term            = "proportional-hazards test (treatment)",
      test            = "cox.zph chi-sq",
      statistic       = zph_chisq,
      df1             = 1, df2 = NA_real_, n = n_event,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = zph_p,
      qualitative     = if_else(ph_ok, "PH assumption met (p>=0.05)",
                                "PH violated — see time-varying refit"),
      source_script   = "code/14_morphology_kaplan.R",
      source_artifact = "output/tables/14_cox_ph_tests.csv"
    )
} else tibble()

# ===========================================================================
# Block 13 — Thermal context vs Cunning et al. 2024 acute ED50 (script 26)
# ===========================================================================
# External benchmark (not a result of this study): published acute-heat ED50
# values for A. pulchra at Mahana, used to place the 31C chronic treatment in
# context (it sits below the acute lethal threshold, i.e. sublethal stress).
thermal_rows <- if (file.exists(file.path(TBL_DIR, "26_thermal_context.csv"))) {
  read_csv(file.path(TBL_DIR, "26_thermal_context.csv"), show_col_types = FALSE) |>
    transmute(
      domain          = "Thermal context (external benchmark)",
      response        = "A. pulchra acute Fv/Fm ED50 (Cunning et al. 2024, Mahana)",
      model_type      = "CBASS reference",
      term            = metric,
      test            = "reference value",
      statistic       = value,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = value, units = "°C or count",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      qualitative = "acute thermal-tolerance benchmark; 31C chronic is sublethal (below acute ED50)",
      source_script = "code/sensitivity/26_thermal_context.R",
      source_artifact = "data/external/cunning2024_apulchra_ed50.csv"
    )
} else tibble()

# ===========================================================================
# Block 14 — Inter-milestone lag (script 14): closure -> regeneration
# ===========================================================================
# Median number of days between a wound closing and the tissue then
# regenerating, per treatment, plus the % of wounds that closed but never
# regenerated (pct_closed_no_regen). CI columns carry the IQR of the lag.
lag_rows <- if (file.exists(file.path(TBL_DIR, "14_milestone_lag_summary.csv"))) {
  read_csv(file.path(TBL_DIR, "14_milestone_lag_summary.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Healing milestone lag",
      response = pair, model_type = "event timing",
      term = paste0(treatment, " (n closed=", n_closed,
                    ", reached both=", n_reached_both, ")"),
      test = "median lag (days)", statistic = median_lag,
      df1 = NA_real_, df2 = NA_real_, n = n_reached_both,
      estimate = median_lag, units = "days",
      pct_change = pct_closed_no_regen, ci_low = iqr_low, ci_high = iqr_high,
      p_value = NA_real_,
      qualitative = sprintf("%s: %.0f%% closed but never regenerated; median lag %s d",
                            treatment, pct_closed_no_regen,
                            ifelse(is.na(median_lag), "NA", as.character(round(median_lag,1)))),
      source_script = "code/14_morphology_kaplan.R",
      source_artifact = "output/tables/14_milestone_lag_summary.csv"
    )
} else tibble()

# ===========================================================================
# Block 15 — Variance partitioning / ICC (script 27)
# ===========================================================================
# Intraclass correlations: the fraction of variance each random effect (tank,
# genet, colony, ...) explains in each model — i.e. where the variability lives.
icc_rows <- if (file.exists(file.path(TBL_DIR, "27_variance_partitioning.csv"))) {
  read_csv(file.path(TBL_DIR, "27_variance_partitioning.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Variance partitioning", response = model,
      model_type = "ICC", term = paste0("ICC[", component, "]"),
      test = "variance fraction", statistic = icc,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = icc, units = "fraction of variance",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      qualitative = sprintf("%s: %.1f%% of variance", component, 100 * icc),
      source_script = "code/sensitivity/27_variance_partitioning.R",
      source_artifact = "output/tables/27_variance_partitioning.csv"
    )
} else tibble()

# ===========================================================================
# Block 16 — Multiple-testing sensitivity (script 28)
# ===========================================================================
# Documents which hypotheses were a priori/confirmatory (reported on raw p) vs
# exploratory (Benjamini-Hochberg FDR-adjusted), with the rationale for each, so
# reviewers can see the multiple-comparisons stance is principled, not post hoc.
mt_rows <- if (file.exists(file.path(TBL_DIR, "28_multiple_testing.csv"))) {
  read_csv(file.path(TBL_DIR, "28_multiple_testing.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Confirmatory / exploratory testing", response = hypothesis,
      model_type = "a priori (raw) / exploratory (BH)", term = test,
      test = "reported p", statistic = p_reported,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = p_reported, units = "p (raw if a priori, BH if exploratory)",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = p_value,
      qualitative = sprintf("%s; %s (rationale: %s)",
                            hypothesis,
                            ifelse(significant, "significant", "n.s."),
                            rationale),
      source_script = "code/sensitivity/28_multiple_testing.R",
      source_artifact = "output/tables/28_multiple_testing.csv"
    )
} else tibble()

# ===========================================================================
# Block 17 — Morphology probability-scale contrasts (script 29)
# ===========================================================================
# Back-transforms the binary morphology GLMMs to the probability scale so the
# heat effect reads as an absolute difference in the probability of trait
# expression (delta_prob = P(28C) - P(31C)) — more interpretable than the
# log-odds in Block 4. Also carries the odds ratio for cross-reference.
probc_rows <- if (file.exists(file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"))) {
  read_csv(file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"),
           show_col_types = FALSE) |>
    filter(is.finite(delta_prob)) |>
    transmute(
      domain = "Morphology", response = str_to_sentence(gsub("_", " ", trait)),
      model_type = "GLMM contrast (prob scale)",
      term = sprintf("28C - 31C at day %d", day),
      test = "delta probability", statistic = delta_prob,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = delta_prob, units = "probability",
      pct_change = 100 * delta_prob, ci_low = NA_real_, ci_high = NA_real_,
      p_value = p_value,
      qualitative = sprintf("P(28C)=%.2f vs P(31C)=%.2f; OR(31vs28)=%.2g",
                            prob_28, prob_31, or_31_vs_28),
      source_script = "code/sensitivity/29_morphology_prob_contrasts.R",
      source_artifact = "output/tables/29_morphology_prob_contrasts.csv"
    )
} else tibble()

# ===========================================================================
# Block 18 — Composite source-patch resilience scores (script 19)
# ===========================================================================
# The cross-response synthesis from script 19: each source patch's mean standardized
# heat sensitivity (lower = more resilient; source patch C wins). NB the ci_low/ci_high
# columns are repurposed here to carry the median sensitivity and the PCA
# displacement (not a confidence interval) so they ride along in the same row.
resilience_rows <- if (file.exists(file.path(TBL_DIR, "19_genet_resilience_summary.csv"))) {
  read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Genet resilience",
      response = "Composite heat sensitivity across responses",
      model_type = "row-max standardized composite",
      term = sprintf("genet=%s", thicket),
      test = "mean standardized sensitivity",
      statistic = mean_sensitivity,
      df1 = NA_real_, df2 = NA_real_, n = n_responses,
      estimate = mean_sensitivity,
      units = "standardized sensitivity",
      pct_change = NA_real_,
      ci_low = median_sensitivity,
      ci_high = pca_displacement,
      p_value = NA_real_,
      qualitative = sprintf(
        "Genet %s: mean sensitivity %.2f, median %.2f, PCA displacement %.2f, rank %.0f",
        thicket, mean_sensitivity, median_sensitivity, pca_displacement,
        rank_overall
      ),
      source_script = "code/19_genet_dashboard.R",
      source_artifact = "output/tables/19_genet_resilience_summary.csv"
    )
} else tibble()

# Same composite, but split by scope (heat-only vs heat-while-wounded) to show
# whether a genet's resilience is general or specific to the wounded condition.
resilience_scope_rows <- if (file.exists(file.path(TBL_DIR, "19c_resilience_decomp_by_scope.csv"))) {
  read_csv(file.path(TBL_DIR, "19c_resilience_decomp_by_scope.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Genet resilience",
      response = "Composite heat sensitivity by scope",
      model_type = "row-max standardized composite",
      term = sprintf("genet=%s; %s", thicket, scope),
      test = "mean standardized sensitivity",
      statistic = mean_sensitivity,
      df1 = NA_real_, df2 = NA_real_, n = n_responses,
      estimate = mean_sensitivity,
      units = "standardized sensitivity",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      qualitative = sprintf("Genet %s, %s: mean sensitivity %.2f across %d responses",
                            thicket, scope, mean_sensitivity, n_responses),
      source_script = "code/19_genet_dashboard.R",
      source_artifact = "output/tables/19c_resilience_decomp_by_scope.csv"
    )
} else tibble()

# ===========================================================================
# Block 19 — Microscope-only morphology validation (script 11)
# ===========================================================================
# The microscope/photo cohort is a separate 16-coral validation dataset, not
# pooled with the main morphology models. Include only the exploratory endpoint
# and first-observed timing checks in the master table so manuscript users can
# cite them as supporting microscope evidence without confusing them for the
# main n=48 morphology inference.
micro_rows <- if (file.exists(file.path(TBL_DIR, "11_microscope_event_tests.csv"))) {
  micro_tests <- read_csv(file.path(TBL_DIR, "11_microscope_event_tests.csv"),
                          show_col_types = FALSE)
  bind_rows(
    micro_tests |>
      transmute(
        domain = "Microscope morphology",
        response = label,
        model_type = "photo-only endpoint Fisher exact",
        term = sprintf("endpoint D%s: 31C - 28C", endpoint_day),
        test = "Fisher exact",
        statistic = NA_real_,
        df1 = 1, df2 = NA_real_, n = n_28 + n_31,
        estimate = risk_diff_31_minus_28,
        units = "proportion difference",
        pct_change = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_,
        p_value = fisher_p,
        qualitative = sprintf(
          "Microscope endpoint: %s at 28C vs %s at 31C (risk diff=%+.2f; exploratory)",
          sprintf("%d/%d", yes_28, n_28),
          sprintf("%d/%d", yes_31, n_31),
          risk_diff_31_minus_28
        ),
        source_script = "code/11_microscope_physio.R",
        source_artifact = "output/tables/11_microscope_event_tests.csv"
      ),
    micro_tests |>
      transmute(
        domain = "Microscope morphology",
        response = label,
        model_type = "photo-only first-observed timing",
        term = "log-rank: 31C vs 28C",
        test = "log-rank chi-sq",
        statistic = logrank_chisq,
        df1 = logrank_df, df2 = NA_real_, n = n_28 + n_31,
        estimate = NA_real_,
        units = NA_character_,
        pct_change = NA_real_,
        ci_low = NA_real_, ci_high = NA_real_,
        p_value = logrank_p,
        qualitative = sprintf(
          "Microscope first-observed timing for %s (exploratory log-rank p=%s)",
          trait, fmt_p(logrank_p)
        ),
        source_script = "code/11_microscope_physio.R",
        source_artifact = "output/tables/11_microscope_event_tests.csv"
      )
  )
} else tibble()

# ===========================================================================
# Combine and write
# ===========================================================================
# Stack every block into one long table (optional blocks contribute nothing
# when empty), round all numeric columns to 4 dp, sort, then add a single
# human-readable `description` string per row and move it to the front.
master <- bind_rows(anova12, genet_rows, r2_rows,
                    morph_anova, morph_fixed,
                    interval_rows, cox_rows, cox_genet_rows, cox_tt,
                    lrt13, pca_load, pca_disp, bw_lm, bw_tank_test,
                    clmm_rows, bw_means_rows, bw_pct_drop, zoox_means_rows,
                    ts_rows, coxph_rows, thermal_rows,
                    lag_rows, icc_rows, mt_rows, probc_rows,
                    resilience_rows, resilience_scope_rows,
                    micro_rows) |>
  mutate(across(c(statistic, estimate, pct_change, ci_low, ci_high, p_value),
                \(x) round(x, 4))) |>
  arrange(domain, response, model_type, term) |>
  # Single human-readable description of WHAT each analysis is — the "what was
  # tested" sentence — so the table reads as: [description] | stat | df | est | p | ...
  mutate(description = sprintf("%s: %s — %s [%s, %s]",
                              domain, response, term, model_type, test)) |>
  relocate(description)

# The tidy, programmatic source of truth — every downstream consumer reads this.
write_csv(master, file.path(TBL_DIR, "20_master_results.csv"))

# ---- Paper-ready formatted version ---------------------------------------
# A second, presentation-formatted copy: pretty column names, p-values via
# fmt_p(), df collapsed to "df1, df2" when a denominator df exists, estimates
# carrying their units, and CIs/%-changes as formatted strings — drop-in for a
# manuscript table.
paper_ready <- master |>
  transmute(
    Description   = description,
    Domain        = domain,
    Response      = response,
    Model         = model_type,
    Term          = term,
    Test          = test,
    Statistic     = fmt_num(statistic, 2),
    df            = if_else(!is.na(df2),
                            sprintf("%s, %s", fmt_num(df1,0), fmt_num(df2,1)),
                            fmt_num(df1, 0)),
    N             = fmt_num(n, 0),
    Estimate      = if_else(!is.na(units),
                            paste(fmt_num(estimate, 3), units),
                            fmt_num(estimate, 3)),
    `% change`    = if_else(!is.na(pct_change),
                            sprintf("%+.1f%%", pct_change),
                            NA_character_),
    `95% CI`      = if_else(!is.na(ci_low),
                            sprintf("[%.2f, %.2f]", ci_low, ci_high),
                            NA_character_),
    p             = fmt_p(p_value),
    Interpretation = qualitative,
    Source         = source_script
  )
write_csv(paper_ready, file.path(TBL_DIR, "20_master_results_paper_ready.csv"))

# ---- Console summary ------------------------------------------------------
# Quick sanity report: total rows and the breakdown by domain and model type,
# so a run confirms every block contributed and nothing silently dropped out.
cat("\n=== Master results table ===\n")
cat("Rows:", nrow(master), "\n")
cat("By domain:\n")
print(table(master$domain))
cat("\nBy model type:\n")
print(table(master$model_type))
cat("\nWrote:\n",
    " - output/tables/20_master_results.csv             (tidy, programmatic)\n",
    " - output/tables/20_master_results_paper_ready.csv (formatted)\n")
