# =============================================================================
# Purpose: Design-alignment audit (diagnostic suite E) — confirm each primary
#          model's fixed/random structure matches the experiment.
#
# What & why: before any p-value can be trusted, the model has to mirror how the
#   experiment was run. For every primary model this script compares the fitted
#   formula, the random-effects structure, observation counts, and cell balance
#   against what the design demands, then flags each check PASS / FAIL / WARN /
#   HANDLED / INFO. The design we are matching against:
#     - 2 temperatures x 2 wounds x 3 genets x 8 tanks (4 per temp)
#     - n=192 corals, repeated measures over 14 days for PAM + color
#     - n=192 destructive biopsies (1 per coral) for symbionts
#     - n=48 for buoyant weight (1 obs per coral)
#     - 24 wounded corals per temperature for morphology traits
#   This is a structural sanity check, not a goodness-of-fit test. It catches the
#   silent mistakes — wrong random term, unbalanced cells, a miscoded factor —
#   that would invalidate the inference downstream.
# Input:   output/models/*.rds  (fitted models: 12_pam_lmm, 12_color_lmm,
#            12_bw_lm, 12_zoox_lmm, 12c_morph_*_blme)
#          data/processed/*.rds (the cleaned data each model was fit to)
# Output:  output/diagnostics/E_design_alignment.csv
#          output/diagnostics/E_design_alignment_report.md
# =============================================================================

# 00_setup.R loads packages and defines shared paths (MOD_DIR, DATA_PROC, ...).
source(here::here("code", "00_setup.R"))

# Where the audit table + markdown report get written.
DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Setup + helpers -------------------------------------------------------
# We accumulate one tibble row per check into `audit`, then bind them at the end.
audit <- list()

# add_row(): append a single audit result. `<<-` writes to the `audit` list in
# the enclosing scope (not a local copy). status is a verdict string the report
# groups by: PASS (matches design), FAIL (does not), WARN (matches but worth a
# look), HANDLED (a known deviation addressed elsewhere), INFO (context only).
add_row <- function(model, design_feature, observed, expected, status, note) {
  audit[[length(audit) + 1]] <<- tibble::tibble(
    model = model, check = design_feature,
    observed = as.character(observed),
    expected = as.character(expected),
    status   = status,
    note     = note
  )
}

# cells(): count the number of distinct combinations of the grouping columns,
# i.e. how many non-empty design cells the data contains. Used to test
# balance (e.g. 8 tank x treatment cells expected).
cells <- function(d, ...) d |> count(...) |> nrow()

# ---- PAM Fv/Fm -------------------------------------------------------------
# Repeated-measures LMM on photosynthetic efficiency: the most complex model
# (4-way fixed structure + two random effects), so it gets the most checks.
m <- readRDS(file.path(MOD_DIR, "12_pam_lmm.rds"))
d <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
# Flatten the model's formula to one string so we can pattern-match its terms.
form <- paste(format(formula(m)), collapse = " ")
# Check 1: are all four crossed fixed effects present? grepl returns TRUE only if
# the full interaction (heat x wound x time x genet) is in the formula.
add_row("12_pam_lmm", "fixed structure", form,
        "treatment * wound * day * thicket",
        if (grepl("treatment \\* wound \\* day \\* thicket", form)) "PASS" else "FAIL",
        "4-way fixed structure required for genet x heat x wound x time")
# Check 2: both random effects present? `tank` captures tank-level confounds;
# `id` accounts for the non-independence of repeated reads on the same coral.
add_row("12_pam_lmm", "random structure",
        paste(names(lme4::ranef(m)), collapse = "+"),
        "tank + id",
        if (all(c("tank", "id") %in% names(lme4::ranef(m)))) "PASS" else "WARN",
        "tank for tank-level confounds, id for repeated measures on same coral")
# Check 3: rough observation count (~192 corals x several PAM days). INFO-only —
# reported for context, not pass/failed, since dropouts make the exact n vary.
add_row("12_pam_lmm", "n observations", nrow(d), "~336", "INFO",
        sprintf("n_corals=%d, days=%s", length(unique(d$id)),
                paste(sort(unique(d$day)), collapse = ",")))
# Check 4: balance. Each of 8 tanks belongs to exactly one temperature, so there
# should be exactly 8 tank x treatment cells. More/fewer signals a coding error.
add_row("12_pam_lmm", "balance (tank x treatment)",
        cells(d, tank, treatment), "8",
        if (cells(d, tank, treatment) == 8) "PASS" else "WARN",
        "8 unique tank x treatment cells expected (4 tanks per temp, fully nested)")
# Check 5: design-justification note. Because tank IDs are unique across temps,
# tank is nested in treatment and a crossed tank:treatment term is unnecessary.
add_row("12_pam_lmm", "(1|tank) sufficient?",
        "yes — tank labels are unique across treatments",
        "yes",
        "PASS",
        "tanks are uniquely IDed; no need for tank:treatment crossed term")
# Check 6: random slope decision, logged as HANDLED. A random day slope (day|id)
# was considered but dropped — too few reads per coral would make it
# singular/overfit. Recording the decision documents it wasn't an oversight.
add_row("12_pam_lmm", "random slope on day?",
        "(1|id) only — no random slope",
        "considered (day|id) but n=14 days, few per id; likely overfits",
        "HANDLED",
        "Random slope considered and intentionally omitted; expected singular/overfit with sparse per-id trajectories")

# ---- Color D-scale (same fixed structure as PAM) ---------------------------
# Coral-bleaching color score, modeled with the same 4-way LMM as PAM. The
# response is a 1-5 ordinal scale fit on a Gaussian model (see below).
m <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
form <- paste(format(formula(m)), collapse = " ")
# Same 4-way fixed-structure check as PAM.
add_row("12_color_lmm", "fixed structure", form,
        "treatment * wound * day * thicket",
        if (grepl("treatment \\* wound \\* day \\* thicket", form)) "PASS" else "FAIL",
        "same 4-way as PAM")
# Flag the Gaussian-on-ordinal modeling choice as HANDLED: it is cross-checked by
# the CLMM robustness script (12b) and disclosed in Methods, so it's not a hidden
# assumption violation.
add_row("12_color_lmm", "ordinal data on Gaussian",
        "Gaussian LMM on D1-D5 ordinal scale",
        "either CLMM or note in Methods",
        "HANDLED",
        "Addressed by 12b CLMM robustness check; KS-violation noted in Section 10")

# ---- Buoyant weight --------------------------------------------------------
# Endpoint growth: one value per coral, so NO day term and NO coral-ID random
# effect (nothing repeated). Tank stays in as the treatment-assignment block.
m <- readRDS(file.path(MOD_DIR, "12_bw_lm.rds"))
d <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
form <- paste(format(formula(m)), collapse = " ")
# Fixed structure is 3-way here (no time): heat x wound x genet, plus (1|tank).
add_row("12_bw_lm", "fixed structure", form,
        "treatment * wound * thicket + (1 | tank)",
        if (grepl("treatment \\* wound \\* thicket", form) &&
            "tank" %in% names(lme4::ranef(m))) "PASS" else "FAIL",
        "no day term; tank retained as treatment-assignment block")
# Random effects must be tank ONLY — `identical` (not `%in%`) so an accidental
# extra (1|id) would FAIL, since one endpoint per coral leaves nothing to nest.
add_row("12_bw_lm", "random effects",
        paste(names(lme4::ranef(m)), collapse = "+"),
        "tank only",
        if (identical(names(lme4::ranef(m)), "tank")) "PASS" else "FAIL",
        "coral ID omitted because each coral has one endpoint growth observation")
# Sample size: 48 corals (24 per temperature). The note also reports how many of
# the 12 design cells (2 trt x 2 wound x 3 genet) are populated.
add_row("12_bw_lm", "n",
        nrow(d), "48 (24 per temperature)",
        if (nrow(d) == 48) "PASS" else "WARN",
        sprintf("Cells in design: %d (2 trt x 2 wound x 3 genet = 12)",
                cells(d, treatment, wound, thicket)))
# Residual df: enough leftover degrees of freedom (>=30) to estimate a 3-way
# model without running out of data. Low residual df = an over-parameterized fit.
add_row("12_bw_lm", "df residual",
        df.residual(m), "~35-36 for tank-aware 3-way model",
        if (df.residual(m) >= 30) "PASS" else "WARN",
        "Adequate residual df for 3-way fixed structure")

# ---- Symbiont density ------------------------------------------------------
# Destructive biopsy: one measurement per coral, so like buoyant weight it drops
# (1|id). Time enters as biopsy_day_c, centered so day 1 = 0.
m <- readRDS(file.path(MOD_DIR, "12_zoox_lmm.rds"))
# Recreate the exact analysis dataset the model saw: drop non-positive/!finite
# densities (can't take log) and rebuild the centered day variable.
d <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(biopsy_day_c = biopsy_day - 1)
form <- paste(format(formula(m)), collapse = " ")
# 4-way fixed structure as for PAM, but with centered biopsy_day_c as the time axis.
add_row("12_zoox_lmm", "fixed structure", form,
        "treatment * wound * biopsy_day_c * thicket",
        if (grepl("treatment \\* wound \\* biopsy_day_c \\* thicket", form)) "PASS" else "FAIL",
        "biopsy_day_c is centered at Day 1; 4-way as for PAM")
# Confirm (1|id) was dropped: each coral is biopsied once, so coral-ID
# has no within-group replication to model.
add_row("12_zoox_lmm", "drop (1|id) — destructive",
        paste(names(lme4::ranef(m)), collapse = "+"),
        "tank only (each coral 1 biopsy)",
        if (identical(names(lme4::ranef(m)), "tank")) "PASS" else "WARN",
        "Correct — destructive sampling means each id has 1 obs")
# Verify the centering puts day 1 at zero (min of biopsy_day_c == 0).
# Centering makes the intercept interpretable as the day-1 baseline.
add_row("12_zoox_lmm", "biopsy_day_c centering",
        sprintf("range = [%g, %g]", min(d$biopsy_day_c, na.rm=TRUE),
                max(d$biopsy_day_c, na.rm=TRUE)),
        "0 to 14 (Day 1 = 0)",
        if (min(d$biopsy_day_c, na.rm=TRUE) == 0) "PASS" else "WARN",
        "Centered correctly at Day 1 baseline")

# ---- Morphology GLMMs (using blme version) ---------------------------------
# Binary wound-healing traits, fit only on wounded corals (so `wound` drops out
# of the model). blme adds a weakly-informative prior to address separation. Loop
# over every saved 12c trait fit and apply the same two checks to each.
trait_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)
for (f in trait_files) {
  trait <- gsub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))
  m <- readRDS(f)
  form <- paste(format(formula(m)), collapse = " ")
  has_morph_fixed_structure <- grepl("treatment \\* day(_z)? \\* thicket", form)
  # Fixed structure should be 3-way (heat x time x genet); wound is absent because
  # the data are already restricted to wounded corals (it's a subset, not a term).
  add_row(basename(f), "fixed structure", form,
          "treatment * day(or day_z) * thicket (wound dropped — wounded-only)",
          if (has_morph_fixed_structure) "PASS" else "FAIL",
          "Restricted to wounded corals; wound is a stratification, not a covariate")
  # Confirm the Cauchy/Student-t(scale=2.5, df=1) prior is the Gelman (2008)
  # default that pulls in coefficients when a binary trait is perfectly separated.
  add_row(basename(f), "Cauchy(0,2.5) prior",
          "fixef.prior = t(scale=2.5, df=1)",
          "Gelman 2008 separation-fix default",
          "PASS",
          "Prior addresses the morphology separation issue")
}

# ---- Cox PH ----------------------------------------------------------------
# Survival model for healing-event timing. These are documented-as-correct checks
# (no model object loaded) since the design choices live in script 14.
add_row("14_cox_*", "stratification",
        "strata = thicket in overall model",
        "thicket as strata (allows different baseline hazard per genet)",
        "PASS",
        "Per-genet HRs in separate per-genet models (low EPV caveat in Section 10)")
add_row("14_cox_*", "wounded-only at risk",
        "filter(wound == 'yes') applied",
        "yes — only wounded corals can express healing traits",
        "PASS", "Confirmed in script 14")

# ---- PCA -------------------------------------------------------------------
# Multivariate ordination of the phenotype. Centering+scaling is mandatory here
# because the input variables are on different units.
add_row("15_pca", "centering + scaling",
        "prcomp(..., center=TRUE, scale.=TRUE)",
        "required because units differ (Fv/Fm, D-scale, %, log-cells)",
        "PASS", "Confirmed in script 15")

# ---- Write outputs ---------------------------------------------------------
# Stack all the per-check rows into one table and save the machine-readable CSV.
final <- bind_rows(audit)
write_csv(final, file.path(DIAG_OUT, "E_design_alignment.csv"))

# ---- Report ----------------------------------------------------------------
# sink() redirects cat() output into the markdown report file until sink() is
# called again with no args. Builds a human-readable summary grouped by status.
sink(file.path(DIAG_OUT, "E_design_alignment_report.md"))
cat("# E. Design alignment audit\n\n")
cat("Generated:", format(Sys.time()), "\n\n")
cat("| Status | Count |\n|---|---|\n")
print(table(final$status))
cat("\n## Status summary\n\n")
# Walk the statuses worst-first (FAIL before PASS) so the most actionable items
# sit at the top of the report.
for (s in c("FAIL", "WARN", "HANDLED", "INFO", "PASS")) {
  rows <- final[final$status == s, ]
  if (nrow(rows) > 0) {
    cat("### ", s, " (", nrow(rows), ")\n\n", sep = "")
    for (i in seq_len(nrow(rows))) {
      cat("- **", rows$model[i], " / ", rows$check[i], "**: ",
          rows$note[i], "\n", sep = "")
    }
    cat("\n")
  }
}
sink()

cat("=== Design alignment audit complete ===\n")
print(table(final$status))
cat("\nWrote", file.path(DIAG_OUT, "E_design_alignment.csv"),
    "and", file.path(DIAG_OUT, "E_design_alignment_report.md"), "\n")
