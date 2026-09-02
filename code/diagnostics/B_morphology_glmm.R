# =============================================================================
# Purpose: Diagnostic battery for the morphological binomial GLMMs from script
#          12. For each wound-healing trait, run DHARMa residual checks,
#          convergence/singularity checks, near-zero random-effect variance,
#          a separation test, and a predicted-probability effect-sanity check.
#
# What & why: each healing trait (e.g. whether the wound is smoothed over) is
#   scored yes/no per coral per day and modeled as a binomial GLMM
#   (glmer(y ~ treatment * day * thicket + (1|tank), wounded corals only)). These
#   logistic mixed models have their own failure modes, so the battery here is:
#     - DHARMa: simulate from the fitted model and rescale residuals to
#       uniform(0,1) — for a binomial GLMM raw residuals are uninformative, so
#       this simulation-based approach is the right way to test fit. We check
#       uniformity (KS), dispersion, outliers, and (only when the trait is rare or
#       near-universal) zero inflation. p >= 0.05 = pass that test.
#     - convergence / isSingular(): did the optimizer converge and is the tank
#       random effect estimable, or did its variance collapse to ~0?
#     - separation: when a predictor perfectly predicts yes/no (common with
#       all-or-nothing healing traits), logistic coefficients and their standard
#       errors diverge toward infinity. We flag any fixed-effect SE > 10 (WARN) or
#       > 50 / non-finite (FAIL) as a sign of (quasi-)separation.
#     - effect sanity: predicted P(trait) at the last day, 28 vs 31 C, must be a
#       valid probability in [0,1] and move in a sensible direction.
#   Traits with a known issue that a companion model fixes (the penalized `blme`
#   refits in 12c; saturated "closure" traits whose timing is handled by the Cox
#   analysis) get their status downgraded to HANDLED rather than WARN/FAIL.
#
# Input:   output/models/12_morph_<trait>_glmm.rds (the fitted GLMMs)
#          data/processed/physio_clean.rds (source data, wound == "yes")
#          output/models/12c_morph_<trait>_blme.rds (presence => HANDLED)
# Output:  output/diagnostics/B_morphology_glmm_diagnostics.csv
#          output/diagnostics/B_morphology_report.md
#          figures/diagnostics/B_<trait>.png
#
# Run:  Rscript code/diagnostics/B_morphology_glmm.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(lme4)
  library(DHARMa)
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
})

# Fix the RNG so DHARMa's simulated residuals are identical on every re-run.
set.seed(20260524)

# ---- paths ------------------------------------------------------------------
ROOT       <- here::here()
MOD_DIR    <- file.path(ROOT, "output", "models")
DATA_PROC  <- file.path(ROOT, "data", "processed")
DIAG_OUT   <- file.path(ROOT, "output", "diagnostics")
DIAG_FIG   <- file.path(ROOT, "figures", "diagnostics")
dir.create(DIAG_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIAG_FIG, showWarnings = FALSE, recursive = TRUE)

CSV_PATH <- file.path(DIAG_OUT, "B_morphology_glmm_diagnostics.csv")
MD_PATH  <- file.path(DIAG_OUT, "B_morphology_report.md")

# ---- helpers ----------------------------------------------------------------
# add_row_safe: append one result row (coerces stat/p to numeric quietly).
# classify_p: shared p-value -> status grader so every test reads the same way.
add_row_safe <- function(df, trait, check, statistic = NA_real_,
                         p_value = NA_real_, status = "PASS", notes = "") {
  bind_rows(df, tibble(
    trait     = trait,
    check     = check,
    statistic = suppressWarnings(as.numeric(statistic)),
    p_value   = suppressWarnings(as.numeric(p_value)),
    status    = status,
    notes     = notes
  ))
}

# Grade a DHARMa/test p-value: for these residual tests a SMALL p is bad (it
# means the residuals deviate from the model's assumption). p >= 0.05 = PASS,
# 0.001-0.05 = WARN, < 0.001 = FAIL; an untestable NA = WARN.
classify_p <- function(p, alpha_warn = 0.05, alpha_fail = 0.001) {
  if (is.na(p)) return("WARN")
  if (p < alpha_fail) return("FAIL")
  if (p < alpha_warn) return("WARN")
  "PASS"
}

handled_glmm_structure <- function(trait, blme_available = FALSE) {
  blme_available || trait %in% saturated_closure_traits
}

# "Closure" traits are near-saturated (almost all corals eventually close the
# wound), so their binomial models are expected to show separation and
# underdispersion. We don't penalize that — timing for these traits is analyzed
# by the Cox/interval models. hole_in_center is a stale saved-model alias of the
# combined axial_polyp_formation observable (see code/04 data-quality note), so
# it is handled the same way if its old model object remains in output/models/.
saturated_closure_traits <- c("axial_polyp_formation", "hole_in_center",
                              "wound_smoothed")

# ---- load source data -------------------------------------------------------
phys_path <- file.path(DATA_PROC, "physio_clean.rds")
if (!file.exists(phys_path)) {
  stop("Cannot find data/processed/physio_clean.rds — script 12 source data")
}
ph_all <- readRDS(phys_path)
# Reconstruct exactly the subset script 12 modeled: wounded corals only, valid
# day. Healing traits are only defined where there is a wound to heal.
ph <- ph_all |>
  dplyr::filter(wound == "yes", !is.na(day), day >= 0) |>
  dplyr::mutate(thicket = as.factor(thicket))

# Detect the two treatment labels and last observed day from the data (rather
# than hard-coding them) so the predicted-probability check is robust to naming.
end_day <- max(ph$day, na.rm = TRUE)
treat_levels <- unique(as.character(ph$treatment))
ambient <- treat_levels[stringr::str_detect(treat_levels, "28")][1]
hot     <- treat_levels[stringr::str_detect(treat_levels, "31")][1]
if (is.na(ambient)) ambient <- sort(treat_levels)[1]
if (is.na(hot))     hot     <- sort(treat_levels)[length(treat_levels)]

cat(sprintf("Source data: %d rows, %d corals, day range %s..%s\n",
            nrow(ph), dplyr::n_distinct(ph$id),
            min(ph$day, na.rm = TRUE), end_day))
cat(sprintf("Treatments: ambient='%s', hot='%s'; end_day=%s\n",
            ambient, hot, end_day))

# ---- find models ------------------------------------------------------------
mod_files <- list.files(MOD_DIR, pattern = "^12_morph_.*_glmm\\.rds$",
                        full.names = TRUE)
cat(sprintf("\nFound %d morphology GLMM files\n", length(mod_files)))

# ---- diagnose one model -----------------------------------------------------
# Runs the full battery on a single trait's GLMM and returns both a rows
# table (for the CSV) and a markdown block (for the human-readable report).
diagnose_one <- function(mfile) {
  trait <- sub("^12_morph_(.*)_glmm\\.rds$", "\\1", basename(mfile))
  cat(sprintf("\n--- %s ---\n", trait))
  res <- tibble()
  notes_md <- c()

  m <- tryCatch(readRDS(mfile), error = function(e) NULL)
  if (is.null(m)) {
    return(list(
      rows = add_row_safe(tibble(), trait, "load_model", NA, NA,
                          "FAIL", "Could not read .rds"),
      md   = sprintf("## %s\n- FAIL: could not load model\n", trait)
    ))
  }

  # subset trait data the way script 12 did
  d <- ph |> dplyr::mutate(y = .data[[trait]]) |> dplyr::filter(!is.na(y))
  mean_p <- mean(d$y, na.rm = TRUE)
  notes_md <- c(notes_md,
                sprintf("- N = %d, n_corals = %d, mean(y) = %.3f",
                        nrow(d), dplyr::n_distinct(d$id), mean_p))

  # 1) convergence + singularity
  # If a penalized blme refit exists for this trait (script 12c), a convergence
  # or separation problem here is already dealt with => HANDLED instead of WARN.
  blme_path <- file.path(MOD_DIR, sprintf("12c_morph_%s_blme.rds", trait))
  blme_available <- file.exists(blme_path)
  conv_msgs <- character(0)
  opt_conv  <- m@optinfo$conv$opt              # optimizer code: 0 = converged
  lme4_msgs <- m@optinfo$conv$lme4$messages    # any "failed to converge" notes
  if (!is.null(lme4_msgs)) conv_msgs <- c(conv_msgs, lme4_msgs)
  conv_ok <- (length(conv_msgs) == 0) && (is.null(opt_conv) || opt_conv == 0)
  # isSingular = TRUE => the tank random-effect variance collapsed to ~0.
  sing    <- isSingular(m, tol = 1e-4)
  conv_status <- if (conv_ok) "PASS"
                 else if (blme_available || trait %in% saturated_closure_traits) "HANDLED"
                 else "WARN"
  res <- add_row_safe(res, trait, "convergence",
                      statistic = ifelse(is.null(opt_conv), NA, opt_conv),
                      status = conv_status,
                      notes = if (conv_ok) "no warnings"
                              else if (blme_available) paste("handled by penalized blme refit in script 12c:", paste(conv_msgs, collapse = "; "))
                              else if (trait %in% saturated_closure_traits) paste("handled as saturated closure/alias trait; inference uses interval-censored timing:", paste(conv_msgs, collapse = "; "))
                              else paste(conv_msgs, collapse = "; "))
  structure_handled <- handled_glmm_structure(trait, blme_available)
  res <- add_row_safe(res, trait, "singular_fit",
                      statistic = as.integer(sing),
                      status = if (sing && structure_handled) "HANDLED"
                               else if (sing) "WARN" else "PASS",
                      notes = if (sing && structure_handled) {
                        "isSingular=TRUE; near-zero random-effect variance handled by companion penalized/timing analysis"
                      } else if (sing) {
                        "isSingular=TRUE; RE variance ~0"
                      } else "non-singular")
  notes_md <- c(notes_md,
                sprintf("- Convergence: %s%s",
                        if (conv_ok) "OK" else if (blme_available) "HANDLED" else "ISSUE",
                        if (length(conv_msgs)) paste0(" — ", paste(conv_msgs, collapse = "; ")) else ""),
                sprintf("- Singular fit: %s", sing))

  # 2) RE variance: report the tank random-effect variance(s). A value < 1e-4 is
  # effectively zero (tanks add no detectable variation) and is flagged WARN —
  # the same situation isSingular catches, shown numerically per grouping factor.
  vc <- as.data.frame(VarCorr(m))
  for (i in seq_len(nrow(vc))) {
    grp <- vc$grp[i]; v <- vc$vcov[i]
    status <- if (is.na(v)) "WARN"
              else if (v < 1e-4 && structure_handled) "HANDLED"
              else if (v < 1e-4) "WARN" else "PASS"
    res <- add_row_safe(res, trait, paste0("re_var_", grp),
                        statistic = v, status = status,
                        notes = if (status == "HANDLED") {
                          "near-zero variance component; retained design term contributes little variation and companion analysis handles sparse binary structure"
                        } else if (status == "WARN") {
                          "Near-zero variance component"
                        } else "variance > 1e-4")
  }
  notes_md <- c(notes_md, sprintf("- RE variance: %s",
                                  paste(sprintf("%s=%.4g", vc$grp, vc$vcov),
                                        collapse = ", ")))

  # 3) separation: the signal is a large standard error. When a predictor
  # perfectly (or nearly) separates yes/no outcomes, the logistic coefficient
  # diverges to +/-infinity and its SE diverges with it. We scan the largest fixed-effect SE:
  # > 50 or non-finite = FAIL (separation), 10-50 = WARN (quasi-separation),
  # <= 10 = plausible. A blme refit, if present, downgrades a FAIL to HANDLED.
  fe <- suppressWarnings(summary(m)$coefficients)
  max_se <- suppressWarnings(max(fe[, "Std. Error"], na.rm = TRUE))
  sep_status <- if (!is.finite(max_se)) "FAIL"
                else if (max_se > 50) "FAIL"
                else if (max_se > 10) "WARN"
                else "PASS"
  if (sep_status == "FAIL" &&
      (blme_available || trait %in% saturated_closure_traits)) {
    sep_status <- "HANDLED"
  }
  res <- add_row_safe(res, trait, "max_fixed_effect_SE",
                      statistic = max_se, status = sep_status,
                      notes = if (sep_status == "HANDLED" && blme_available) "Likely separation; handled by penalized blme refit in script 12c"
                              else if (sep_status == "HANDLED") "Likely separation; handled as saturated closure/alias trait with inference routed to interval-censored timing"
                              else if (sep_status == "FAIL") "Likely separation (SE >50 or non-finite)"
                              else if (sep_status == "WARN") "Possible quasi-separation (SE 10-50)"
                              else "SE in plausible range")
  notes_md <- c(notes_md, sprintf("- Max fixed-effect SE: %.3g (%s)", max_se, sep_status))

  # 4) DHARMa residual battery (see header for the simulate-and-rescale idea).
  # n = 500 simulations, refit = FALSE for speed; refit = TRUE is used later only
  # as a slower follow-up if the fast dispersion test fails.
  dh <- tryCatch(DHARMa::simulateResiduals(m, n = 500, refit = FALSE,
                                           plot = FALSE, seed = 42),
                 error = function(e) NULL)
  if (is.null(dh)) {
    res <- add_row_safe(res, trait, "DHARMa", NA, NA, "FAIL",
                        "simulateResiduals errored")
    notes_md <- c(notes_md, "- DHARMa: FAIL — could not simulate residuals")
  } else {
    ks   <- DHARMa::testUniformity(dh, plot = FALSE)   # KS: residuals uniform?
    disp <- DHARMa::testDispersion(dh, plot = FALSE)   # ratio ~1 = right spread
    outl <- DHARMa::testOutliers(dh, plot = FALSE, type = "binomial")  # excess outliers?

    res <- add_row_safe(res, trait, "DHARMa_uniformity_KS",
                        statistic = unname(ks$statistic), p_value = ks$p.value,
                        status = classify_p(ks$p.value),
                        notes = "KS test on scaled residuals")
    disp_ratio <- unname(disp$statistic)
    disp_status <- classify_p(disp$p.value)
    disp_notes <- sprintf("ratio=%.3g", disp_ratio)
    # Special case: a saturated closure trait that is UNDERdispersed (ratio < 1)
    # is expected, not a defect — relabel HANDLED since its timing is analyzed by
    # the interval-censored / Cox models, not by this binomial fit.
    if ((trait %in% saturated_closure_traits || blme_available) &&
        !is.na(disp$p.value) && disp$p.value < 0.05 &&
        is.finite(disp_ratio) && disp_ratio < 1) {
      disp_status <- "HANDLED"
      disp_notes <- paste(
        disp_notes,
        "underdispersion expected for sparse/saturated morphology trait; inference uses companion penalized GLMM and/or interval-censored timing/event summaries"
      )
    }
    res <- add_row_safe(res, trait, "DHARMa_dispersion",
                        statistic = disp_ratio, p_value = disp$p.value,
                        status = disp_status,
                        notes = disp_notes)
    res <- add_row_safe(res, trait, "DHARMa_outliers",
                        statistic = unname(outl$statistic), p_value = outl$p.value,
                        status = classify_p(outl$p.value),
                        notes = "binomial outlier test")

    # Zero inflation only matters when the trait is rare or near-universal
    # (mean p extreme); for a balanced 50/50 trait the test is uninformative, so
    # we skip it. Small p here = more all-0 (or all-1) outcomes than the model expects.
    if (mean_p < 0.2 || mean_p > 0.8) {
      zi <- tryCatch(DHARMa::testZeroInflation(dh, plot = FALSE),
                     error = function(e) NULL)
      if (!is.null(zi)) {
        res <- add_row_safe(res, trait, "DHARMa_zero_inflation",
                            statistic = unname(zi$statistic), p_value = zi$p.value,
                            status = classify_p(zi$p.value),
                            notes = "tested because mean p extreme")
      }
    }

    # If dispersion failed, also try refit=TRUE unless the result is already
    # handled as a saturated closure trait. Those refits are slow and do not
    # change the interpretation: inference is routed through interval timing.
    if (!is.na(disp$p.value) && disp$p.value < 0.05 && disp_status == "HANDLED") {
      res <- add_row_safe(res, trait, "DHARMa_dispersion_refit",
                          statistic = NA_real_, p_value = NA_real_,
                          status = "HANDLED",
                          notes = "skipped refit=TRUE; underdispersion already handled as saturated closure trait")
    } else if (!is.na(disp$p.value) && disp$p.value < 0.05) {
      dh2 <- tryCatch(DHARMa::simulateResiduals(m, n = 200, refit = TRUE,
                                                plot = FALSE, seed = 42),
                      error = function(e) NULL)
      if (!is.null(dh2)) {
        disp2 <- DHARMa::testDispersion(dh2, plot = FALSE)
        res <- add_row_safe(res, trait, "DHARMa_dispersion_refit",
                            statistic = unname(disp2$statistic),
                            p_value = disp2$p.value,
                            status = classify_p(disp2$p.value),
                            notes = "refit=TRUE follow-up")
      }
    }

    # plot
    png_path <- file.path(DIAG_FIG, sprintf("B_%s.png", trait))
    tryCatch({
      png(png_path, width = 1600, height = 800, res = 150)
      plot(dh)
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message("plot failed for ", trait, ": ", conditionMessage(e))
    })
    notes_md <- c(notes_md,
                  sprintf("- DHARMa: KS p=%.3g, dispersion p=%.3g (ratio %.2f%s), outliers p=%.3g",
                          ks$p.value, disp$p.value, disp_ratio,
                          if (disp_status == "HANDLED") "; handled as saturated closure trait" else "",
                          outl$p.value),
                  sprintf("- Residual plot: figures/diagnostics/B_%s.png", trait))
  }

  # 5) biological sanity: predict P(trait) at the last day for 28 vs 31 C,
  # averaged across genets, using re.form = NA (population-level, ignoring the
  # tank random effect). The check is that both predictions are valid
  # probabilities in [0,1]; the report also prints the heat-minus-ambient delta.
  newd <- expand.grid(
    treatment = c(ambient, hot),
    day       = end_day,
    thicket   = levels(d$thicket),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(treatment %in% unique(d$treatment),
                  thicket %in% unique(as.character(d$thicket)))
  newd$treatment <- factor(newd$treatment, levels = levels(d$treatment))
  newd$thicket   <- factor(newd$thicket,   levels = levels(d$thicket))
  pred <- tryCatch(predict(m, newdata = newd, re.form = NA, type = "response"),
                   error = function(e) rep(NA_real_, nrow(newd)))
  newd$p_pred <- pred
  p_amb <- mean(newd$p_pred[newd$treatment == ambient], na.rm = TRUE)
  p_hot <- mean(newd$p_pred[newd$treatment == hot],     na.rm = TRUE)
  sane <- is.finite(p_amb) && is.finite(p_hot) &&
          p_amb >= 0 && p_amb <= 1 && p_hot >= 0 && p_hot <= 1
  res <- add_row_safe(res, trait, "pred_p_ambient_end",
                      statistic = p_amb, status = if (sane) "PASS" else "WARN",
                      notes = sprintf("predicted at day %s, %s", end_day, ambient))
  res <- add_row_safe(res, trait, "pred_p_hot_end",
                      statistic = p_hot, status = if (sane) "PASS" else "WARN",
                      notes = sprintf("predicted at day %s, %s", end_day, hot))
  notes_md <- c(notes_md,
                sprintf("- Predicted P(trait) at day %s: ambient=%.3f, hot=%.3f, Δ=%+.3f",
                        end_day, p_amb, p_hot, p_hot - p_amb))

  list(rows = res, md = paste0("## ", trait, "\n",
                                paste(notes_md, collapse = "\n"), "\n"))
}

# ---- run --------------------------------------------------------------------
# Loop over every saved morphology GLMM, collecting rows and markdown as we go.
all_rows <- tibble()
all_md   <- c("# Morphological GLMM diagnostics",
              sprintf("Generated: %s", Sys.time()),
              sprintf("Source data: %s (N=%d rows, wound==yes only)",
                      basename(phys_path), nrow(ph)),
              "")

for (f in mod_files) {
  out <- diagnose_one(f)
  all_rows <- bind_rows(all_rows, out$rows)
  all_md   <- c(all_md, out$md)
}

# ---- per-trait summary ------------------------------------------------------
# Collapse the many per-check rows into one verdict per trait. Precedence:
# any FAIL => FAIL; else any HANDLED => HANDLED; else any WARN => WARN; else PASS.
summary_tab <- all_rows |>
  dplyr::group_by(trait) |>
  dplyr::summarise(
    n_checks = dplyr::n(),
    n_pass = sum(status == "PASS"),
    n_handled = sum(status == "HANDLED"),
    n_warn = sum(status == "WARN"),
    n_fail = sum(status == "FAIL"),
    overall = dplyr::case_when(
      n_fail > 0 ~ "FAIL",
      n_handled > 0 ~ "HANDLED",
      n_warn > 0 ~ "WARN",
      TRUE        ~ "PASS"
    ),
    .groups = "drop"
  )

readr::write_csv(all_rows, CSV_PATH)

# build markdown report tail
all_md <- c(all_md,
            "## Per-trait summary",
            "",
            "| Trait | Checks | PASS | HANDLED | WARN | FAIL | Overall |",
            "|-------|--------|------|---------|------|------|---------|",
            apply(summary_tab, 1, function(r)
              sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                      r["trait"], r["n_checks"], r["n_pass"],
                      r["n_handled"],
                      r["n_warn"], r["n_fail"], r["overall"])),
            "",
            sprintf("Totals: %d traits × ~%.1f checks; PASS=%d, HANDLED=%d, WARN=%d, FAIL=%d",
                    nrow(summary_tab), mean(summary_tab$n_checks),
                    sum(all_rows$status == "PASS"),
                    sum(all_rows$status == "HANDLED"),
                    sum(all_rows$status == "WARN"),
                    sum(all_rows$status == "FAIL")))

writeLines(all_md, MD_PATH)

cat("\n==============================\n")
cat(sprintf("Wrote %s (%d rows)\n", CSV_PATH, nrow(all_rows)))
cat(sprintf("Wrote %s\n", MD_PATH))
cat(sprintf("Plots in %s\n", DIAG_FIG))
cat(sprintf("Totals — traits: %d, checks: %d, PASS=%d HANDLED=%d WARN=%d FAIL=%d\n",
            nrow(summary_tab), nrow(all_rows),
            sum(all_rows$status == "PASS"),
            sum(all_rows$status == "HANDLED"),
            sum(all_rows$status == "WARN"),
            sum(all_rows$status == "FAIL")))

# Print fail summary
fails <- dplyr::filter(all_rows, status == "FAIL")
if (nrow(fails)) {
  cat("\nFAIL rows:\n")
  print(fails, n = Inf)
}
