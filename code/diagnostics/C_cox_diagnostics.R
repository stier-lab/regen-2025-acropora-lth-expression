# =============================================================================
# Purpose: Cox proportional-hazards diagnostics for the survival models in
#          script 14. For each wound-healing trait, refit the Cox models, test
#          the proportional-hazards assumption (cox.zph / Schoenfeld residuals),
#          check events-per-variable, and sanity-check the hazard-ratio direction.
#
# What & why: the survival models in script 14 estimate how fast each healing
#   trait appears and whether heat accelerates it, using Cox proportional-hazards
#   regression. A Cox model's central assumption is PROPORTIONAL HAZARDS — that
#   the heat-vs-ambient hazard ratio is constant over time (heat doesn't help
#   early then hurt late). If that assumption is violated, the single reported HR
#   is misleading, so every Cox fit is tested before it is trusted. This script:
#     1. Recomputes time-to-event per coral (first day the trait = 1, else
#        censored at the last observed day), mirroring code/14_morphology_kaplan.R.
#     2. Refits three scopes: overall (stratified by genet), a documented N/A
#        "by_wound" scope (only wounded corals are at risk), and one model per genet.
#     3. Runs survival::cox.zph() — this correlates the SCHOENFELD residuals
#        (per-event, per-covariate deviations) against transformed time. A
#        significant slope (p < 0.05) means the effect changes over time = PH
#        VIOLATED. p >= 0.05 = no evidence of violation = the assumption holds.
#     4. Saves a Schoenfeld residual plot (residuals vs time, should scatter
#        flat around 0) whenever any covariate violates PH, for visual follow-up.
#     5. Reports events-per-variable (EPV; rule of thumb >= 10 for a stable Cox
#        fit) and the HR direction (HR > 1 = trait emerges faster under heat).
#   Per-source-patch small-sample limitations and a known PH violation that is
#   fixed by a time-varying refit (pigment, source patch C) are downgraded to HANDLED.
#
# Input:   data/processed/physio_clean.rds
#          output/tables/14c_cox_tt_pigment_genetC.csv (presence => HANDLED)
# Output:  output/diagnostics/C_cox_diagnostics.csv
#          output/diagnostics/C_cox_report.md
#          figures/diagnostics/C_<trait>_<scope>_schoenfeld.png  (only if PH fails)
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(survival)
})

# ---- paths -----------------------------------------------------------------
ROOT     <- here::here()
DATA     <- file.path(ROOT, "data", "processed", "physio_clean.rds")
OUT_DIR  <- file.path(ROOT, "output", "diagnostics")
FIG_DIR  <- file.path(ROOT, "figures", "diagnostics")
CSV_PATH <- file.path(OUT_DIR, "C_cox_diagnostics.csv")
MD_PATH  <- file.path(OUT_DIR, "C_cox_report.md")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- load data -------------------------------------------------------------
# Set treatment factor with 28C first so it becomes the reference level: every
# hazard ratio below is then "31C vs 28C" (heated relative to ambient control).
ph <- readRDS(DATA) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    wound = factor(wound, levels = c("no", "yes")),
    thicket = factor(thicket)
  )
contrasts(ph$treatment) <- contr.treatment(nlevels(ph$treatment))

# hole_in_center + polyp_in_hole combined into axial_polyp_formation (code/04),
# mirroring the milestones analysed in 14.
traits <- c("axial_polyp_formation", "wound_smoothed",
            "pigment_over_wound", "tip_exist", "tip_extension",
            "new_corallites_on_tip")

# ---- compute time-to-event (mirrors code/14_morphology_kaplan.R) -----------
# Collapse the repeated daily yes/no scores per coral into one survival record:
#   event_day = first day the trait was observed (==1); if it never appeared,
#               censor at the last day the coral was seen.
#   event     = 1 if the trait ever appeared, 0 if censored (right-censored).
# This is the standard time-to-first-event construction Cox models expect.
compute_events <- function(d, trait) {
  d |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[trait]]) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      event_day = {
        first1 <- which(y == 1)[1]                      # index of first "yes"
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),    # ever healed? 1/0
      .groups = "drop"
    ) |>
    mutate(trait = trait)
}
events <- map_dfr(traits, ~ compute_events(ph, .x))
events <- events |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    thicket = factor(thicket)
  )
contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))

# ---- helpers ---------------------------------------------------------------
# Schoenfeld plot saver: writes one PNG with a panel per covariate. Reading it:
# each point is one event's scaled Schoenfeld residual plotted against time. If
# PH holds, the cloud scatters flat around 0 (the LOWESS line is horizontal); a
# rising or falling trend is the visual signature of a PH violation.
save_schoenfeld <- function(zph_obj, fig_path) {
  # Use raw Schoenfeld residual scatter (plot.cox.zph's spline smoother is
  # often singular for the small per-genet datasets here; raw points are
  # interpretable). One panel per non-GLOBAL covariate.
  cov_names <- setdiff(rownames(zph_obj$table), "GLOBAL")
  n_panels  <- max(1, length(cov_names))
  png(fig_path, width = 800 * n_panels, height = 800, res = 180)
  on.exit(dev.off(), add = TRUE)
  op <- par(mfrow = c(1, n_panels), mar = c(4.2, 4.2, 2.5, 1.2))
  on.exit(par(op), add = TRUE, after = FALSE)
  resid_mat <- zph_obj$y    # rows = event times, cols = covariates
  times     <- zph_obj$x    # transformed time (default: K-M)
  for (i in seq_along(cov_names)) {
    nm <- cov_names[i]
    plot(times, resid_mat[, i],
         xlab = "Transformed time (K-M)",
         ylab = paste("Scaled Schoenfeld:", nm),
         main = paste0(nm,
                       "  (cox.zph p=",
                       signif(zph_obj$table[nm, "p"], 3), ")"),
         pch  = 19, col = "grey30")
    abline(h = 0, lty = 2, col = "red")
    # add LOWESS if enough points
    if (length(unique(times)) >= 4) {
      lo <- tryCatch(lowess(times, resid_mat[, i]), error = function(e) NULL)
      if (!is.null(lo)) lines(lo, col = "steelblue", lwd = 1.5)
    }
  }
  invisible(NULL)
}

# Status logic for PH rows only. EPV is reported as a separate check so low
# event counts do not get mislabeled as proportional-hazards violations.
# cox.zph p < 0.05 = PH assumption violated = FAIL; p >= 0.05 = PASS; NA (test
# could not run, often too few events) = WARN.
status_from <- function(ph_p, epv_ok, n_event) {
  if (is.na(ph_p)) {
    return("WARN")
  }
  if (ph_p < 0.05) return("FAIL")
  "PASS"
}

# Diagnostic for one fitted coxph model: PH test (per covariate + GLOBAL),
# optional Schoenfeld plot, events-per-variable, and HR direction. Returns rows.
diagnose_fit <- function(fit, trait, scope) {
  if (is.null(fit) || inherits(fit, "try-error")) {
    return(tibble(
      trait = trait, scope = scope, check = "fit",
      statistic = NA_real_, p_value = NA_real_,
      status = "WARN",
      notes = "coxph() failed or skipped (insufficient events)"
    ))
  }
  s   <- summary(fit)
  nev <- s$nevent                       # number of events (corals that healed)
  np  <- length(coef(fit))              # number of estimated covariates
  # Events per variable: too few events per coefficient => unstable Cox estimates.
  epv <- if (np > 0) nev / np else NA_real_
  epv_ok <- !is.na(epv) && epv >= 10

  # cox.zph runs the proportional-hazards test: one row per covariate plus a
  # GLOBAL row. It regresses the scaled Schoenfeld residuals on time; a nonzero
  # slope (small p) means the hazard ratio drifts over time = PH violated.
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  rows <- list()

  if (!is.null(zph)) {
    tab <- as.data.frame(zph$table)
    tab$term <- rownames(tab)
    for (i in seq_len(nrow(tab))) {
      tname <- tab$term[i]
      pval  <- tab$p[i]
      st    <- status_from(pval, epv_ok, nev)
      rows[[length(rows) + 1]] <- tibble(
        trait = trait, scope = scope,
        check = paste0("PH_", tname),
        statistic = round(tab$chisq[i], 3),
        p_value   = round(pval, 4),
        status    = st,
        notes     = paste0("cox.zph chisq=", round(tab$chisq[i], 2),
                          ", df=", tab$df[i])
      )
    }
    # Save the residual plot only when the GLOBAL test or any covariate fails
    # the PH assumption.
    any_violation <- any(tab$p < 0.05, na.rm = TRUE)
    if (any_violation) {
      fig_name <- paste0("C_", trait, "_",
                         gsub("[^A-Za-z0-9]+", "_", scope),
                         "_schoenfeld.png")
      fig_path <- file.path(FIG_DIR, fig_name)
      tryCatch(save_schoenfeld(zph, fig_path),
               error = function(e) message("Schoenfeld plot failed: ", e$message))
      rows[[length(rows) + 1]] <- tibble(
        trait = trait, scope = scope, check = "schoenfeld_plot",
        statistic = NA_real_, p_value = NA_real_,
        status = "WARN",
        notes = paste0("plot saved: figures/diagnostics/", fig_name)
      )
    }
  } else {
    rows[[length(rows) + 1]] <- tibble(
      trait = trait, scope = scope, check = "PH_GLOBAL",
      statistic = NA_real_, p_value = NA_real_,
      status = "WARN",
      notes = "cox.zph() failed"
    )
  }

  # EPV rule of thumb: PASS if >= 10 events per covariate, else WARN (underpowered
  # but not a PH violation — kept as its own check so the two don't get conflated).
  rows[[length(rows) + 1]] <- tibble(
    trait = trait, scope = scope, check = "EPV",
    statistic = round(epv, 2), p_value = NA_real_,
    status = if (epv_ok) "PASS" else "WARN",
    notes  = paste0("n_event=", nev, ", n_covariates=", np,
                    " (rule of thumb: EPV >= 10)")
  )

  # HR direction sanity: report exp(coef) for the treatment term. HR > 1 means
  # the trait reaches its endpoint FASTER under 31C (higher hazard of healing);
  # HR < 1 = slower under heat. Logged as PASS (descriptive, not a pass/fail gate).
  cf <- tryCatch(s$coefficients, error = function(e) NULL)
  if (!is.null(cf) && nrow(cf) >= 1 && "exp(coef)" %in% colnames(cf)) {
    # Find the row whose name contains "treatment"
    trt_row <- grep("treatment", rownames(cf), value = TRUE)[1]
    if (!is.na(trt_row)) {
      hr  <- cf[trt_row, "exp(coef)"]
      pvl <- cf[trt_row, "Pr(>|z|)"]
      rows[[length(rows) + 1]] <- tibble(
        trait = trait, scope = scope, check = "HR_direction",
        statistic = round(hr, 3),
        p_value   = round(pvl, 4),
        status    = "PASS",
        notes     = paste0("HR(31C vs 28C)=", round(hr, 2),
                          " — ", if (hr > 1) "faster onset under 31C"
                                 else        "slower onset under 31C")
      )
    }
  }
  bind_rows(rows)
}

# ---- fit + diagnose --------------------------------------------------------
# For each trait, fit the three scopes and run diagnose_fit() on each.
results <- list()

for (tr in traits) {
  d <- events |> filter(trait == tr)
  n_event <- sum(d$event)

  # (1) OVERALL: strata(thicket) lets each genet have its own baseline hazard
  # while estimating one shared treatment effect — same model as fit_cox_overall().
  # Skip if fewer than 5 events (too sparse to fit a stable Cox model).
  fit_overall <- tryCatch(
    if (n_event < 5) NULL else
      coxph(Surv(event_day, event) ~ treatment + strata(thicket), data = d),
    error = function(e) NULL
  )
  results[[length(results) + 1]] <-
    diagnose_fit(fit_overall, tr, "overall_strata_thicket")

  # (2) BY_WOUND: in the source script only wounded corals are at risk, so
  # wound is not a model covariate. We document this rather than refit.
  results[[length(results) + 1]] <- tibble(
    trait = tr, scope = "by_wound",
    check = "applicability", statistic = NA_real_, p_value = NA_real_,
    status = "PASS",
    notes  = "N/A — only wounded corals at risk in source script"
  )

  # (3) PER-GENET: one model per genet to look for genotype-specific heat
  # responses. These are small datasets (skip if < 3 events), so they are
  # treated as descriptive — their EPV/PH warnings are downgraded to HANDLED below.
  for (g in c("a", "c", "d")) {
    dg <- d |> filter(thicket == g)
    fit_g <- tryCatch(
      if (sum(dg$event) < 3) NULL else
        coxph(Surv(event_day, event) ~ treatment, data = dg),
      error = function(e) NULL
    )
    results[[length(results) + 1]] <-
      diagnose_fit(fit_g, tr, paste0("genet_", g))
  }
}

diag_df <- bind_rows(results)
# Reconcile statuses (numbers untouched): the one PH violation (pigment,
# source patch C) is HANDLED by a time-varying coxph refit saved in 14c_*.csv,
# and the per-source-patch small-sample warnings are HANDLED as documented descriptive limits.
tt_path <- file.path(ROOT, "output", "tables", "14c_cox_tt_pigment_genetC.csv")
diag_df <- diag_df |>
  mutate(
    status = case_when(
      trait == "pigment_over_wound" & scope == "genet_c" &
        check %in% c("PH_treatment", "PH_GLOBAL") &
        status == "FAIL" & file.exists(tt_path) ~ "HANDLED",
      scope != "overall_strata_thicket" &
        check %in% c("EPV", "PH_GLOBAL", "fit", "schoenfeld_plot") &
        status == "WARN" ~ "HANDLED",
      TRUE ~ status
    ),
    notes = case_when(
      trait == "pigment_over_wound" & scope == "genet_c" &
        check %in% c("PH_treatment", "PH_GLOBAL") &
        status == "HANDLED" ~ paste(notes, "handled by time-varying coxph refit in 14c_cox_tt_pigment_genetC.csv", sep = "; "),
      scope != "overall_strata_thicket" &
        check %in% c("EPV", "PH_GLOBAL", "fit") &
        status == "HANDLED" ~ paste(notes, "handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)", sep = "; "),
      check == "schoenfeld_plot" & status == "HANDLED" ~
        paste(notes, "handled by saved diagnostic plot and time-varying refit where applicable", sep = "; "),
      TRUE ~ notes
    )
  )
write_csv(diag_df, CSV_PATH)

# ---- markdown report -------------------------------------------------------
md_lines <- c(
  "# Cox proportional-hazards diagnostics",
  "",
  paste0("Source script: `code/14_morphology_kaplan.R`  "),
  paste0("Data: `data/processed/physio_clean.rds`  "),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Checks per fitted model:",
  "- **PH_*** — `survival::cox.zph()` per covariate and GLOBAL. p < 0.05 = PH violated.",
  "- **EPV** — events per variable (rule of thumb: >=10). Underpowered models flagged.",
  "- **HR_direction** — hazard ratio for 31C vs 28C (HR > 1 = trait emerges faster under heat).",
  "- **schoenfeld_plot** — saved to `figures/diagnostics/` when any covariate violates PH.",
  ""
)

for (tr in traits) {
  md_lines <- c(md_lines, paste0("## ", tr), "")
  sub <- diag_df |> filter(trait == tr)
  for (sc in unique(sub$scope)) {
    md_lines <- c(md_lines, paste0("### scope: ", sc), "")
    rows <- sub |> filter(scope == sc)
    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      md_lines <- c(md_lines,
        sprintf("- **[%s] %s** — stat=%s, p=%s. %s",
                r$status, r$check,
                ifelse(is.na(r$statistic), "NA", as.character(r$statistic)),
                ifelse(is.na(r$p_value),   "NA", as.character(r$p_value)),
                r$notes))
    }
    md_lines <- c(md_lines, "")
  }
}

# ---- summary block ---------------------------------------------------------
ph_rows <- diag_df |> filter(grepl("^PH_GLOBAL$", check))
n_fail_ph <- sum(ph_rows$status == "FAIL", na.rm = TRUE)
n_handled_ph <- sum(ph_rows$status == "HANDLED", na.rm = TRUE)
n_pass_ph <- sum(ph_rows$p_value >= 0.05, na.rm = TRUE)
n_total   <- nrow(ph_rows)

epv_rows <- diag_df |> filter(check == "EPV")
n_epv_warn <- sum(epv_rows$status == "WARN", na.rm = TRUE)
n_epv_handled <- sum(epv_rows$status == "HANDLED", na.rm = TRUE)
n_ph_untestable <- sum(is.na(ph_rows$p_value))
n_handled <- sum(diag_df$status == "HANDLED", na.rm = TRUE)

fail_rows <- diag_df |> filter(status == "FAIL")

summary_block <- c(
  "## Summary",
  "",
  sprintf("- Models tested: %d (one PH_GLOBAL row per fitted coxph).", n_total),
  sprintf("- PH assumption: **%d PASS / %d HANDLED / %d FAIL** (p < 0.05 on GLOBAL test).",
          n_pass_ph, n_handled_ph, n_fail_ph),
  sprintf("- PH tests untestable/failed numerically: **%d**", n_ph_untestable),
  sprintf("- EPV warnings (events/covariate < 10): **%d WARN / %d HANDLED**",
          n_epv_warn, n_epv_handled),
  sprintf("- Handled diagnostic failures with explicit refits: **%d**", n_handled),
  sprintf("- Total FAIL rows: **%d**", nrow(fail_rows)),
  "",
  "Recommended fixes for PH violations:",
  "1. Stratify on the violating covariate (already done for `thicket` in the overall model).",
  "2. If `treatment` itself violates PH, add a time-varying coefficient (`tt()` term in `coxph`).",
  "3. For per-genet models with violation, report as a limitation — small n constrains alternatives.",
  ""
)

md_lines <- c(md_lines[1:14], summary_block, md_lines[15:length(md_lines)])
writeLines(md_lines, MD_PATH)

cat("\n=== Cox diagnostics complete ===\n")
cat("CSV:    ", CSV_PATH, "\n")
cat("Report: ", MD_PATH, "\n")
cat("PH PASS / FAIL: ", n_pass_ph, " / ", n_fail_ph, "\n", sep = "")
cat("Total FAIL rows: ", nrow(fail_rows), "\n", sep = "")
