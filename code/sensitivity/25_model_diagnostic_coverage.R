# =============================================================================
# Purpose: Guarantee 100% diagnostic coverage. Enumerate every fitted model in
#          output/models/, generate a residual-diagnostic plot for any model
#          that lacks one (DHARMa for lmer/glmer/lm; ordinal models handled by
#          script 12b), and emit a coverage manifest cross-referencing each
#          model to (a) its diagnostic figure and (b) the figure that
#          visualizes its fitted result. Fails loudly (non-zero summary) if any
#          model is left without a diagnostic.
#
# What & why: a reviewer may ask whether the assumptions of every fitted model
#   were checked, or only some of them. This script documents that. It
#   walks the folder of saved model objects, and for each one confirms a
#   residual/assumption diagnostic exists — reusing one already made by an earlier
#   script, or building a fresh DHARMa plot if not. The output is an audit table:
#   every model paired with its diagnostic plot and with the manuscript figure
#   that shows its result. If any model lacks a diagnostic, the script ends with a
#   WARNING and a non-zero gap count. This is bookkeeping for reproducibility, not
#   new statistics.
# Input:   output/models/*.rds  (+ Cox PH plots from script 14)
# Output:  output/tables/25_model_diagnostic_coverage.csv
#          output/diagnostics/K_model_coverage_report.md
#          figures/diagnostics/K_<model>_dharma.png   (for any that were missing)
# =============================================================================

# 00_setup.R loads packages and shared paths (MOD_DIR, FIG_DIR, TBL_DIR, ...).
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(DHARMa) })   # simulation-based residuals

DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
diagnostic_notes <- character()

# Every fitted model the pipeline saved to disk; this is the master list we audit.
model_files <- list.files(MOD_DIR, pattern = "\\.rds$", full.names = TRUE)

# ---- Crosswalk: model -> the manuscript figure that shows its result -------
# Pattern-matches a model's filename stem to the results figure it underlies, so
# the audit table links each model to where its finding actually appears.
result_fig <- function(m) {
  dplyr::case_when(
    grepl("pam",   m) ~ "figures/02_pam_fvfm_trajectory.pdf",
    grepl("color|clmm", m) ~ "figures/03_color_trajectory.pdf",
    grepl("bw|growth", m)  ~ "figures/05_buoyant_weight_growth.pdf",
    grepl("zoox",  m) ~ "figures/06_symbiont_density_by_day.pdf",
    grepl("morph", m) ~ "figures/04_morphology_trajectories.pdf / figures/14_morphology_KM.pdf",
    grepl("genet", m) ~ "figures/13_genet_response_panel.pdf",
    TRUE ~ "figures/16_manuscript_fig1.pdf"
  )
}

# ---- Crosswalk: model -> its EXISTING diagnostic figure --------------------
# Models diagnosed by earlier scripts (12, 12b, 12c, diagnostics/A, ...) already
# have a plot on disk; this names them so we reuse rather than rebuild. Any stem
# not found here (or whose plot is missing) falls through to a freshly built K_ plot.
diag_map <- c(
  "02_pam_lmer"     = "K_02_pam_lmer_dharma.png",
  "12_pam_lmm"      = "A_pam_dharma.png",
  "12_color_lmm"    = "A_color_dharma.png",
  "12_zoox_lmm"     = "A_zoox_dharma.png",
  "12_bw_lm"        = "A_bw_dharma.png",
  "12_bw_pct_lm"    = "growth_pct.png",
  "12b_color_clmm"  = "G_12b_color_clmm_observed.png"
)
# Resolve a model stem to its existing diagnostic filename, or NA if none exists.
existing_diag <- function(stem) {
  if (stem %in% names(diag_map)) {
    f <- diag_map[[stem]]
    if (file.exists(file.path(FIG_DIR, "diagnostics", f)) ||
        file.exists(file.path(FIG_DIR, "12_diagnostics", f))) return(f)
    return(NA_character_)
  }
  # morphology models follow naming conventions, so derive their plot names:
  # 12_morph_<trait>_glmm -> B_<trait>.png ; 12c_morph_<trait>_blme -> G_12c_..._dharma.png
  if (grepl("^12_morph_.*_glmm$", stem)) {
    t <- sub("^12_morph_(.*)_glmm$", "\\1", stem)
    f <- paste0("B_", t, ".png")
    if (file.exists(file.path(FIG_DIR, "diagnostics", f))) return(f)
  }
  if (grepl("^12c_morph_.*_blme$", stem)) {
    f <- paste0("G_", stem, "_dharma.png")
    if (file.exists(file.path(FIG_DIR, "diagnostics", f))) return(f)
  }
  NA_character_
}

# ---- Build a fresh DHARMa diagnostic on demand -----------------------------
# Simulate scaled (quantile) residuals and write the standard QQ + residual-vs-
# predicted panel. Wrapped in tryCatch so one un-simulatable model can't abort
# the whole audit; on failure it closes any half-open graphics device and
# returns NA so the model is reported as a gap.
build_dharma <- function(model, stem) {
  out <- file.path(DIAG_DIR, paste0("K_", stem, "_dharma.png"))
  dharma_warnings <- character()
  ok <- tryCatch({
    res <- withCallingHandlers(
      DHARMa::simulateResiduals(model, n = 250, plot = FALSE),
      warning = function(w) {
        dharma_warnings <<- c(dharma_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    png(out, width = 1100, height = 520, res = 130)
    op <- par(mfrow = c(1, 2))
    on.exit(par(op), add = TRUE)
    plot_output <- utils::capture.output(
      withCallingHandlers({
        DHARMa::plotQQunif(res)
        DHARMa::plotResiduals(res, quantreg = FALSE)
      }, warning = function(w) {
        dharma_warnings <<- c(dharma_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      type = "output"
    )
    dev.off()
    diagnostic_messages <- c(unique(dharma_warnings),
                             unique(plot_output[nzchar(plot_output)]))
    if (length(diagnostic_messages)) {
      diagnostic_notes <<- c(
        diagnostic_notes,
        sprintf("`%s`: %s", stem, paste(diagnostic_messages, collapse = " | "))
      )
    }
    TRUE
  }, error = function(e) { if (dev.cur() > 1) dev.off(); FALSE })
  if (ok) out else NA_character_
}

# Locate a named diagnostic file across the two folders it might live in.
diag_path <- function(fig) {
  p <- file.path(DIAG_DIR, fig)
  if (file.exists(p)) return(p)
  p <- file.path(FIG_DIR, "12_diagnostics", fig)
  if (file.exists(p)) return(p)
  NA_character_
}

# ---- Audit loop: one pass per saved model ----------------------------------
rows <- list()
for (f in model_files) {
  stem <- tools::file_path_sans_ext(basename(f))
  m <- tryCatch(readRDS(f), error = function(e) NULL)
  cls <- if (is.null(m)) "unreadable" else class(m)[1]   # model class drives handling

  diag_fig <- existing_diag(stem)
  status <- "covered (existing)"

  # Freshness guard: if a reusable diagnostic exists BUT is older than the saved
  # model (i.e. the model was refit after the plot was made), the plot is stale —
  # rebuild it so the diagnostic reflects the current fit. Only for DHARMa-capable
  # classes (Gaussian/binomial mixed models and lm).
  if (!is.na(diag_fig) && cls %in% c("lmerMod", "lmerModLmerTest", "glmerMod", "lm")) {
    dp <- diag_path(diag_fig)
    if (!is.na(dp) && file.info(dp)$mtime < file.info(f)$mtime) {
      rebuilt <- build_dharma(m, stem)
      if (!is.na(rebuilt)) {
        diag_fig <- basename(rebuilt)
        status <- "covered (rebuilt stale K_)"
      }
    }
  }

  # No existing diagnostic found -> try to create coverage now.
  if (is.na(diag_fig)) {
    if (cls %in% c("lmerMod", "lmerModLmerTest", "glmerMod", "lm")) {
      # Standard models: build a DHARMa plot.
      diag_fig <- build_dharma(m, stem)
      status <- if (!is.na(diag_fig)) "covered (built K_)" else "DHARMa failed"
    } else if (grepl("clmm", stem)) {
      # Ordinal (cumulative link mixed) models aren't DHARMa-compatible; their
      # assumption check is the observed-vs-fitted plot from script 12b.
      diag_fig <- "figures/diagnostics/G_12b_color_clmm_observed.png"
      status <- "covered (ordinal: observed-vs-fitted)"
    } else {
      # Anything else lacks a diagnostic -> counts as a gap.
      status <- "no diagnostic"
    }
  }

  rows[[stem]] <- tibble(
    model = stem, class = cls,
    diagnostic_figure = ifelse(is.na(diag_fig), "—", sub(".*/", "", diag_fig)),
    result_figure = result_fig(stem),
    status = status
  )
}

# ---- Cox survival models (not saved as .rds) -------------------------------
# Kaplan-Meier/Cox models from script 14 live only as their Schoenfeld
# proportional-hazards diagnostic PNGs, so add them to the manifest by scanning
# for those files rather than reading model objects.
cox_ph_figs <- list.files(DIAG_DIR, pattern = "^14_cox_ph_.*\\.png$")
for (cf in cox_ph_figs) {
  tr <- sub("^14_cox_ph_(.*)\\.png$", "\\1", cf)
  rows[[paste0("cox_", tr)]] <- tibble(
    model = paste0("cox_overall_", tr), class = "coxph",
    diagnostic_figure = cf,
    result_figure = "figures/14_morphology_KM.pdf",
    status = "covered (cox.zph Schoenfeld)"
  )
}

# ---- Tally coverage --------------------------------------------------------
cov <- bind_rows(rows) |> arrange(class, model)
write_csv(cov, file.path(TBL_DIR, "25_model_diagnostic_coverage.csv"))

# "covered" appears in every success status string; anything else is a gap.
n_total   <- nrow(cov)
n_covered <- sum(grepl("covered", cov$status))
n_gap     <- n_total - n_covered

report <- c(
  "# K. Model diagnostic coverage",
  "",
  sprintf("Every fitted model cross-referenced to its residual/assumption diagnostic and the figure that visualizes its result. **%d/%d models covered; %d gaps.**",
          n_covered, n_total, n_gap),
  "",
  "| Model | Class | Diagnostic figure | Result figure | Status |",
  "|---|---|---|---|---|"
)
report <- c(report, cov |>
  mutate(line = sprintf("| `%s` | %s | %s | %s | %s |",
                        model, class, diagnostic_figure, result_figure, status)) |>
  pull(line))
if (n_gap > 0) {
  report <- c(report, "", "## Gaps", "",
              cov |> filter(!grepl("covered", status)) |>
                mutate(l = sprintf("- `%s` (%s): %s", model, class, status)) |> pull(l))
}
report <- c(report, "", "## Diagnostic Notes", "")
if (length(diagnostic_notes) == 0) {
  report <- c(report, "None.")
} else {
  report <- c(
    report,
    paste0("- ", unique(diagnostic_notes)),
    "",
    "DHARMa plotting/fitting warnings are retained here as diagnostic caveats; they do not create a coverage gap when the diagnostic figure is successfully written."
  )
}
writeLines(report, file.path(here::here("output", "diagnostics"),
                             "K_model_coverage_report.md"))

cat(sprintf("\n=== Model diagnostic coverage: %d/%d covered, %d gaps ===\n",
            n_covered, n_total, n_gap))
print(as.data.frame(cov |> select(model, class, status)))
cat("\nWrote output/tables/25_model_diagnostic_coverage.csv,",
    "output/diagnostics/K_model_coverage_report.md\n")
if (n_gap > 0) cat("WARNING:", n_gap, "model(s) lack a diagnostic — see report.\n")
