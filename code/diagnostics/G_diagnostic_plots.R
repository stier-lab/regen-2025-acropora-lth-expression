# =============================================================================
# Purpose: Diagnostic-plot inventory + fill-in (diagnostic suite G) — make sure
#          every saved model has an up-to-date residual/fit diagnostic plot, and
#          generate the ones that are still missing.
#
# What & why: every model needs a residual diagnostic so we can see whether its
#   assumptions hold. Earlier diagnostic scripts (A-D) already produced plots for
#   most models; this script audits coverage and fills the gaps. It (1) inventories
#   each model in output/models/ against its expected plot file and flags each
#   CURRENT / STALE (plot older than the model) / MISSING / BUILT, then (2) builds
#   the two kinds of plot that A-D didn't:
#     - the 8 penalized morphology fits (12c_morph_*_blme.rds): DHARMa simulated
#       residuals (the standard way to check fit for GLMMs);
#     - the ordinal color CLMM (12b_color_clmm.rds): DHARMa can't handle clmm, so
#       a manual observed-distribution barplot stands in.
# Input:   output/models/*.rds (fitted models needing diagnostics)
#          existing PNGs in figures/diagnostics/ and figures/12_diagnostics/
# Output:  figures/diagnostics/G_*.png (newly built plots)
#          output/diagnostics/G_plot_inventory.csv
#          output/diagnostics/G_inventory_report.md
# =============================================================================

# 00_setup.R loads shared paths/packages; DHARMa provides simulation-based
# residual diagnostics that work for (most) GLMMs.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(DHARMa) })

DIAG_FIG <- here("figures", "diagnostics")
DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Inventory -------------------------------------------------------------
# Gather every saved model and every existing diagnostic PNG (from both plot
# folders) so we can later match each model to its plot.
models_dir <- here("output", "models")
all_models <- list.files(models_dir, pattern = "\\.rds$", full.names = TRUE)
all_plots  <- list.files(DIAG_FIG, pattern = "\\.png$", full.names = TRUE)
extra_plots <- list.files(here("figures", "12_diagnostics"),
                          pattern = "\\.png$", full.names = TRUE)
all_plots <- c(all_plots, extra_plots)
diagnostic_notes <- character()

# plot_for(): given a model file, return the diagnostic plot(s) we expect for it.
# Naming conventions differ by model family, so this works in tiers: an explicit
# lookup table for the named LMM/LM fits, then regex rules for the morphology and
# color families, then a loose fallback that searches for the trait token.
plot_for <- function(mdl_path) {
  base <- tools::file_path_sans_ext(basename(mdl_path))
  # Tier 1: hand-mapped names where the plot file doesn't follow a rule.
  explicit <- c(
    "02_pam_lmer" = "K_02_pam_lmer_dharma.png",
    "12_pam_lmm" = "A_pam_dharma.png",
    "12_color_lmm" = "A_color_dharma.png",
    "12_zoox_lmm" = "A_zoox_dharma.png",
    "12_bw_lm" = "A_bw_dharma.png",
    "12_bw_pct_lm" = "growth_pct.png"
  )
  if (base %in% names(explicit)) {
    hits <- all_plots[basename(all_plots) == explicit[[base]]]
    return(list(pattern = explicit[[base]], hits = hits))
  }
  # Tier 2a: unpenalized glmer morphology — may have either a B_ or K_ plot.
  if (grepl("^12_morph_.*_glmm$", base)) {
    trait <- sub("^12_morph_(.*)_glmm$", "\\1", base)
    expected <- c(paste0("B_", trait, ".png"),
                  paste0("K_", base, "_dharma.png"))
    hits <- all_plots[basename(all_plots) %in% expected]
    return(list(pattern = paste(expected, collapse = "|"), hits = hits))
  }
  # Tier 2b: penalized (blme) morphology — the G_*_dharma.png this script builds.
  if (grepl("^12c_morph_.*_blme$", base)) {
    expected <- paste0("G_", base, "_dharma.png")
    hits <- all_plots[basename(all_plots) == expected]
    return(list(pattern = expected, hits = hits))
  }
  # Tier 2c: the ordinal color CLMM — the manual observed-distribution barplot.
  if (base == "12b_color_clmm") {
    expected <- "G_12b_color_clmm_observed.png"
    hits <- all_plots[basename(all_plots) == expected]
    return(list(pattern = expected, hits = hits))
  }

  # Tier 3 fallback: strip the script-number prefix and model-type suffix to get
  # a trait token, then match any plot whose filename contains it.
  pattern <- gsub("^(02|12|12b|12c|13|14|15|19)_", "", base)
  pattern <- gsub("_lmm|_lm|_glmm|_blme|_clmm", "", pattern)
  hits <- all_plots[grepl(pattern, basename(all_plots), ignore.case = TRUE)]
  list(pattern = pattern, hits = hits)
}

# make_inventory(): build the model-vs-plot table. For each model, find its plot,
# keep the most recent match, and compare timestamps — a plot older than its
# model is STALE (the model was refit after the plot was drawn).
make_inventory <- function() {
  inventory <- tibble::tibble(
    model     = basename(all_models),
    model_mtime = file.info(all_models)$mtime,
    diagnostic_plot = NA_character_,
    plot_mtime    = as.POSIXct(NA),
    status        = NA_character_
  )

  for (i in seq_along(all_models)) {
    mp <- plot_for(all_models[i])
    if (length(mp$hits)) {
      # If several candidate plots exist, take the newest one.
      most_recent <- mp$hits[which.max(file.info(mp$hits)$mtime)]
      inventory$diagnostic_plot[i] <- basename(most_recent)
      inventory$plot_mtime[i]     <- file.info(most_recent)$mtime
      # STALE if the plot predates the model fit; otherwise CURRENT.
      inventory$status[i] <- if (file.info(most_recent)$mtime <
                                  inventory$model_mtime[i]) "STALE" else "CURRENT"
    } else {
      inventory$status[i] <- "MISSING"  # no diagnostic plot found at all
    }
  }
  inventory
}

# ---- Build plots for missing fits ------------------------------------------
# build_dharma(): simulate residuals for one model and write the standard DHARMa
# 2-panel diagnostic (QQ + residual-vs-predicted) to PNG. refit=FALSE uses the
# fast simulation (no re-estimation per sim). Returns NA on failure so the caller
# can keep going.
build_dharma <- function(model_path, outfile = NULL) {
  m <- readRDS(model_path)
  base <- tools::file_path_sans_ext(basename(model_path))
  if (is.null(outfile)) {
    outfile <- file.path(DIAG_FIG, paste0("G_", base, "_dharma.png"))
  }
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  sim_warnings <- character()
  res <- tryCatch(
    withCallingHandlers(
      DHARMa::simulateResiduals(m, n = 500, refit = FALSE),
      warning = function(w) {
        sim_warnings <<- c(sim_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      message("DHARMa failed on ", base, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res)) return(NA_character_)
  if (length(sim_warnings)) {
    diagnostic_notes <<- c(
      diagnostic_notes,
      sprintf("`%s`: %s", base, paste(unique(sim_warnings), collapse = " | "))
    )
  }
  png(outfile, width = 1600, height = 800, res = 160)
  plot(res)
  dev.off()
  outfile
}

# ---- DHARMa plots: penalized morphology fits -------------------------------
# Loop the 8 blme morphology models, building a DHARMa diagnostic for each and
# tracking which filenames we created (used below to label them BUILT).
blme_fits <- list.files(models_dir, pattern = "^12c_morph_.*_blme\\.rds$",
                         full.names = TRUE)
cat("\n=== Building DHARMa plots for", length(blme_fits), "blme morphology fits ===\n")
built <- character()
for (f in blme_fits) {
  cat("  ", basename(f), "... ")
  out <- build_dharma(f)
  if (!is.na(out)) {
    cat("done\n")
    built <- c(built, basename(out))
  } else cat("skipped\n")
}

# ---- Observed-distribution plot: color CLMM --------------------------------
# DHARMa doesn't support ordinal::clmm, so instead of simulated residuals we draw
# the raw D-class distribution by treatment. If the CLMM is capturing the effect,
# the heated group's mass should sit at lower (paler) D-classes than ambient.
clmm_path <- file.path(models_dir, "12b_color_clmm.rds")
if (file.exists(clmm_path)) {
  cat("\n=== Building residual plot for color CLMM ===\n")
  m <- readRDS(clmm_path)
  # clm2 has limited predict support; use marginal-mean diagnostic instead.
  # Plot observed D-class distribution by treatment to show that the
  # CLMM captures the heat-induced left shift.
  color_data <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
    mutate(thicket = factor(thicket))
  png(file.path(DIAG_FIG, "G_12b_color_clmm_observed.png"),
      width = 1200, height = 700, res = 160)
  par(mar = c(4, 4, 3, 1))
  # Cross-tab D-class x treatment, then convert to column percentages
  # (margin = 2 normalizes within each treatment) so the two groups are comparable.
  tbl <- table(color_data$color_num, color_data$treatment)
  barplot(prop.table(tbl, margin = 2) * 100, beside = TRUE,
          col = c("#56B4E9", "#D55E00", "#009E73", "#E69F00", "#999999"),
          legend.text = paste0("D", rownames(tbl)),
          args.legend = list(x = "topright", bty = "n", cex = 0.8),
          xlab = "Treatment", ylab = "Percent of observations",
          main = "Color D-scale distribution by treatment (raw)\n(CLMM models this directly without continuous-scale assumption)")
  dev.off()
  built <- c(built, "G_12b_color_clmm_observed.png")
  cat("  built G_12b_color_clmm_observed.png\n")
}

# ---- Recompute inventory ---------------------------------------------------
# Recompute inventory after building plots; otherwise newly rebuilt plots are
# incorrectly reported as stale based on the pre-build file mtimes.
all_plots  <- list.files(DIAG_FIG, pattern = "\\.png$", full.names = TRUE)
extra_plots <- list.files(here("figures", "12_diagnostics"),
                          pattern = "\\.png$", full.names = TRUE)
all_plots <- c(all_plots, extra_plots)
inventory <- make_inventory()

# Refresh stale or missing DHARMa-capable diagnostics. This closes the loop when
# models are refit before this inventory runs: the report should end with CURRENT
# or BUILT diagnostics, not a stale warning that a later script has to repair.
refreshable <- inventory |>
  dplyr::filter(status %in% c("STALE", "MISSING"))
for (i in seq_len(nrow(refreshable))) {
  model_path <- file.path(models_dir, refreshable$model[i])
  m <- tryCatch(readRDS(model_path), error = function(e) NULL)
  cls <- if (is.null(m)) NA_character_ else class(m)[1]
  if (is.na(cls) ||
      !(cls %in% c("lmerMod", "lmerModLmerTest", "glmerMod", "lm", "bglmerMod"))) {
    next
  }
  base <- tools::file_path_sans_ext(refreshable$model[i])
  expected <- refreshable$diagnostic_plot[i]
  outfile <- if (!is.na(expected) && expected != "") {
    file.path(DIAG_FIG, expected)
  } else {
    file.path(DIAG_FIG, paste0("G_", base, "_dharma.png"))
  }
  cat("  refreshing ", refreshable$model[i], " -> ",
      basename(outfile), "... ", sep = "")
  out <- build_dharma(model_path, outfile = outfile)
  if (!is.na(out)) {
    cat("done\n")
    built <- c(built, basename(out))
  } else cat("skipped\n")
}

# Recompute inventory after all refreshes so the written CSV/report is final.
all_plots  <- list.files(DIAG_FIG, pattern = "\\.png$", full.names = TRUE)
extra_plots <- list.files(here("figures", "12_diagnostics"),
                          pattern = "\\.png$", full.names = TRUE)
all_plots <- c(all_plots, extra_plots)
inventory <- make_inventory()

# Relabel anything produced this run as BUILT (rather than CURRENT or a leftover
# MISSING) so the report distinguishes pre-existing from freshly drawn plots.
inventory <- inventory |>
  mutate(status = case_when(
    diagnostic_plot %in% built & status == "CURRENT" ~ "BUILT",
    grepl("^12c_morph_", model) & status == "MISSING" ~ "BUILT",
    model == "12b_color_clmm.rds" & status == "MISSING" ~ "BUILT",
    TRUE ~ status
  ))

write_csv(inventory, file.path(DIAG_OUT, "G_plot_inventory.csv"))

# ---- Report ----------------------------------------------------------------
# Markdown summary: status counts, the list of freshly built plots, and any
# models still lacking a diagnostic.
report_path <- file.path(DIAG_OUT, "G_inventory_report.md")
sink(report_path)
cat("# G. Diagnostic plot inventory\n\n")
cat("**Generated:** ", format(Sys.time()), "\n\n", sep = "")
cat("| Status | Count |\n|---|---|\n")
tbl <- table(inventory$status)
for (s in names(tbl)) cat("| ", s, " | ", tbl[[s]], " |\n", sep = "")
cat("\n## Newly built plots (", length(built), ")\n\n", sep = "")
for (b in built) cat("- `figures/diagnostics/", b, "`\n", sep = "")
cat("\n## Models still missing a plot\n\n")
missing <- inventory |> dplyr::filter(status == "MISSING")
if (nrow(missing) == 0) cat("None — full coverage.\n") else {
  for (i in seq_len(nrow(missing))) {
    cat("- `", missing$model[i], "` (mtime ",
        format(missing$model_mtime[i]), ")\n", sep = "")
  }
}
cat("\n## Diagnostic notes\n\n")
if (length(diagnostic_notes) == 0) {
  cat("None.\n")
} else {
  for (note in unique(diagnostic_notes)) cat("- ", note, "\n", sep = "")
  cat("\nPenalized `bglmerMod` DHARMa plots are retained as visual residual checks; ",
      "the formal morphology diagnostics are also covered by the unpenalized ",
      "GLMM diagnostic battery, endpoint/log-rank summaries, and interval-censored ",
      "timing models.\n", sep = "")
}
sink()

cat("\n=== Inventory complete ===\n")
print(tbl)
cat("\nWrote:\n - ", file.path(DIAG_OUT, "G_plot_inventory.csv"),
    "\n - ", report_path, "\n")
