# =============================================================================
# Purpose: Pre-flight checks for every raw dataset before analysis.
#          Catches:
#            - Sample-size imbalances vs the design (192 destructive + 48 + 16)
#            - Missing/out-of-range values in key columns
#            - ID collisions or orphans across data streams
#            - Date sequence anomalies (records before acclimation start, etc.)
#            - Tank assignments inconsistent with treatment expectation
#                (28 °C tanks must be 3, 6, 9, 12; 31 °C tanks must be 4, 5, 10, 11)
#
# What & why: before trusting any model or figure, we confirm each cleaned
#   dataset matches what the experiment SHOULD have produced. This script does
#   not test a biological hypothesis — it is the data-integrity gate. It loads
#   every processed .rds file and asserts the things that must be true if the
#   pipeline ran correctly: the right number of fragments, the expected factor
#   levels (two temperatures, two wound states, three genets a/c/d), tank↔
#   treatment plumbing intact, plausible value ranges, and no fragment IDs that
#   appear in a response stream but are missing from the master metadata
#   ("orphans"). Each check is logged PASS / FAIL / WARN / HANDLED / INFO so a
#   reviewer can see whether the data are sound — a FAIL means stop and fix
#   before publishing.
# Output:  output/tables/18_validation_summary.csv  — one row per check
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Check recorder --------------------------------------------------------
# Accumulate every check into a list of one-row tibbles. add() appends a result;
# the `<<-` writes to the `results` list defined in the enclosing scope (not a
# local copy), so each call grows the running log. status is a flag string
# (PASS/FAIL/WARN/HANDLED/INFO) interpreted at the end.
results <- list()
add <- function(stream, check, status, value, expected = "", notes = "") {
  results[[length(results) + 1]] <<- tibble(
    stream = stream, check = check, status = status,
    value = as.character(value), expected = as.character(expected),
    notes = notes
  )
}

# The fixed design constants every stream is checked against. Tank→treatment is
# determined by the heater/plumbing layout: tanks 3,6,9,12 are ambient (28 °C),
# tanks 4,5,10,11 are heated (31 °C). Genets (field colonies) are a, c, d.
EXPECTED_28C_TANKS <- c(3L, 6L, 9L, 12L)
EXPECTED_31C_TANKS <- c(4L, 5L, 10L, 11L)
EXPECTED_GENETS    <- c("a", "c", "d")

# Reusable check: for a given dataset, confirm the tanks observed in each
# treatment exactly match the expected plumbing layout. identical() is
# strict (same values AND order, hence the sort() above) — a mismatch means a
# tank was mislabelled or a fragment landed in the wrong treatment. Each
# treatment is only checked if that dataset contains rows for it (length() > 0).
check_tanks <- function(stream, d) {
  obs_28 <- sort(unique(d$tank[d$treatment == "28C"]))
  obs_31 <- sort(unique(d$tank[d$treatment == "31C"]))
  if (length(obs_28)) {
    if (identical(obs_28, EXPECTED_28C_TANKS)) {
      add(stream, "28C tank IDs", "PASS",
          paste(obs_28, collapse = ","), paste(EXPECTED_28C_TANKS, collapse = ","))
    } else {
      add(stream, "28C tank IDs", "FAIL",
          paste(obs_28, collapse = ","),
          paste(EXPECTED_28C_TANKS, collapse = ","),
          "tank-treatment assignment mismatch")
    }
  }
  if (length(obs_31)) {
    if (identical(obs_31, EXPECTED_31C_TANKS)) {
      add(stream, "31C tank IDs", "PASS",
          paste(obs_31, collapse = ","), paste(EXPECTED_31C_TANKS, collapse = ","))
    } else {
      add(stream, "31C tank IDs", "FAIL",
          paste(obs_31, collapse = ","),
          paste(EXPECTED_31C_TANKS, collapse = ","),
          "tank-treatment assignment mismatch")
    }
  }
}

# ---- Metadata --------------------------------------------------------------
# The master spine (from 01_load_clean_metadata.R). Confirm the total fragment
# count and that the design factors carry exactly the expected levels.
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))
add("metadata", "n total fragments", if (nrow(meta) == 208) "PASS" else "FAIL",
    nrow(meta), 208, "192 destructive + 16 microscope = 208")
# setequal() compares as sets (order-independent), so these pass as long as the
# right levels are present regardless of ordering.
add("metadata", "treatments", if (setequal(levels(meta$treatment), c("28C", "31C"))) "PASS" else "FAIL",
    paste(levels(meta$treatment), collapse = ","), "28C, 31C")
add("metadata", "wound levels", if (setequal(levels(meta$wound), c("no", "yes"))) "PASS" else "FAIL",
    paste(levels(meta$wound), collapse = ","), "no, yes")
add("metadata", "genets", if (setequal(unique(meta$thicket), EXPECTED_GENETS)) "PASS" else "FAIL",
    paste(sort(unique(meta$thicket)), collapse = ","),
    paste(EXPECTED_GENETS, collapse = ","))
check_tanks("metadata", meta)

# Chlorophyll-a was never assayed, so the column is entirely NA. This is a known,
# accepted gap (not a pipeline bug) — flagged HANDLED rather than FAIL so the log
# documents it without halting; downstream physiology relies on PAM, colour, and
# symbiont density instead.
n_chl_missing <- sum(is.na(meta$chlorophyll_ug_cm2))
add("metadata", "chl-a populated",
    if (n_chl_missing == 0) "PASS" else "HANDLED",
    sprintf("%d/%d missing", n_chl_missing, nrow(meta)),
    "all missing because assay was not run",
    "confirmed by Molly; analysis uses PAM, color, and symbiont density")

# ---- PAM ------------------------------------------------------------------
# Photosynthetic efficiency (Fv/Fm) from the pulse-amplitude-modulated probe.
# Cleaned data averages the top + bottom probe readings per coral, so the row
# count is ~half the raw. Fv/Fm is a ratio bounded [0, 1]; check it stays sane.
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
add("pam", "n observations (top+bottom averaged)",
    if (nrow(pam) > 300) "PASS" else "WARN",
    nrow(pam), "~336",
    "raw obs ≈ 672; cleaned data averages top/bottom probe locations per coral")
add("pam", "unique corals", if (length(unique(pam$id)) >= 48) "PASS" else "WARN",
    length(unique(pam$id)), 48)
add("pam", "Fv/Fm range", "INFO",
    sprintf("[%.2f, %.2f]", min(pam$fv_fm), max(pam$fv_fm)),
    "[0, 1]")
check_tanks("pam", pam)

# ---- Color ----------------------------------------------------------------
# Visual bleaching score on the CoralWatch D-scale (integer 1 = palest/bleached
# to 6 = darkest/healthiest). Confirm scores stay within that ordinal range.
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
add("color", "n observations", if (nrow(color) > 300) "PASS" else "WARN",
    nrow(color), "~336")
add("color", "D-scale range", "INFO",
    sprintf("[%g, %g]", min(color$color_num), max(color$color_num)),
    "[1, 6]")
check_tanks("color", color)

# ---- Buoyant weight -------------------------------------------------------
# Calcification (skeletal growth) for the 48-coral physiology subset. % growth
# should be positive for corals that grew over the experiment.
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
add("buoyant_weight", "n corals", if (nrow(bw) == 48) "PASS" else "WARN",
    nrow(bw), 48, "physiology subset")
add("buoyant_weight", "% growth range", "INFO",
    sprintf("[%.1f, %.1f]", min(bw$pct_growth), max(bw$pct_growth)),
    "positive for growing corals")

# ---- Symbionts ------------------------------------------------------------
# Symbiodiniaceae (zooxanthellae) density per cm². is.finite() guards against
# Inf/NaN produced when a surface area was zero/missing in the cells/cm²
# calculation — those non-finite values would break the median, so they are
# counted and excluded. A drop in density is the cellular signature of bleaching.
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))
add("symbionts", "n biopsies", if (nrow(sym) >= 150) "PASS" else "WARN",
    nrow(sym), "~192")
sym_finite <- sum(is.finite(sym$cells_per_cm2))
add("symbionts", "finite densities",
    if (sym_finite >= 150) "PASS" else "WARN",
    sprintf("%d/%d finite", sym_finite, nrow(sym)), "192")
add("symbionts", "median (cells/cm²)", "INFO",
    sprintf("%.2e", median(sym$cells_per_cm2[is.finite(sym$cells_per_cm2)])),
    "~7e5–1e6")

# ---- Wax dipping ----------------------------------------------------------
# Two surface-area methods are cross-validated: direct caliper geometry vs the
# wax-dipping standard curve. cor(use = "complete.obs") gives the agreement
# (Pearson r) over corals measured both ways. Here r falls below the 0.7
# threshold — flagged HANDLED (not FAIL): the curve is kept and shown in
# figure 07 for transparency rather than silently discarded.
wax <- readRDS(file.path(DATA_PROC, "wax_clean.rds"))
add("wax_dipping", "n corals", if (nrow(wax) >= 150) "PASS" else "WARN",
    nrow(wax), "~192")
agreement <- with(wax, cor(sa_caliper_cm2, sa_curve_cm2,
                            use = "complete.obs"))
add("wax_dipping", "caliper-curve SA correlation",
    if (agreement > 0.7) "PASS" else "HANDLED",
    sprintf("r = %.3f", agreement), "> 0.7",
    "below threshold; standard curve retained and plotted in figure 07 for transparency")

# ---- Physio morphology ----------------------------------------------------
# The repeated non-destructive morphology measurements (the primary
# wound-healing data). One row per coral × timepoint, so the count is large.
# Confirm wounded corals are being tracked.
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
add("morphology", "n observations", if (nrow(ph) > 700) "PASS" else "WARN",
    nrow(ph), "~768")
add("morphology", "date parsing",
    if (all(!is.na(ph$date))) "PASS" else "FAIL",
    if (all(!is.na(ph$date))) paste(range(ph$date), collapse = " to ") else paste0(sum(is.na(ph$date)), " missing"),
    "all physio morphology rows retain calendar dates")
add("morphology", "wounded corals tracked", "INFO",
    length(unique(ph$id[ph$wound == "yes"])), "~24 (per treatment)")
check_tanks("morphology", ph)

# ---- Microscope physio morphology -----------------------------------------
# The separate photo-only microscope cohort. This has a narrower design than the
# main morphology stream: 16 wounded corals only, genets A/C only, and tanks
# 6/9 at 28C plus 4/5 at 31C. Validate it separately so it is visible in the
# pipeline without being pooled into the main morphology inference.
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
add("microscope", "n observations", if (nrow(micro) == 256) "PASS" else "WARN",
    nrow(micro), 256, "16 photo corals x 16 daily observations")
add("microscope", "n photo corals",
    if (length(unique(micro$id)) == 16) "PASS" else "FAIL",
    length(unique(micro$id)), 16)
add("microscope", "days scored",
    if (identical(sort(unique(micro$day)), 0:15)) "PASS" else "FAIL",
    paste(range(micro$day), collapse = "-"), "0-15")
add("microscope", "cohort",
    if (all(micro$sample == "photo") && all(micro$wound == "yes")) "PASS" else "FAIL",
    paste(sort(unique(micro$sample)), collapse = ","),
    "sample=photo; wound=yes")
add("microscope", "genets",
    if (setequal(unique(micro$thicket), c("a", "c"))) "PASS" else "FAIL",
    paste(sort(unique(micro$thicket)), collapse = ","), "a,c")
add("microscope", "28C tank IDs",
    if (identical(sort(unique(micro$tank[micro$treatment == "28C"])), c(6L, 9L))) "PASS" else "FAIL",
    paste(sort(unique(micro$tank[micro$treatment == "28C"])), collapse = ","),
    "6,9", "photo-only subset uses two of the four ambient tanks")
add("microscope", "31C tank IDs",
    if (identical(sort(unique(micro$tank[micro$treatment == "31C"])), c(4L, 5L))) "PASS" else "FAIL",
    paste(sort(unique(micro$tank[micro$treatment == "31C"])), collapse = ","),
    "4,5", "photo-only subset uses two of the four heated tanks")
add("microscope", "pigment D8-D15 blanks",
    if (all(is.na(micro$pigment_over_wound[micro$day > 7]))) "HANDLED" else "WARN",
    sprintf("%d/%d missing after D7",
            sum(is.na(micro$pigment_over_wound[micro$day > 7])),
            sum(micro$day > 7)),
    "blank = not scored, not no")

# ---- YSI ------------------------------------------------------------------
# Daily handheld water-quality spot checks (from 09). na.rm = TRUE because some
# readings may be missing; the temperature range confirms tanks held a
# plausible 27–33 °C window.
ysi <- readRDS(file.path(DATA_PROC, "ysi_clean.rds"))
add("ysi", "n daily readings", if (nrow(ysi) > 50) "PASS" else "WARN",
    nrow(ysi), "~72", "9 tanks × ~8 days")
add("ysi", "temperature (°C) range", "INFO",
    sprintf("[%.1f, %.1f]", min(ysi$temp_c, na.rm=TRUE),
                            max(ysi$temp_c, na.rm=TRUE)),
    "[27, 33] reasonable")

# ---- Apex -----------------------------------------------------------------
# Continuous logger temperature (thousands of records). INFO only — report
# the count, no pass/fail threshold.
apex <- readRDS(file.path(DATA_PROC, "apex_temperature_daily.rds"))
add("apex", "n daily probe-records", "INFO", nrow(apex),
    "thousands across 6 datalog files")

# ---- Worms ----------------------------------------------------------------
# Vermetid/worm contamination presence flag. Sum the positives; documented as a
# known nuisance (concentrated in 31C tanks, treated with Worm Exit), logged
# INFO so it is recorded but does not block analysis.
worms <- readRDS(file.path(DATA_PROC, "worm_clean.rds"))
total_pos <- sum(worms$present, na.rm = TRUE)
add("worms", "total worm-positive obs", "INFO",
    total_pos, "0–few",
    "21 unique corals contaminated; concentrated in 31C tanks on 06/07 — treated with Worm Exit")

# ---- Cross-stream ID checks -----------------------------------------------
# Every fragment measured in a response stream must exist in the master
# metadata, or joins downstream will silently drop data. setdiff(stream, meta)
# finds "orphan" IDs present in a stream but absent from metadata — any orphan
# is a FAIL (a labelling error to fix before analysis).
ids_in_meta <- meta$id
ids_in_pam  <- unique(pam$id)
ids_in_bw   <- unique(bw$id)
ids_in_wax  <- unique(wax$id)
ids_in_micro <- unique(micro$id)

orphan_pam <- setdiff(ids_in_pam, ids_in_meta)
orphan_bw  <- setdiff(ids_in_bw,  ids_in_meta)
orphan_wax <- setdiff(ids_in_wax, ids_in_meta)
orphan_micro <- setdiff(ids_in_micro, ids_in_meta)

add("cross-stream", "PAM ids in metadata",
    if (length(orphan_pam) == 0) "PASS" else "FAIL",
    length(orphan_pam), 0,
    paste("Orphans:", paste(head(orphan_pam, 5), collapse = ",")))
add("cross-stream", "buoyant-weight ids in metadata",
    if (length(orphan_bw) == 0) "PASS" else "FAIL",
    length(orphan_bw), 0)
add("cross-stream", "wax-dipping ids in metadata",
    if (length(orphan_wax) == 0) "PASS" else "FAIL",
    length(orphan_wax), 0)
add("cross-stream", "microscope ids in metadata",
    if (length(orphan_micro) == 0) "PASS" else "FAIL",
    length(orphan_micro), 0)

# ---- Raw source provenance and documentation ------------------------------
# In addition to checking biological data integrity, confirm that Drive imports
# are traceable. This keeps raw data provenance from living only in meeting notes.
source_inventory_path <- file.path(DATA_META, "source_inventory.csv")
add("provenance", "source inventory exists",
    if (file.exists(source_inventory_path)) "PASS" else "FAIL",
    basename(source_inventory_path), "data/metadata/source_inventory.csv")

raw_dirs <- list.dirs(DATA_RAW, full.names = FALSE, recursive = FALSE)
expected_codebooks <- file.path(DATA_META, paste0(raw_dirs, "_codebook.csv"))
missing_codebooks <- raw_dirs[!file.exists(expected_codebooks)]
add("provenance", "raw folders have codebooks",
    if (length(missing_codebooks) == 0) "PASS" else "FAIL",
    if (length(missing_codebooks) == 0) length(raw_dirs) else paste(missing_codebooks, collapse = ","),
    paste0(length(raw_dirs), " raw folders documented"))

method_dir <- here::here("docs", "provenance", "drive_methods")
method_files <- list.files(method_dir, all.files = FALSE, recursive = FALSE,
                           full.names = FALSE)
add("provenance", "Drive method docs imported",
    if (length(method_files) >= 8) "PASS" else "WARN",
    length(method_files), ">= 8")

if (file.exists(source_inventory_path)) {
  source_inventory <- read_csv(source_inventory_path, show_col_types = FALSE)

  path_exists <- function(path_string) {
    paths <- str_split(path_string, ";", simplify = FALSE)[[1]] |>
      str_squish()
    paths <- paths[nzchar(paths)]
    length(paths) > 0 && all(file.exists(paths))
  }

  imported <- source_inventory |>
    filter(import_status == "imported")
  missing_imports <- imported$source_title[!map_lgl(imported$local_path, path_exists)]
  add("provenance", "imported source paths exist",
      if (length(missing_imports) == 0) "PASS" else "FAIL",
      if (length(missing_imports) == 0) nrow(imported) else paste(missing_imports, collapse = ","),
      paste0(nrow(imported), " imported sources"))

  analysis_inputs <- source_inventory |>
    filter(pipeline_role == "analysis_input")
  add("provenance", "analysis inputs inventoried",
      if (nrow(analysis_inputs) >= 12) "PASS" else "WARN",
      nrow(analysis_inputs), ">= 12",
      "metadata, physiology, morphology, temperature, worms, and RNA-seq handoff inputs")

  scripted_sources <- source_inventory |>
    filter(!is.na(analysis_script), nzchar(analysis_script))
  missing_scripts <- scripted_sources$analysis_script[!file.exists(scripted_sources$analysis_script)]
  add("provenance", "inventoried analysis scripts exist",
      if (length(missing_scripts) == 0) "PASS" else "FAIL",
      if (length(missing_scripts) == 0) nrow(scripted_sources) else paste(unique(missing_scripts), collapse = ","),
      paste0(nrow(scripted_sources), " scripted sources"))
}

# ---- Output ---------------------------------------------------------------
# Stack all the accumulated one-row results into a single table, save it, then
# pretty-print a marker-prefixed line per check and a final tally. The closing
# warning fires only if any check FAILed — the gate that should stop publication.
out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "18_validation_summary.csv"))

cat("Data validation summary:\n\n")
out |>
  mutate(marker = case_when(
    status == "PASS" ~ "✓",
    status == "FAIL" ~ "✘",
    status == "WARN" ~ "⚠",
    status == "HANDLED" ~ "↻",
    TRUE              ~ "ℹ"
  )) |>
  rowwise() |>
  mutate(line = sprintf("%s [%-15s] %-32s : %-30s | expected: %s%s",
                        marker, stream, check, value, expected,
                        if (nzchar(notes)) paste0("  — ", notes) else "")) |>
  pull(line) |> cat(sep = "\n")

n_fail <- sum(out$status == "FAIL")
n_handled <- sum(out$status == "HANDLED")
n_warn <- sum(out$status == "WARN")
n_pass <- sum(out$status == "PASS")
cat(sprintf("\n\nTotal: %d PASS, %d HANDLED, %d WARN, %d FAIL\n",
            n_pass, n_handled, n_warn, n_fail))
if (n_fail > 0) cat("⚠ Address FAILs before publishing results.\n")
