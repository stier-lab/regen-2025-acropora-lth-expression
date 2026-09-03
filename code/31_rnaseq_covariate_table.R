# =============================================================================
# Purpose: HANDOFF CONVENIENCE for the RNA-seq analysis (Shreya / Bay lab) —
#          build one tidy table with a row per RNA-seq library, joined to the
#          phenotype covariates, so the gene-expression data can be modeled
#          against organismal traits without re-deriving the joins.
#
#          Coding is harmonized to the phenotype convention (data dictionary):
#          wound U/W -> no/yes; source patch A/C/D -> lowercase a/c/d;
#          Temp_C -> 28C/31C.
#          Join key: Fragment_ID == id.
#
#          Covariate types:
#            - design        : treatment, wound, source patch, tank, day, plate
#                              (per library)
#            - per-fragment   : symbiont density at that biopsy (where measured)
#            - per-source     : resilience score / PCA displacement / rank
#                               (source-level; identical for all libraries from
#                               a source patch)
#          (Per-coral regeneration outcome is NOT attached: RNA-seq fragments were
#          destructively biopsied at D1/D3/D10, before the regeneration trajectory.)
#
# What & why: a coral fragment that was sequenced is the SAME individual we
#   measured phenotypes on. Asking whether gene expression tracks the organism's
#   physiology requires expression and phenotype side by side, keyed to the
#   same fragment. This script does that bookkeeping ONCE: it takes the list of
#   RNA-seq libraries and glues on (a) the experimental design each fragment was
#   under, (b) its symbiont density at biopsy, and (c) its source-level resilience
#   scores. It writes TWO files — a harmonized join (recoded to match the rest of
#   the phenotype analysis) and a raw, un-recoded lookup (so the gene-expression side can
#   choose her own factor levels / reference categories). It fits no model; it
#   hands the expression analyst a ready-to-merge covariate table.
#
# Input:   data/raw/plate_layout/Selected_Samples.csv
#          data/processed/symbiont_chl_clean.rds
#          output/tables/19_genet_resilience_summary.csv
# Output:  output/tables/31_rnaseq_phenotype_covariates.csv      (harmonized)
#          output/tables/31_rnaseq_library_lookup_raw.csv        (raw, un-recoded)
# =============================================================================

# 00_setup.R loads packages and defines shared paths (DATA_RAW, DATA_PROC,
# TBL_DIR for output/tables, ...).
source(here::here("code", "00_setup.R"))

# ---- RNA-seq library rows, as they appear in the source plate layout ---------
# Selected_Samples.csv is the wet-lab plate map; it lists every well, of several
# sample types. We keep only SampleType == "gene" (the RNA-seq libraries) — the
# other rows (e.g. SNP/DNA samples) are not part of the expression covariate table.
sel_raw <- read_csv(file.path(DATA_RAW, "plate_layout", "Selected_Samples.csv"),
                    show_col_types = FALSE) |>
  filter(SampleType == "gene")

# ---- RAW, un-recoded lookup --------------------------------------------------
# library_id <-> fragment_id <-> the original design fields, verbatim from
# Selected_Samples.csv. We impose NO recoding here (no 28C/31C, no a/c/d, no
# no/yes) so Shreya can set her own factor levels and reference categories. The
# harmonized table below is the convenience join; this is the source of truth
# for "what was the raw design value for this library."
raw_lookup <- sel_raw |>
  transmute(library_id  = Sample_ID,
            fragment_id = Fragment_ID,
            Temp_C, Day, Wound, Tank, Genotype, Plate, WoundOrder)
write_csv(raw_lookup, file.path(TBL_DIR, "31_rnaseq_library_lookup_raw.csv"))

# ---- RNA-seq library list (144 'gene' libraries), harmonized to phenotype coding
# transmute() keeps ONLY the columns it names. Here we both rename and recode so
# the design fields read identically to the rest of the phenotype pipeline (so a
# later join aligns without further recoding). `id` is the integer Fragment_ID — the join key used
# below to attach per-fragment symbiont density.
libs <- sel_raw |>
  transmute(
    library_id = Sample_ID,
    id         = as.integer(Fragment_ID),
    treatment  = paste0(as.integer(Temp_C), "C"),     # 28C / 31C
    wound      = if_else(Wound == "W", "yes", "no"),  # U/W -> no/yes
    genet      = str_to_lower(Genotype),              # legacy name: source patch A/C/D -> a/c/d
    tank       = as.integer(Tank),
    day        = as.integer(Day),
    plate      = as.integer(Plate)
  )

# ---- Per-fragment symbiont density (destructive biopsy; one row per coral) ---
# The symbiont table can have repeats per fragment; distinct(id) keeps one row
# per coral so the left_join below stays one-row-per-library (no fan-out).
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  distinct(id, .keep_all = TRUE) |>
  transmute(id, symbiont_cells_per_cm2 = cells_per_cm2)

# ---- Per-source resilience covariates (from the dashboard) ------------------
# These three numbers describe the field source patch (A/C/D), not a verified
# genetic individual, so every library sharing a source label gets the same value
# once joined by `genet`. The column name is preserved for backward compatibility
# with earlier handoff files.
res <- read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"),
                show_col_types = FALSE) |>
  transmute(genet                       = thicket,
            genet_mean_heat_sensitivity = round(mean_sensitivity, 3),
            genet_pca_displacement      = round(pca_displacement, 3),
            genet_resilience_rank       = rank_overall)

# Assemble: start from one row per library, attach per-fragment symbiont density
# (by id) and per-source resilience (by the legacy `genet` column). left_join keeps every library even
# if a covariate is missing (NA) — e.g. fragments with no symbiont measurement.
covariates <- libs |>
  left_join(sym, by = "id") |>
  left_join(res, by = "genet") |>
  arrange(plate, library_id)

write_csv(covariates, file.path(TBL_DIR, "31_rnaseq_phenotype_covariates.csv"))

cat("\n=== RNA-seq phenotype covariate table ===\n")
cat(sprintf("Libraries: %d  (design cells: %d temp x %d wound x %d source x %d day)\n",
            nrow(covariates),
            n_distinct(covariates$treatment), n_distinct(covariates$wound),
            n_distinct(covariates$genet), n_distinct(covariates$day)))
cat(sprintf("With a symbiont-density value: %d / %d\n",
            sum(!is.na(covariates$symbiont_cells_per_cm2)), nrow(covariates)))
cat("Wrote output/tables/31_rnaseq_phenotype_covariates.csv (harmonized)\n")
cat("Wrote output/tables/31_rnaseq_library_lookup_raw.csv (raw, un-recoded)\n")
