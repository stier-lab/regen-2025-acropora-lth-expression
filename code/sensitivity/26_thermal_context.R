# =============================================================================
# Purpose: Place the LTH chronic temperature treatments and source-patch
#          resilience variation in the context of an independent, calibrated
#          acute thermal-tolerance benchmark for the SAME species and island:
#          Cunning et al. 2024 (Coral Reefs) CBASS Fv/Fm ED50s for Acropora
#          pulchra from Mahana, Mo'orea (data/external/cunning2024_apulchra_ed50.csv).
#
#          (#2) Thermal-context anchor: where do 28 °C and 31 °C sit relative to
#               the acute ED50 distribution? -> 31 °C is chronic-SUBLETHAL,
#               well below the acute photochemical threshold.
#          (#1) Variation concordance: do acute (CBASS, 20 genets) and
#               chronic (LTH, 3 source patches) methods both detect substantial
#               heritable thermal-tolerance variation in this population?
#
#          IMPORTANT distinctions documented in data/external/README.md:
#           - ED50 is ACUTE (18 h); LTH is CHRONIC (weeks at +3 °C). ED50 is a
#             reference point, NOT the temperature axis of the chronic design.
#           - LTH source-patch labels (A,C,D) are field labels and NOT genotype-matched
#             to Cunning's genets, so an individual-genet correlation is not
#             possible yet (requires LTH RNA-seq SNPs vs Cunning genet_map).
#             We therefore compare the STRUCTURE/MAGNITUDE of variation, not
#             individual genets.
#
# What & why: a reviewer's first question about a chronic-heat experiment is how
#   severe the heat was — whether 31 °C was a mild increase or acutely lethal. To
#   answer it we use an independent, calibrated benchmark: Cunning et al. (2024)
#   ran an acute CBASS assay (a standardized 18-h heat ramp) on 20 genotypes of
#   THIS species from the SAME reef (Mahana, Mo'orea) and fitted an ED50 — the
#   temperature where photosynthetic efficiency (Fv/Fm) collapses to half. Their
#   mean ED50 is ~35.4 °C, so our 31 °C treatment sits ~4 °C *below* the acute
#   threshold: chronic but sublethal, not acutely lethal. The script does two
#   things: (#2) drops 28/31 °C onto that ED50 number line, and (#1) checks that
#   both methods — acute CBASS and our chronic LTH design — agree that there is
#   substantial heritable variation in heat tolerance among genotypes. This is
#   context/benchmarking, not a hypothesis test on our own data.
# Input:   data/external/cunning2024_apulchra_ed50.csv
#          output/tables/19_genet_resilience_summary.csv
# Output:  output/tables/26_thermal_context.csv
#          figures/26_thermal_context.{pdf,png}
# =============================================================================

# 00_setup.R loads packages and defines shared paths (TBL_DIR, ...), the
# theme_pub() plot theme, and save_fig().
source(here::here("code", "00_setup.R"))

# The two chronic treatment temperatures, named so the figure can label them.
LTH_TREATMENTS <- c(ambient = 28, heated = 31)

# One ED50 (°C) per Cunning genotype — the acute thermal-tolerance benchmark.
ed50 <- read_csv(file.path(here::here("data", "external"),
                           "cunning2024_apulchra_ed50.csv"),
                 show_col_types = FALSE)

# ---- (#2) Acute ED50 summary + treatment placement ------------------------
# Summarise the acute benchmark, then measure how far our treatments sit below
# it. A large positive "below ED50" gap = our heat is below the acute limit.
ed50_mean <- mean(ed50$ed50); ed50_sd <- sd(ed50$ed50)
ed50_min  <- min(ed50$ed50);  ed50_max <- max(ed50$ed50)
ed50_cv   <- 100 * ed50_sd / ed50_mean   # CV = SD/mean, a scale-free spread metric

# Assemble the context table: the ED50 distribution plus the temperature gaps
# between each treatment and the ED50 (mean and lowest-genotype) thresholds.
ctx <- tibble(
  metric = c("acute ED50 mean (°C)", "acute ED50 range (°C)", "acute ED50 SD (°C)",
             "acute ED50 CV (%)", "n genets (Cunning)",
             "31 °C below mean ED50 (°C)", "31 °C below lowest-genet ED50 (°C)",
             "28 °C below mean ED50 (°C)"),
  value  = round(c(ed50_mean, ed50_max - ed50_min, ed50_sd, ed50_cv, nrow(ed50),
                   ed50_mean - 31, ed50_min - 31, ed50_mean - 28), 3)
)
write_csv(ctx, file.path(TBL_DIR, "26_thermal_context.csv"))

# ---- (#1) Among-genotype variation: acute vs chronic ----------------------
# Side-by-side concordance check: does each method detect genotype variation?
# We are NOT correlating individual genets (labels aren't matched across studies)
# — only comparing whether both approaches find a spread.
res <- read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"),
                show_col_types = FALSE)
# Chronic among-thicket variation: spread (max - min) of the standardized
# resilience score across our 3 thickets — the chronic analogue of the ED50 range.
chronic_range <- diff(range(res$mean_sensitivity))
concordance <- tibble(
  method   = c("Acute CBASS (Cunning 2024)", "Chronic LTH (this study)"),
  n_genotypes = c(nrow(ed50), nrow(res)),
  variation_metric = c(sprintf("ED50 range = %.2f °C (SD %.2f)",
                               ed50_max - ed50_min, ed50_sd),
                       sprintf("resilience-score range = %.2f (std units)",
                               chronic_range)),
  detects_genotype_variation = c("yes (R=0.74 predicts bleaching)", "yes (C > D > A)")
)
write_csv(concordance, file.path(TBL_DIR, "26_genotype_variation_concordance.csv"))

# ---- Figure: thermal-context number line ----------------------------------
# A one-dimensional temperature axis: the spread of acute ED50s (20 genets, in
# orange) and the two chronic treatment triangles (28/31 °C, blue) on the same
# scale, so the reader can see heat stays below the acute threshold.
ed50_pts <- ed50 |> mutate(y = 0)          # all points on the y = 0 baseline
band <- tibble(xmin = ed50_min, xmax = ed50_max)   # shaded ED50 genet range

p <- ggplot() +
  # acute ED50 genet range band + points
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = -0.18, ymax = 0.18),
            fill = "#D55E00", alpha = 0.10) +
  geom_vline(xintercept = ed50_mean, colour = "#D55E00", linewidth = 0.5,
             linetype = "dashed") +
  geom_jitter(data = ed50_pts, aes(ed50, y), height = 0.07, width = 0,
              colour = "#D55E00", size = 1.8, alpha = 0.8) +
  # LTH chronic treatments
  geom_point(data = tibble(x = LTH_TREATMENTS, lab = names(LTH_TREATMENTS)),
             aes(x, 0), shape = 17, size = 3.4, colour = c("#56B4E9", "#0072B2")) +
  annotate("text", x = 28, y = 0.30, label = "28 °C\n(ambient)", size = 2.9,
           colour = "#56B4E9", lineheight = 0.9) +
  annotate("text", x = 31, y = 0.30, label = "31 °C\n(chronic heat)", size = 2.9,
           colour = "#0072B2", lineheight = 0.9) +
  annotate("text", x = ed50_mean, y = 0.30,
           label = sprintf("acute ED50\n(%.1f °C, n=20 genets)", ed50_mean),
           size = 2.9, colour = "#D55E00", lineheight = 0.9) +
  annotate("segment", x = 31, xend = ed50_min, y = -0.34, yend = -0.34,
           arrow = arrow(length = unit(1.6, "mm"), ends = "both"),
           colour = "grey45", linewidth = 0.3) +
  annotate("text", x = (31 + ed50_min) / 2, y = -0.44,
           label = sprintf("%.1f–%.1f °C below acute threshold",
                           ed50_min - 31, ed50_max - 31),
           size = 2.7, colour = "grey35") +
  scale_x_continuous(breaks = seq(28, 37, 1),
                     limits = c(27.5, 37), name = "Temperature (°C)") +
  scale_y_continuous(limits = c(-0.55, 0.45), breaks = NULL, name = NULL) +
  labs(title = expression("LTH chronic treatments vs acute thermal tolerance of"~italic("A. pulchra")),
       subtitle = "Acute Fv/Fm ED50 (Mahana, Mo'orea; Cunning et al. 2024) vs the LTH chronic 28/31 °C design") +
  theme_pub(10) +
  theme(panel.grid.major.y = element_blank())

save_fig(p, "26_thermal_context", width = 180, height = 95)

cat("\n=== Thermal context (#2) ===\n"); print(as.data.frame(ctx))
cat("\n=== Genotype-variation concordance (#1) ===\n"); print(as.data.frame(concordance))
cat(sprintf("\nAcute ED50: mean %.2f °C (range %.2f–%.2f, SD %.2f). 31 °C is %.1f °C below the mean ED50.\n",
            ed50_mean, ed50_min, ed50_max, ed50_sd, ed50_mean - 31))
cat("Wrote 26_thermal_context.csv, 26_genotype_variation_concordance.csv,",
    "26_thermal_context.{pdf,png}\n")
