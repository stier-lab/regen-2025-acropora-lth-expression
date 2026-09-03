# =============================================================================
# Purpose: PRELIMINARY RNA-seq SNP cluster integration.
#
#          Rachael Bay sent PRELIM_LTH_genoclusters.csv on 2026-09-02 with
#          sample-level genetic PCs, coverage, and preliminary genotype clusters
#          for the 144 RNA-seq libraries. Her note explicitly says these are
#          preliminary clustering results from an incomplete SNP set, with some
#          low-coverage samples possibly assigned to their own cluster. This
#          script therefore treats the file as an exploratory covariate layer,
#          not as final genotype identity.
#
# What & why: code/31 builds the stable RNA-seq phenotype covariate handoff.
#   This script joins Rachael's preliminary SNP layer onto that table by coral
#   fragment ID, audits the join, summarizes how source-patch labels map to
#   SNP clusters, checks treatment/wound/day balance by cluster, and runs the
#   only direct sample-level phenotype check currently available for the RNA-seq
#   libraries: symbiont density from the same destructive biopsies. PAM, color,
#   growth, and morphology were measured on different fragments, so they are
#   flagged as not directly joinable to this SNP file.
#
# Input:   data/raw/rnaseq/PRELIM_LTH_genoclusters.csv
#          output/tables/31_rnaseq_phenotype_covariates.csv
#          output/tables/19_genet_resilience_summary.csv
# Output:  output/tables/32_prelim_snp_join_audit.csv
#          output/tables/32_prelim_snp_response_joinability.csv
#          output/tables/32_prelim_snp_rnaseq_covariates.csv
#          output/tables/32_prelim_snp_cluster_summary.csv
#          output/tables/32_prelim_snp_source_resilience_overlay.csv
#          output/tables/32_prelim_snp_design_balance.csv
#          output/tables/32_prelim_snp_design_association_tests.csv
#          output/tables/32_prelim_snp_symbiont_summaries.csv
#          output/tables/32_prelim_snp_symbiont_model_summary.csv
#          output/tables/32_prelim_snp_symbiont_nested_tests.csv
#          output/tables/32_prelim_snp_symbiont_contrasts.csv
#          figures/32_prelim_snp_structure.{pdf,png}
#          figures/32_prelim_snp_design_balance.{pdf,png}
#          figures/32_prelim_snp_symbiont_heat_effects.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

snp_path <- file.path(DATA_RAW, "rnaseq", "PRELIM_LTH_genoclusters.csv")
cov_path <- file.path(TBL_DIR, "31_rnaseq_phenotype_covariates.csv")

if (!file.exists(snp_path)) {
  cat("\n=== Preliminary SNP integration skipped ===\n")
  cat("Missing data/raw/rnaseq/PRELIM_LTH_genoclusters.csv\n")
} else {
  if (!file.exists(cov_path)) {
    source(here::here("code", "31_rnaseq_covariate_table.R"))
  }

  snp_raw <- read_csv(snp_path, show_col_types = FALSE) |>
    clean_names()

  required_cols <- c("sample", "id", "pc1", "pc2", "pc3", "pc4",
                     "thicket", "cov", "cluster")
  missing_cols <- setdiff(required_cols, names(snp_raw))
  if (length(missing_cols) > 0) {
    stop("Preliminary SNP file is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  snp <- snp_raw |>
    transmute(
      prelim_snp_sample_label = sample,
      prelim_snp_library_id   = id,
      id_from_sample          = as.integer(str_extract(sample, "^[0-9]+")),
      id_from_library         = as.integer(str_extract(id, "^[0-9]+")),
      prelim_snp_pc1          = pc1,
      prelim_snp_pc2          = pc2,
      prelim_snp_pc3          = pc3,
      prelim_snp_pc4          = pc4,
      prelim_snp_source_thicket = str_to_lower(thicket),
      prelim_snp_cov          = cov,
      prelim_snp_cov_million  = cov / 1e6,
      prelim_snp_cluster      = as.character(cluster)
    ) |>
    mutate(
      prelim_snp_coverage_band = case_when(
        prelim_snp_cov < 5e6  ~ "<5M",
        prelim_snp_cov < 10e6 ~ "5-10M",
        TRUE                  ~ ">=10M"
      ),
      prelim_snp_lowest_decile_cov =
        prelim_snp_cov <= quantile(prelim_snp_cov, 0.10, na.rm = TRUE)
    )

  id_mismatch <- snp |>
    filter(is.na(id_from_sample) | is.na(id_from_library) |
             id_from_sample != id_from_library)
  if (nrow(id_mismatch) > 0) {
    stop("SNP sample labels and library IDs disagree for fragment IDs.")
  }

  snp <- snp |>
    transmute(
      id = id_from_library,
      prelim_snp_sample_label,
      prelim_snp_library_id,
      prelim_snp_pc1,
      prelim_snp_pc2,
      prelim_snp_pc3,
      prelim_snp_pc4,
      prelim_snp_source_thicket,
      prelim_snp_cov,
      prelim_snp_cov_million,
      prelim_snp_coverage_band,
      prelim_snp_lowest_decile_cov,
      prelim_snp_cluster
    )

  covariates <- read_csv(cov_path, show_col_types = FALSE) |>
    rename(
      source_thicket = genet,
      source_mean_heat_sensitivity = genet_mean_heat_sensitivity,
      source_pca_displacement = genet_pca_displacement,
      source_resilience_rank = genet_resilience_rank
    )

  joined <- covariates |>
    left_join(snp, by = "id")

  if (any(is.na(joined$prelim_snp_cluster))) {
    stop("Some RNA-seq covariate rows did not receive a preliminary SNP cluster.")
  }
  if (any(joined$source_thicket != joined$prelim_snp_source_thicket)) {
    stop("Source patch label in the SNP file disagrees with code/31 covariates.")
  }

  cluster_n <- joined |>
    count(prelim_snp_cluster, name = "cluster_n")
  small_clusters <- cluster_n |>
    filter(cluster_n < 6) |>
    pull(prelim_snp_cluster)
  small_label <- paste0("small clusters (",
                        paste(sort(as.integer(small_clusters)), collapse = "/"),
                        ")")
  cluster_order <- cluster_n |>
    arrange(as.integer(prelim_snp_cluster)) |>
    mutate(label = if_else(prelim_snp_cluster %in% small_clusters,
                           small_label,
                           paste0("cluster ", prelim_snp_cluster))) |>
    pull(label) |>
    unique()

  joined <- joined |>
    mutate(
      prelim_snp_cluster_analysis = if_else(
        prelim_snp_cluster %in% small_clusters,
        small_label,
        paste0("cluster ", prelim_snp_cluster)
      ),
      prelim_snp_cluster_analysis =
        factor(prelim_snp_cluster_analysis, levels = cluster_order),
      source_label = factor(str_to_upper(source_thicket), levels = c("A", "D", "C")),
      treatment = factor(treatment, levels = c("28C", "31C")),
      wound = factor(wound, levels = c("no", "yes")),
      day_f = factor(day, levels = c(1, 3, 10))
    ) |>
    arrange(plate, library_id)

  write_csv(joined, file.path(TBL_DIR, "32_prelim_snp_rnaseq_covariates.csv"))

  # ---- Join and response availability audits -------------------------------
  phenotype_sources <- tibble(
    response_variable = c("symbiont density", "PAM Fv/Fm", "color-card score",
                          "skeletal growth", "gross morphology milestones",
                          "microscope morphology milestones"),
    processed_file = c("symbiont_chl_clean.rds", "pam_clean.rds",
                       "color_clean.rds", "buoyant_weight_clean.rds",
                       "physio_clean.rds", "microscope_physio_clean.rds")
  )

  response_joinability <- phenotype_sources |>
    rowwise() |>
    mutate(
      processed_path = file.path(DATA_PROC, processed_file),
      processed_exists = file.exists(processed_path),
      processed_unique_ids = if (processed_exists) {
        length(unique(readRDS(processed_path)$id))
      } else {
        NA_integer_
      },
      matched_rnaseq_ids = if (processed_exists) {
        sum(unique(joined$id) %in% unique(readRDS(processed_path)$id))
      } else {
        NA_integer_
      },
      direct_sample_level_join = matched_rnaseq_ids == nrow(joined),
      current_use = case_when(
        response_variable == "symbiont density" & direct_sample_level_join
          ~ "direct SNP-by-response check",
        matched_rnaseq_ids == 0
          ~ "not directly joinable; use source/treatment summaries unless these fragments are genotyped",
        TRUE
          ~ "partial overlap; audit before modeling"
      )
    ) |>
    ungroup() |>
    select(-processed_path)

  write_csv(response_joinability,
            file.path(TBL_DIR, "32_prelim_snp_response_joinability.csv"))

  join_audit <- tibble(
    check = c(
      "SNP rows",
      "Unique SNP fragment IDs",
      "SNP sample label agrees with SNP library ID",
      "RNA-seq covariate rows",
      "All SNP IDs present in code/31 covariates",
      "All code/31 covariate IDs present in SNP file",
      "Source patch labels agree after join",
      "Rows with same-fragment symbiont density"
    ),
    value = c(
      as.character(nrow(snp)),
      as.character(n_distinct(snp$id)),
      as.character(nrow(id_mismatch) == 0),
      as.character(nrow(covariates)),
      as.character(all(snp$id %in% covariates$id)),
      as.character(all(covariates$id %in% snp$id)),
      as.character(all(joined$source_thicket == joined$prelim_snp_source_thicket)),
      as.character(sum(!is.na(joined$symbiont_cells_per_cm2)))
    ),
    status = c(
      if_else(nrow(snp) == 144, "PASS", "WARN"),
      if_else(n_distinct(snp$id) == nrow(snp), "PASS", "FAIL"),
      if_else(nrow(id_mismatch) == 0, "PASS", "FAIL"),
      if_else(nrow(covariates) == 144, "PASS", "WARN"),
      if_else(all(snp$id %in% covariates$id), "PASS", "FAIL"),
      if_else(all(covariates$id %in% snp$id), "PASS", "FAIL"),
      if_else(all(joined$source_thicket == joined$prelim_snp_source_thicket),
              "PASS", "FAIL"),
      if_else(sum(!is.na(joined$symbiont_cells_per_cm2)) == nrow(joined),
              "PASS", "WARN")
    ),
    notes = c(
      "Expected 144 RNA-seq libraries.",
      "One SNP row per RNA-seq fragment.",
      "Fragment number parsed from both SNP sample fields.",
      "Expected code/31 RNA-seq handoff table.",
      "Join uses fragment ID, not the differing library label strings.",
      "Join uses fragment ID, not the differing library label strings.",
      "Preliminary SNP source-patch label matches the wet-lab source-patch label.",
      "This is the only direct same-fragment response variable currently available."
    )
  )

  write_csv(join_audit, file.path(TBL_DIR, "32_prelim_snp_join_audit.csv"))

  # ---- Cluster summaries and design balance --------------------------------
  source_cluster_counts <- joined |>
    count(source_thicket, prelim_snp_cluster, name = "n") |>
    group_by(source_thicket) |>
    mutate(pct_source = n / sum(n)) |>
    ungroup()

  shared_clusters <- source_cluster_counts |>
    group_by(prelim_snp_cluster) |>
    summarise(n_sources = n_distinct(source_thicket), .groups = "drop") |>
    filter(n_sources > 1) |>
    pull(prelim_snp_cluster)

  cluster_summary <- joined |>
    group_by(prelim_snp_cluster) |>
    summarise(
      n_samples = n(),
      source_patches = paste(sort(unique(str_to_upper(source_thicket))),
                             collapse = "/"),
      n_source_patches = n_distinct(source_thicket),
      median_cov_million = median(prelim_snp_cov_million, na.rm = TRUE),
      min_cov_million = min(prelim_snp_cov_million, na.rm = TRUE),
      max_cov_million = max(prelim_snp_cov_million, na.rm = TRUE),
      n_below_5m = sum(prelim_snp_cov < 5e6),
      n_below_10m = sum(prelim_snp_cov < 10e6),
      pct_below_10m = mean(prelim_snp_cov < 10e6),
      .groups = "drop"
    ) |>
    arrange(as.integer(prelim_snp_cluster))

  write_csv(cluster_summary,
            file.path(TBL_DIR, "32_prelim_snp_cluster_summary.csv"))

  source_resilience_overlay <- joined |>
    mutate(shared_cluster = prelim_snp_cluster %in% shared_clusters) |>
    group_by(source_thicket) |>
    summarise(
      n_rnaseq_libraries = n(),
      n_prelim_snp_clusters = n_distinct(prelim_snp_cluster),
      dominant_cluster = names(which.max(table(prelim_snp_cluster))),
      dominant_cluster_n = max(as.integer(table(prelim_snp_cluster))),
      dominant_cluster_pct = dominant_cluster_n / n(),
      n_libraries_in_shared_ad_clusters = sum(shared_cluster),
      pct_libraries_in_shared_ad_clusters = mean(shared_cluster),
      median_cov_million = median(prelim_snp_cov_million, na.rm = TRUE),
      n_below_10m = sum(prelim_snp_cov < 10e6),
      .groups = "drop"
    ) |>
    left_join(
      read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"),
               show_col_types = FALSE) |>
        transmute(source_thicket = thicket,
                  source_mean_heat_sensitivity = mean_sensitivity,
                  source_pca_displacement = pca_displacement,
                  source_resilience_rank = rank_overall),
      by = "source_thicket"
    ) |>
    arrange(source_resilience_rank)

  write_csv(source_resilience_overlay,
            file.path(TBL_DIR, "32_prelim_snp_source_resilience_overlay.csv"))

  design_balance <- joined |>
    mutate(across(c(treatment, wound, day_f, tank, plate), as.character)) |>
    select(source_thicket, prelim_snp_cluster, prelim_snp_cluster_analysis,
           treatment, wound, day = day_f, tank, plate) |>
    pivot_longer(c(treatment, wound, day, tank, plate),
                 names_to = "design_factor", values_to = "level") |>
    count(source_thicket, prelim_snp_cluster, prelim_snp_cluster_analysis,
          design_factor, level, name = "n") |>
    arrange(design_factor, prelim_snp_cluster_analysis, level, source_thicket)

  write_csv(design_balance,
            file.path(TBL_DIR, "32_prelim_snp_design_balance.csv"))

  assoc_test <- function(dat, variable, scope) {
    cluster_group <- droplevels(dat$prelim_snp_cluster_analysis)
    design_level <- droplevels(as.factor(dat[[variable]]))
    tab <- table(cluster_group, design_level)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(
        scope = scope,
        design_factor = variable,
        method = "not estimable",
        statistic = NA_real_,
        p.value = NA_real_,
        n = nrow(dat),
        note = "Fewer than two cluster groups or fewer than two design levels."
      ))
    }

    tst <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 9999))
    tibble(
      scope = scope,
      design_factor = variable,
      method = "Pearson chi-square with simulated p-value (9999 replicates)",
      statistic = unname(tst$statistic),
      p.value = tst$p.value,
      n = nrow(dat),
      note = "Monte Carlo chi-square check; descriptive because clusters are preliminary and post hoc."
    )
  }

  design_tests <- bind_rows(
    lapply(c("treatment", "wound", "day_f", "plate"),
           \(v) assoc_test(joined, v, "all sources")),
    lapply(c("treatment", "wound", "day_f", "plate"),
           \(v) assoc_test(filter(joined, source_thicket == "a"),
                           v, "source A only")),
    lapply(c("treatment", "wound", "day_f", "plate"),
           \(v) assoc_test(filter(joined, source_thicket == "d"),
                           v, "source D only"))
  ) |>
    mutate(design_factor = dplyr::recode(design_factor, day_f = "day"))

  write_csv(design_tests,
            file.path(TBL_DIR, "32_prelim_snp_design_association_tests.csv"))

  # ---- Direct sample-level response check: symbiont density -----------------
  sym <- joined |>
    filter(!is.na(symbiont_cells_per_cm2), symbiont_cells_per_cm2 > 0) |>
    mutate(
      symbiont_million_per_cm2 = symbiont_cells_per_cm2 / 1e6,
      log10_symbiont_cells_per_cm2 = log10(symbiont_cells_per_cm2)
    )

  sym_summaries <- sym |>
    group_by(source_thicket, prelim_snp_cluster_analysis, treatment, wound, day_f) |>
    summarise(
      n = n(),
      mean_symbiont_million_per_cm2 =
        mean(symbiont_million_per_cm2, na.rm = TRUE),
      se_symbiont_million_per_cm2 =
        sd(symbiont_million_per_cm2, na.rm = TRUE) / sqrt(n()),
      mean_log10_symbiont_cells_per_cm2 =
        mean(log10_symbiont_cells_per_cm2, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(source_thicket, prelim_snp_cluster_analysis, treatment, wound, day_f)

  write_csv(sym_summaries,
            file.path(TBL_DIR, "32_prelim_snp_symbiont_summaries.csv"))

  m_design <- lm(log10_symbiont_cells_per_cm2 ~ treatment + wound + day_f,
                 data = sym)
  m_source <- lm(log10_symbiont_cells_per_cm2 ~
                   treatment * source_thicket + wound + day_f,
                 data = sym)
  m_cluster <- lm(log10_symbiont_cells_per_cm2 ~
                    treatment * prelim_snp_cluster_analysis + wound + day_f,
                  data = sym)
  m_pc <- lm(log10_symbiont_cells_per_cm2 ~
               treatment + wound + day_f +
               prelim_snp_pc1 + prelim_snp_pc2 +
               prelim_snp_pc3 + prelim_snp_pc4 +
               treatment:prelim_snp_pc1 + treatment:prelim_snp_pc2 +
               treatment:prelim_snp_pc3 + treatment:prelim_snp_pc4,
             data = sym)

  model_summary <- bind_rows(
    glance(m_design)  |> mutate(model = "design only"),
    glance(m_source)  |> mutate(model = "design + source_thicket + treatment:source_thicket"),
    glance(m_cluster) |> mutate(model = "design + prelim_snp_cluster + treatment:prelim_snp_cluster"),
    glance(m_pc)      |> mutate(model = "design + prelim_snp_PCs + treatment:prelim_snp_PCs")
  ) |>
    select(model, r.squared, adj.r.squared, sigma, statistic, p.value,
           df, logLik, AIC, BIC, deviance, df.residual, nobs)

  nested_tests <- bind_rows(
    tidy(anova(m_design, m_source)) |>
      mutate(comparison = "source_thicket terms beyond treatment/wound/day",
             step = row_number()),
    tidy(anova(m_design, m_cluster)) |>
      mutate(comparison = "preliminary SNP cluster terms beyond treatment/wound/day",
             step = row_number()),
    tidy(anova(m_design, m_pc)) |>
      mutate(comparison = "preliminary SNP PC terms beyond treatment/wound/day",
             step = row_number())
  ) |>
    select(comparison, step, everything())

  write_csv(model_summary,
            file.path(TBL_DIR, "32_prelim_snp_symbiont_model_summary.csv"))
  write_csv(nested_tests,
            file.path(TBL_DIR, "32_prelim_snp_symbiont_nested_tests.csv"))

  source_contrasts <- emmeans(m_source, ~ treatment | source_thicket) |>
    contrast(method = "pairwise") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    transmute(
      contrast_basis = "source_thicket",
      group = str_to_upper(source_thicket),
      contrast,
      estimate_log10_28c_minus_31c = estimate,
      SE, df,
      lower_CL = lower.CL,
      upper_CL = upper.CL,
      p.value,
      note = "Positive estimate means lower adjusted symbiont density at 31C."
    )

  cluster_contrasts <- emmeans(m_cluster,
                               ~ treatment | prelim_snp_cluster_analysis) |>
    contrast(method = "pairwise") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    transmute(
      contrast_basis = "prelim_snp_cluster",
      group = as.character(prelim_snp_cluster_analysis),
      contrast,
      estimate_log10_28c_minus_31c = estimate,
      SE, df,
      lower_CL = lower.CL,
      upper_CL = upper.CL,
      p.value,
      note = "Exploratory; small preliminary clusters have wide uncertainty."
    )

  symbiont_contrasts <- bind_rows(source_contrasts, cluster_contrasts)

  write_csv(symbiont_contrasts,
            file.path(TBL_DIR, "32_prelim_snp_symbiont_contrasts.csv"))

  # ---- Figures --------------------------------------------------------------
  coverage_shapes <- c("<5M" = 4, "5-10M" = 1, ">=10M" = 16)

  cluster_centroids <- joined |>
    group_by(prelim_snp_cluster) |>
    summarise(
      prelim_snp_pc1 = median(prelim_snp_pc1, na.rm = TRUE),
      prelim_snp_pc2 = median(prelim_snp_pc2, na.rm = TRUE),
      .groups = "drop"
    )

  p_pca <- ggplot(joined,
                  aes(prelim_snp_pc1, prelim_snp_pc2,
                      colour = source_thicket,
                      shape = prelim_snp_coverage_band)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.25) +
    geom_point(size = 2.5, alpha = 0.78, stroke = 0.7) +
    geom_label(data = cluster_centroids,
               aes(x = prelim_snp_pc1, y = prelim_snp_pc2,
                   label = prelim_snp_cluster),
               inherit.aes = FALSE,
               linewidth = 0.18, fill = "white", colour = "grey20",
               size = 3.0) +
    scale_colour_manual(values = PAL_GENO, name = "Source patch",
                        labels = c(a = "A", c = "C", d = "D")) +
    scale_shape_manual(values = coverage_shapes, name = "Coverage") +
    labs(x = "Preliminary SNP axis 1",
         y = "Preliminary SNP axis 2",
         title = "A. DNA-marker summary",
         subtitle = "Numbers mark cluster centers") +
    theme_pub(9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 8),
          legend.position = "bottom")

  p_cluster_counts <- source_cluster_counts |>
    mutate(prelim_snp_cluster =
             factor(prelim_snp_cluster,
                    levels = sort(unique(as.integer(prelim_snp_cluster))))) |>
    ggplot(aes(prelim_snp_cluster, n, fill = source_thicket)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.22) +
    geom_text(aes(label = n),
              position = position_stack(vjust = 0.5),
              size = 2.7, colour = "white") +
    scale_fill_manual(values = PAL_GENO, name = "Source patch",
                      labels = c(a = "A", c = "C", d = "D")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Preliminary SNP cluster",
         y = "RNA-seq libraries",
         title = "B. Cluster membership",
         subtitle = "C is coherent; A/D split and share clusters") +
    theme_pub(9) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 8),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom")

  p_structure <- p_pca + p_cluster_counts +
    plot_layout(widths = c(1.25, 1), guides = "collect") +
    plot_annotation(
      title = "Preliminary SNP clusters revise the source-patch interpretation",
      subtitle = str_wrap(
        "Rachael's 2026-09-02 file is preliminary. It supports source patch C as one coherent cluster, while source-patch labels A and D include multiple clusters and share some clusters.",
        width = 105
      ),
      caption = str_wrap(
        "Coverage bands are descriptive only; Rachael flagged low-coverage values as lower-confidence but did not provide a final quality-control threshold. These clusters should be replaced by the final SNP set when available.",
        width = 105
      )
    ) &
    theme_pub(9) &
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  save_fig(p_structure, "32_prelim_snp_structure", width = 190, height = 105)

  design_plot_data <- design_balance |>
    filter(design_factor %in% c("treatment", "wound", "day")) |>
    group_by(prelim_snp_cluster_analysis, design_factor, level) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(
      design_factor = factor(design_factor,
                             levels = c("treatment", "wound", "day"),
                             labels = c("Temperature", "Wound", "Day")),
      level = factor(level, levels = c("28C", "31C", "no", "yes", "1", "3", "10"))
    )

  p_design <- ggplot(design_plot_data,
                     aes(level, prelim_snp_cluster_analysis, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = n), size = 2.7, colour = "grey15") +
    facet_wrap(~ design_factor, scales = "free_x", nrow = 1) +
    scale_fill_gradient(low = "#F1F5F9", high = "#1F3864",
                        name = "n libraries") +
    labs(x = NULL,
         y = "Preliminary SNP cluster group",
         title = "Preliminary SNP clusters are after-the-fact design groups",
         subtitle = str_wrap(
           "The RNA-seq design is balanced by source patch, temperature, wound state, and day. SNP clusters split A and D after the fact, so cluster-level treatment cells are uneven and should be treated cautiously.",
           width = 105
         )) +
    theme_pub(9) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0),
          legend.position = "bottom")

  save_fig(p_design, "32_prelim_snp_design_balance", width = 175, height = 95)

  contrast_plot <- symbiont_contrasts |>
    filter(contrast == "28C - 31C") |>
    mutate(
      group = factor(group, levels = c("A", "D", "C",
                                       as.character(cluster_order))),
      contrast_basis = factor(contrast_basis,
                              levels = c("source_thicket",
                                         "prelim_snp_cluster"),
                              labels = c("Source patch",
                                         "Preliminary SNP cluster"))
    )

  p_sym <- ggplot(contrast_plot,
                  aes(estimate_log10_28c_minus_31c, group)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey55", linewidth = 0.35) +
    geom_segment(aes(x = lower_CL, xend = upper_CL, y = group, yend = group),
                 linewidth = 0.55, colour = "grey35") +
    geom_point(size = 2.7, colour = "#1F3864") +
    facet_grid(contrast_basis ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Heat effect on symbiont density\n(log scale; 28 °C - 31 °C)",
         y = NULL,
         title = "Direct DNA-marker check: symbiont density",
         subtitle = str_wrap(
           "Only symbiont density was measured on the same 144 RNA-seq fragments. Positive values mean fewer symbionts at 31 °C after comparing fragments with the same wound state and sampling day. Horizontal lines show exploratory uncertainty.",
           width = 105
         ),
         caption = str_wrap(
           "The comparison includes wound state and RNA-seq sampling day. Cluster estimates use preliminary after-the-fact groups; clusters with fewer than six RNA-seq libraries are combined.",
           width = 105
         )) +
    theme_pub(9) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
          strip.text.y = element_text(face = "bold"))

  save_fig(p_sym, "32_prelim_snp_symbiont_heat_effects",
           width = 175, height = 125)

  cat("\n=== Preliminary SNP integration ===\n")
  print(join_audit)
  cat("\n=== Preliminary SNP cluster summary ===\n")
  print(cluster_summary)
  cat("\n=== Direct response availability ===\n")
  print(response_joinability)
  cat("\nWrote preliminary SNP covariate tables and figures with prefix 32_prelim_snp_*\n")
}
