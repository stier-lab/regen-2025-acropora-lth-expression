# =============================================================================
# Purpose: Formal trade-off tests and model diagnostics.
#
#          This script turns the exploratory trade-off plots into a reproducible
#          audit of what can and cannot be tested. It keeps three evidence levels
#          separate:
#            1) same-fragment growth/regeneration tests,
#            2) time-to-event regeneration tests with growth as a covariate,
#            3) source-patch-only heat-tolerance screens.
#
#          Source-patch screens are deliberately not treated as formal regression
#          tests because there are only three source patches. SNP trade-off tests
#          are also limited by joinability: current SNP calls join directly only
#          to symbiont density, not to the growth or morphology fragments.
#
# Input:   data/processed/buoyant_weight_clean.rds
#          data/processed/physio_clean.rds
#          output/tables/12_genet_treatment_effects.csv
#          output/tables/13_genet_emmeans.csv
#          output/tables/32_prelim_snp_response_joinability.csv
# Output:  output/tables/37_tradeoff_scope_inventory.csv
#          output/tables/37_tradeoff_join_audit.csv
#          output/tables/37_tradeoff_endpoint_results.csv
#          output/tables/37_tradeoff_robust_growth_results.csv
#          output/tables/37_tradeoff_discrete_transition_results.csv
#          output/tables/37_tradeoff_time_to_event_results.csv
#          output/tables/37_tradeoff_heat_penalty_leave_one_trait_out.csv
#          output/tables/37_tradeoff_source_level_screen.csv
#          output/tables/37_tradeoff_model_diagnostics.csv
#          output/tables/37_tradeoff_main_results.csv
#          figures/37_tradeoff_diagnostics/*.png
# =============================================================================

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages({
  library(survival)
  library(brglm2)
  library(robustbase)
  library(sandwich)
  library(lmtest)
})

DIAG_DIR <- file.path(FIG_DIR, "37_tradeoff_diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ---------------------------------------------------------------
clean_p <- function(x) {
  ifelse(is.na(x), NA_real_, as.numeric(x))
}

plain_status <- function(p, alpha = 0.05) {
  case_when(
    is.na(p) ~ "not available",
    p >= alpha ~ "passed",
    TRUE ~ "flagged"
  )
}

exact_perm_mean_diff <- function(y, event) {
  keep <- is.finite(y) & !is.na(event)
  y <- y[keep]
  event <- as.integer(event[keep])
  n <- length(y)
  n_event <- sum(event == 1)
  n_nonevent <- sum(event == 0)
  if (n_event == 0 || n_nonevent == 0) {
    return(tibble(
      raw_yes_minus_no = NA_real_,
      permutation_p = NA_real_,
      n_event = n_event,
      n_nonevent = n_nonevent
    ))
  }

  observed <- mean(y[event == 1]) - mean(y[event == 0])
  combos <- utils::combn(seq_len(n), n_event)
  diffs <- apply(combos, 2, function(idx) {
    mean(y[idx]) - mean(y[-idx])
  })

  tibble(
    raw_yes_minus_no = observed,
    permutation_p = mean(abs(diffs) >= abs(observed) - 1e-12),
    n_event = n_event,
    n_nonevent = n_nonevent
  )
}

safe_ci <- function(estimate, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  c(estimate - z * se, estimate + z * se)
}

fit_lm_growth_cost <- function(dat, event_var, event_label, model_label,
                               include_treatment = TRUE) {
  dat2 <- dat |>
    mutate(event_status = factor(.data[[event_var]],
                                 levels = c(0, 1),
                                 labels = c("no", "yes"))) |>
    filter(!is.na(event_status), is.finite(pct_growth))

  rhs <- if (include_treatment) {
    "treatment + thicket + event_status"
  } else {
    "thicket + event_status"
  }
  model <- lm(as.formula(paste("pct_growth ~", rhs)), data = dat2)

  emm <- emmeans(model, "event_status")
  contrast_tbl <- summary(
    contrast(emm, method = list("yes_minus_no" = c(-1, 1))),
    infer = TRUE
  ) |>
    as_tibble()

  aug <- broom::augment(model)
  shapiro_p <- if (nrow(aug) >= 3) {
    stats::shapiro.test(aug$.resid)$p.value
  } else {
    NA_real_
  }
  cooks <- stats::cooks.distance(model)
  cook_threshold <- 4 / nrow(dat2)

  p1 <- ggplot(aug, aes(.fitted, .resid)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "Fitted growth", y = "Residual",
         title = paste(model_label, "residuals")) +
    theme_pub(9)

  p2 <- ggplot(aug, aes(sample = .std.resid)) +
    stat_qq(size = 2, alpha = 0.8) +
    stat_qq_line(linewidth = 0.35, colour = "grey35") +
    labs(x = "Theoretical quantile", y = "Standardized residual",
         title = paste(model_label, "Q-Q")) +
    theme_pub(9)

  p3 <- tibble(obs = seq_along(cooks), cooks = as.numeric(cooks)) |>
    ggplot(aes(obs, cooks)) +
    geom_hline(yintercept = cook_threshold, linewidth = 0.3,
               linetype = "dashed", colour = "grey45") +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "Observation", y = "Cook's distance",
         title = paste(model_label, "influence")) +
    theme_pub(9)

  diag_plot <- p1 + p2 + p3 + plot_layout(nrow = 1)
  diag_path <- file.path(DIAG_DIR, paste0(model_label, "_lm_diagnostics.png"))
  ggsave(diag_path, diag_plot, width = 210, height = 72,
         units = "mm", dpi = 300)

  result <- contrast_tbl |>
    transmute(
      model = model_label,
      analysis_type = "growth endpoint",
      event = event_label,
      n = nrow(dat2),
      n_event = sum(dat2$event_status == "yes"),
      n_nonevent = sum(dat2$event_status == "no"),
      estimate = estimate,
      lower_ci = lower.CL,
      upper_ci = upper.CL,
      p_value = p.value,
      estimate_scale = "percentage-point growth difference: event yes minus no",
      tradeoff_direction = "negative estimate would support a growth cost of regeneration",
      diagnostic_figure = diag_path
    )

  diagnostics <- tibble(
    model = model_label,
    diagnostic = c("Shapiro-Wilk residual normality",
                   "maximum Cook's distance",
                   "observations above 4/n Cook's threshold"),
    value = c(shapiro_p, max(cooks, na.rm = TRUE),
              sum(cooks > cook_threshold, na.rm = TRUE)),
    p_value = c(shapiro_p, NA_real_, NA_real_),
    status = c(plain_status(shapiro_p),
               if_else(max(cooks, na.rm = TRUE) < 1, "passed", "flagged"),
               if_else(sum(cooks > cook_threshold, na.rm = TRUE) <= 1,
                       "passed", "flagged")),
    diagnostic_figure = diag_path
  )

  influence <- tibble(
    model = model_label,
    id = dat2$id,
    treatment = if ("treatment" %in% names(dat2)) as.character(dat2$treatment) else NA_character_,
    thicket = if ("thicket" %in% names(dat2)) as.character(dat2$thicket) else NA_character_,
    pct_growth = dat2$pct_growth,
    event_status = as.character(dat2$event_status),
    cooks_distance = as.numeric(cooks),
    cooks_threshold = cook_threshold,
    flagged = cooks > cook_threshold
  )

  model_formula <- formula(model)
  leave_one_out <- map_dfr(seq_len(nrow(dat2)), \(i) {
    loo_dat <- dat2[-i, ]
    if (n_distinct(loo_dat$event_status) < 2) {
      return(tibble(
        model = model_label,
        omitted_id = dat2$id[[i]],
        omitted_event_status = as.character(dat2$event_status[[i]]),
        omitted_growth = dat2$pct_growth[[i]],
        estimate = NA_real_,
        p_value = NA_real_,
        status = "not estimable after omission"
      ))
    }
    loo_model <- tryCatch(lm(model_formula, data = loo_dat),
                          error = function(e) NULL)
    if (is.null(loo_model)) {
      return(tibble(
        model = model_label,
        omitted_id = dat2$id[[i]],
        omitted_event_status = as.character(dat2$event_status[[i]]),
        omitted_growth = dat2$pct_growth[[i]],
        estimate = NA_real_,
        p_value = NA_real_,
        status = "model failed after omission"
      ))
    }
    loo_contrast <- summary(
      contrast(emmeans(loo_model, "event_status"),
               method = list("yes_minus_no" = c(-1, 1))),
      infer = TRUE
    ) |>
      as_tibble()
    tibble(
      model = model_label,
      omitted_id = dat2$id[[i]],
      omitted_event_status = as.character(dat2$event_status[[i]]),
      omitted_growth = dat2$pct_growth[[i]],
      estimate = loo_contrast$estimate[[1]],
      p_value = loo_contrast$p.value[[1]],
      status = "estimated"
    )
  })

  list(model = model, result = result, diagnostics = diagnostics,
       influence = influence, leave_one_out = leave_one_out)
}

fit_robust_growth_cost <- function(dat, event_var, event_label, model_label,
                                   include_treatment = TRUE) {
  dat2 <- dat |>
    mutate(
      event_bin = as.integer(.data[[event_var]] == 1),
      event_status = factor(event_bin, levels = c(0, 1),
                            labels = c("no", "yes"))
    ) |>
    filter(!is.na(event_bin), is.finite(pct_growth))

  rhs <- if (include_treatment) {
    "treatment + thicket + event_bin"
  } else {
    "thicket + event_bin"
  }
  model_formula <- as.formula(paste("pct_growth ~", rhs))
  model <- robustbase::lmrob(model_formula, data = dat2, setting = "KS2014")
  coefs <- summary(model)$coefficients
  event_row <- coefs["event_bin", ]
  df_res <- max(model$df.residual, 1)
  ci <- event_row[["Estimate"]] +
    c(-1, 1) * qt(0.975, df = df_res) * event_row[["Std. Error"]]

  weights_rob <- stats::weights(model, type = "robustness")
  low_weight_count <- sum(weights_rob < 0.5, na.rm = TRUE)
  model_converged <- isTRUE(model$converged)

  robust_aug <- tibble(
    obs = seq_along(weights_rob),
    fitted = fitted(model),
    residual = residuals(model),
    std_residual = as.numeric(scale(residuals(model))),
    weight = weights_rob
  )

  p1 <- ggplot(robust_aug, aes(fitted, residual)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 2, alpha = 0.8) +
    labs(x = "Fitted growth", y = "Robust residual",
         title = paste(model_label, "residuals")) +
    theme_pub(9)

  p2 <- ggplot(robust_aug, aes(sample = std_residual)) +
    stat_qq(size = 2, alpha = 0.8) +
    stat_qq_line(linewidth = 0.35, colour = "grey35") +
    labs(x = "Theoretical quantile", y = "Standardized robust residual",
         title = paste(model_label, "Q-Q")) +
    theme_pub(9)

  p3 <- ggplot(robust_aug, aes(obs, weight)) +
    geom_hline(yintercept = 0.5, linewidth = 0.3,
               linetype = "dashed", colour = "grey45") +
    geom_point(size = 2, alpha = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Observation", y = "Robustness weight",
         title = paste(model_label, "weights")) +
    theme_pub(9)

  diag_path <- file.path(DIAG_DIR, paste0(model_label, "_lmrob_diagnostics.png"))
  ggsave(diag_path, p1 + p2 + p3 + plot_layout(nrow = 1),
         width = 210, height = 72, units = "mm", dpi = 300)

  result <- tibble(
    model = model_label,
    analysis_type = "robust growth endpoint",
    event = event_label,
    n = nrow(dat2),
    n_event = sum(dat2$event_bin == 1),
    n_nonevent = sum(dat2$event_bin == 0),
    estimate = event_row[["Estimate"]],
    lower_ci = ci[[1]],
    upper_ci = ci[[2]],
    p_value = clean_p(event_row[["Pr(>|t|)"]]),
    estimate_scale = "robust percentage-point growth difference: event yes minus no",
    tradeoff_direction = "negative estimate would support a growth cost of regeneration",
    diagnostic_figure = diag_path
  )

  diagnostics <- tibble(
    model = model_label,
    diagnostic = c("robust LM convergence",
                   "minimum robustness weight",
                   "observations with robustness weight < 0.5"),
    value = c(as.numeric(model_converged),
              min(weights_rob, na.rm = TRUE),
              low_weight_count),
    p_value = NA_real_,
    status = c(if_else(model_converged, "passed", "flagged"),
               if_else(min(weights_rob, na.rm = TRUE) > 0.05,
                       "passed", "flagged"),
               if_else(low_weight_count <= 1, "passed", "flagged")),
    diagnostic_figure = diag_path
  )

  list(model = model, result = result, diagnostics = diagnostics)
}

fit_brglm_regeneration <- function(dat, event_var, event_label, model_label,
                                   include_treatment = TRUE) {
  dat2 <- dat |>
    filter(!is.na(.data[[event_var]]), is.finite(pct_growth)) |>
    mutate(z_growth = as.numeric(scale(pct_growth)))

  rhs <- if (include_treatment) {
    "z_growth + treatment + thicket"
  } else {
    "z_growth + thicket"
  }
  model <- glm(
    as.formula(paste(event_var, "~", rhs)),
    data = dat2,
    family = binomial,
    method = brglmFit,
    type = "AS_mean"
  )

  coefs <- summary(model)$coefficients
  growth_row <- coefs["z_growth", ]
  ci <- safe_ci(growth_row[["Estimate"]], growth_row[["Std. Error"]])
  model_converged <- isTRUE(model$converged)

  diag_path <- file.path(DIAG_DIR, paste0(model_label, "_brglm_dharma.png"))
  dharma <- tryCatch({
    set.seed(42)
    sim <- suppressWarnings(
      DHARMa::simulateResiduals(model, n = 1000, plot = FALSE)
    )
    png(diag_path, width = 1800, height = 1500, res = 220)
    suppressWarnings(plot(sim))
    dev.off()
    list(
      uniformity_p = clean_p(DHARMa::testUniformity(sim)$p.value),
      dispersion_p = clean_p(DHARMa::testDispersion(sim)$p.value),
      outlier_p = clean_p(DHARMa::testOutliers(sim)$p.value),
      ok = TRUE
    )
  }, error = function(e) {
    list(uniformity_p = NA_real_, dispersion_p = NA_real_,
         outlier_p = NA_real_, ok = FALSE)
  })

  result <- tibble(
    model = model_label,
    analysis_type = "regeneration endpoint",
    event = event_label,
    n = nrow(dat2),
    n_event = sum(dat2[[event_var]] == 1, na.rm = TRUE),
    n_nonevent = sum(dat2[[event_var]] == 0, na.rm = TRUE),
    estimate = exp(growth_row[["Estimate"]]),
    lower_ci = exp(ci[[1]]),
    upper_ci = exp(ci[[2]]),
    p_value = clean_p(growth_row[["Pr(>|z|)"]]),
    estimate_scale = "odds ratio for event per 1 SD higher growth",
    tradeoff_direction = "odds ratio below 1 would support slower/lower regeneration in faster-growing fragments",
    diagnostic_figure = diag_path
  )

  diagnostics <- tibble(
    model = model_label,
    diagnostic = c("bias-reduced GLM convergence",
                   "DHARMa uniformity screen",
                   "DHARMa dispersion screen",
                   "DHARMa outlier screen"),
    value = c(as.numeric(model_converged),
              dharma$uniformity_p, dharma$dispersion_p, dharma$outlier_p),
    p_value = c(NA_real_,
                dharma$uniformity_p, dharma$dispersion_p, dharma$outlier_p),
    status = c(if_else(model_converged, "passed", "flagged"),
               plain_status(dharma$uniformity_p),
               plain_status(dharma$dispersion_p),
               plain_status(dharma$outlier_p)),
    diagnostic_figure = diag_path
  )

  list(model = model, result = result, diagnostics = diagnostics)
}

first_event_table <- function(phys, bw, trait) {
  phys |>
    filter(wound == "yes") |>
    group_by(id, treatment, thicket, tank) |>
    summarise(
      max_day = max(day, na.rm = TRUE),
      event_day = {
        event_days <- day[!is.na(.data[[trait]]) & .data[[trait]] == 1]
        if (length(event_days) > 0) min(event_days, na.rm = TRUE) else NA_real_
      },
      event = !is.na(event_day),
      time = if_else(event, event_day, max_day),
      .groups = "drop"
    ) |>
    left_join(select(bw, id, pct_growth), by = "id") |>
    mutate(
      treatment = factor(treatment, levels = c("28C", "31C")),
      thicket = factor(thicket, levels = c("a", "c", "d")),
      z_growth = as.numeric(scale(pct_growth))
    )
}

fit_cox_timing <- function(phys, bw, trait, event_label) {
  dat <- first_event_table(phys, bw, trait)
  model_label <- paste0("cox_", trait, "_growth")

  model <- survival::coxph(
    survival::Surv(time, event) ~ treatment + z_growth + strata(thicket),
    data = dat,
    ties = "efron"
  )
  sm <- summary(model)
  coefs <- as_tibble(sm$coefficients, rownames = "term")
  conf <- as_tibble(sm$conf.int, rownames = "term")
  zph <- survival::cox.zph(model)
  zph_tbl <- as_tibble(zph$table, rownames = "term")

  diag_path <- file.path(DIAG_DIR, paste0(model_label, "_cox_zph.png"))
  schoenfeld <- residuals(model, type = "schoenfeld")
  if (is.null(dim(schoenfeld))) {
    schoenfeld_df <- tibble(
      event_time = as.numeric(names(schoenfeld)),
      term = names(coef(model))[1],
      residual = as.numeric(schoenfeld)
    )
  } else {
    schoenfeld_df <- as_tibble(schoenfeld, rownames = "event_time") |>
      mutate(event_time = as.numeric(event_time)) |>
      pivot_longer(-event_time, names_to = "term", values_to = "residual")
  }
  cox_diag_plot <- ggplot(schoenfeld_df, aes(event_time, residual)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~ term, scales = "free_y") +
    labs(x = "Event day",
         y = "Scaled Schoenfeld residual",
         title = paste(model_label, "proportional-hazards check")) +
    theme_pub(9)
  ggsave(diag_path, cox_diag_plot, width = 150, height = 90,
         units = "mm", dpi = 300)

  infl <- residuals(model, type = "dfbeta")
  max_dfbeta <- max(abs(infl), na.rm = TRUE)
  global_ph_p <- zph_tbl |>
    filter(term == "GLOBAL") |>
    pull(p)
  growth_ph_p <- zph_tbl |>
    filter(term == "z_growth") |>
    pull(p)

  result <- coefs |>
    filter(term == "z_growth") |>
    left_join(conf |> select(term, `lower .95`, `upper .95`), by = "term") |>
    transmute(
      model = model_label,
      analysis_type = "time to regeneration",
      event = event_label,
      n = nrow(dat),
      n_event = sum(dat$event),
      n_nonevent = sum(!dat$event),
      estimate = `exp(coef)`,
      lower_ci = `lower .95`,
      upper_ci = `upper .95`,
      p_value = clean_p(`Pr(>|z|)`),
      estimate_scale = "hazard ratio for event per 1 SD higher growth",
      tradeoff_direction = "hazard ratio below 1 would support slower regeneration in faster-growing fragments",
      diagnostic_figure = diag_path
    )

  diagnostics <- tibble(
    model = model_label,
    diagnostic = c("global proportional hazards",
                   "growth proportional hazards",
                   "maximum absolute dfbeta"),
    value = c(global_ph_p, growth_ph_p, max_dfbeta),
    p_value = c(global_ph_p, growth_ph_p, NA_real_),
    status = c(plain_status(global_ph_p),
               plain_status(growth_ph_p),
               if_else(max_dfbeta < 1, "passed", "flagged")),
    diagnostic_figure = diag_path
  )

  list(model = model, result = result, diagnostics = diagnostics, data = dat)
}

make_transition_table <- function(phys, bw, trait) {
  phys |>
    filter(wound == "yes") |>
    arrange(id, day) |>
    group_by(id, treatment, thicket, tank) |>
    mutate(
      value = replace_na(as.integer(.data[[trait]]), 0L),
      prior_event = lag(cummax(value), default = 0L),
      prev_day = lag(day),
      interval_days = day - prev_day,
      day_end = day
    ) |>
    filter(!is.na(prev_day),
           interval_days > 0,
           prior_event == 0,
           !is.na(value)) |>
    mutate(event = as.integer(value == 1)) |>
    ungroup() |>
    left_join(select(bw, id, pct_growth), by = "id") |>
    mutate(
      treatment = factor(treatment, levels = c("28C", "31C")),
      thicket = factor(thicket, levels = c("a", "c", "d")),
      z_growth = as.numeric(scale(pct_growth)),
      day_end_c = as.numeric(scale(day_end)),
      log_interval = log(interval_days)
    )
}

fit_discrete_transition <- function(phys, bw, trait, event_label) {
  dat <- make_transition_table(phys, bw, trait)
  model_label <- paste0("cloglog_transition_", trait, "_growth")

  model <- glm(
    event ~ treatment + z_growth + day_end_c + thicket +
      offset(log_interval),
    data = dat,
    family = binomial(link = "cloglog"),
    method = brglmFit,
    type = "AS_mean"
  )

  cluster_test <- tryCatch({
    lmtest::coeftest(model,
                     vcov. = sandwich::vcovCL(model, cluster = dat$id))
  }, error = function(e) NULL)

  if (!is.null(cluster_test) && "z_growth" %in% rownames(cluster_test)) {
    growth_est <- cluster_test["z_growth", "Estimate"]
    growth_se <- cluster_test["z_growth", "Std. Error"]
    growth_p <- clean_p(cluster_test["z_growth", "Pr(>|z|)"])
    se_note <- "cluster-robust standard error by fragment ID"
  } else {
    coefs <- summary(model)$coefficients
    growth_est <- coefs["z_growth", "Estimate"]
    growth_se <- coefs["z_growth", "Std. Error"]
    growth_p <- clean_p(coefs["z_growth", "Pr(>|z|)"])
    se_note <- "model-based standard error; cluster-robust estimate unavailable"
  }
  ci <- safe_ci(growth_est, growth_se)

  diag_path <- file.path(DIAG_DIR, paste0(model_label, "_dharma.png"))
  dharma <- tryCatch({
    set.seed(42)
    sim <- suppressWarnings(
      DHARMa::simulateResiduals(model, n = 1000, plot = FALSE)
    )
    png(diag_path, width = 1800, height = 1500, res = 220)
    suppressWarnings(plot(sim))
    dev.off()
    list(
      uniformity_p = clean_p(DHARMa::testUniformity(sim)$p.value),
      dispersion_p = clean_p(DHARMa::testDispersion(sim)$p.value),
      outlier_p = clean_p(DHARMa::testOutliers(sim)$p.value),
      ok = TRUE
    )
  }, error = function(e) {
    list(uniformity_p = NA_real_, dispersion_p = NA_real_,
         outlier_p = NA_real_, ok = FALSE)
  })

  binned_path <- file.path(DIAG_DIR, paste0(model_label, "_binned.png"))
  pred_df <- dat |>
    mutate(
      predicted = predict(model, type = "response"),
      residual = event - predicted,
      pred_bin = ntile(predicted, pmin(5, n()))
    )
  binned <- pred_df |>
    group_by(pred_bin) |>
    summarise(
      predicted = mean(predicted),
      observed = mean(event),
      n = n(),
      .groups = "drop"
    )
  p_pred <- ggplot(binned, aes(predicted, observed)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.35,
                colour = "grey45") +
    geom_point(aes(size = n), alpha = 0.85) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Mean predicted event chance",
         y = "Observed event frequency",
         title = paste(event_label, "calibration")) +
    theme_pub(9) +
    theme(legend.position = "none")
  p_resid <- ggplot(pred_df, aes(day_end, residual)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
    geom_point(alpha = 0.45, size = 1.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                linewidth = 0.4, colour = "black") +
    labs(x = "Day", y = "Observed - predicted",
         title = paste(event_label, "residuals over time")) +
    theme_pub(9)
  p_events <- pred_df |>
    count(day_end, event) |>
    mutate(event = factor(event, levels = c(0, 1),
                          labels = c("not yet", "transition"))) |>
    ggplot(aes(day_end, n, fill = event)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_manual(values = c("not yet" = "grey80",
                                 "transition" = "#D55E00")) +
    labs(x = "Day", y = "At-risk intervals",
         fill = NULL, title = paste(event_label, "event timing")) +
    theme_pub(9)
  ggsave(binned_path, p_pred + p_resid + p_events + plot_layout(nrow = 1),
         width = 220, height = 72, units = "mm", dpi = 300)

  max_hat <- max(stats::hatvalues(model), na.rm = TRUE)
  model_converged <- isTRUE(model$converged)

  result <- tibble(
    model = model_label,
    analysis_type = "discrete-time transition",
    event = event_label,
    n = nrow(dat),
    n_event = sum(dat$event == 1, na.rm = TRUE),
    n_nonevent = sum(dat$event == 0, na.rm = TRUE),
    estimate = exp(growth_est),
    lower_ci = exp(ci[[1]]),
    upper_ci = exp(ci[[2]]),
    p_value = growth_p,
    estimate_scale = paste(
      "daily transition hazard ratio per 1 SD higher growth;",
      se_note
    ),
    tradeoff_direction = "hazard ratio below 1 would support slower regeneration in faster-growing fragments",
    diagnostic_figure = paste(diag_path, binned_path, sep = " / ")
  )

  diagnostics <- tibble(
    model = model_label,
    diagnostic = c("bias-reduced complementary log-log convergence",
                   "DHARMa uniformity screen",
                   "DHARMa dispersion screen",
                   "DHARMa outlier screen",
                   "maximum hat value"),
    value = c(as.numeric(model_converged),
              dharma$uniformity_p, dharma$dispersion_p, dharma$outlier_p,
              max_hat),
    p_value = c(NA_real_,
                dharma$uniformity_p, dharma$dispersion_p, dharma$outlier_p,
                NA_real_),
    status = c(if_else(model_converged, "passed", "flagged"),
               plain_status(dharma$uniformity_p),
               plain_status(dharma$dispersion_p),
               plain_status(dharma$outlier_p),
               if_else(max_hat < 0.5, "passed", "flagged")),
    diagnostic_figure = paste(diag_path, binned_path, sep = " / ")
  )

  list(model = model, result = result, diagnostics = diagnostics, data = dat)
}

source_slope <- function(x, y) {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  unname(coef(lm(y ~ x))[2])
}

# ---- Valid same-fragment table --------------------------------------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  select(id,
         treatment_growth = treatment,
         wound_growth = wound,
         thicket_growth = thicket,
         tank_growth = tank,
         pct_growth)

phys <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))

phys15 <- phys |>
  filter(day == 15) |>
  select(id,
         treatment_morph = treatment,
         wound_morph = wound,
         thicket_morph = thicket,
         tank_morph = tank,
         axial_polyp_formation,
         wound_smoothed,
         tip_extension,
         new_corallites_on_tip)

growth_morph <- inner_join(bw, phys15, by = "id") |>
  mutate(
    treatment_match = as.character(treatment_growth) ==
      as.character(treatment_morph),
    wound_match = as.character(wound_growth) == as.character(wound_morph),
    thicket_match = thicket_growth == thicket_morph,
    tank_match = tank_growth == tank_morph
  )

join_audit <- tibble(
  check = c(
    "growth rows",
    "Day-15 morphology rows",
    "joined rows",
    "unique joined IDs",
    "treatment-label mismatches",
    "wound-label mismatches",
    "source-patch-label mismatches",
    "tank-label mismatches"
  ),
  value = c(
    nrow(bw),
    nrow(phys15),
    nrow(growth_morph),
    n_distinct(growth_morph$id),
    sum(!growth_morph$treatment_match),
    sum(!growth_morph$wound_match),
    sum(!growth_morph$thicket_match),
    sum(!growth_morph$tank_match)
  ),
  status = c(
    "ok",
    "ok",
    if_else(nrow(growth_morph) == nrow(bw), "ok", "review"),
    if_else(n_distinct(growth_morph$id) == nrow(growth_morph), "ok", "fail"),
    if_else(sum(!growth_morph$treatment_match) == 0, "ok", "fail"),
    if_else(sum(!growth_morph$wound_match) <= 1, "known mismatch", "review"),
    if_else(sum(!growth_morph$thicket_match) == 0, "ok", "fail"),
    if_else(sum(!growth_morph$tank_match) == 0, "ok", "fail")
  ),
  notes = c(
    "One growth row per physiology fragment.",
    "One Day-15 morphology row per physiology fragment.",
    "Joined by coral fragment ID.",
    "Duplicate IDs would invalidate same-fragment trade-off tests.",
    "Must be zero for same-fragment treatment tests.",
    "Morphology wound labels define the wounded-only regeneration screen.",
    "Must be zero for source-patch adjustment.",
    "Must be zero for treatment-tank consistency."
  )
)

write_csv(join_audit, file.path(TBL_DIR, "37_tradeoff_join_audit.csv"))

if (any(!growth_morph$treatment_match) ||
    any(!growth_morph$thicket_match) ||
    any(!growth_morph$tank_match) ||
    n_distinct(growth_morph$id) != nrow(growth_morph)) {
  stop("Trade-off join audit failed on treatment, source patch, tank, or IDs.")
}

wounded <- growth_morph |>
  filter(wound_morph == "yes") |>
  mutate(
    treatment = factor(treatment_morph, levels = c("28C", "31C")),
    thicket = factor(thicket_morph, levels = c("a", "c", "d")),
    source_label = str_to_upper(as.character(thicket))
  )

wounded_heat <- wounded |>
  filter(treatment == "31C")

snp_joinability <- read_csv(
  file.path(TBL_DIR, "32_prelim_snp_response_joinability.csv"),
  show_col_types = FALSE
)

scope_inventory <- tibble(
  candidate = c(
    "growth vs Day-15 new corallites",
    "growth vs Day-15 tip extension",
    "growth vs timing of wound smoothing",
    "growth vs timing of tip extension",
    "growth vs timing of new corallites",
    "growth vs daily regeneration transitions",
    "ambient growth vs growth heat drop",
    "ambient growth vs non-growth heat penalty",
    "condition heat penalty vs heated regeneration",
    "SNP cluster vs growth or regeneration"
  ),
  evidence_level = c(
    rep("same fragment", 6),
    rep("source-patch summary", 3),
    "not directly joinable"
  ),
  valid_analysis = c(
    "linear model, robust regression, and heated-only exact permutation sensitivity",
    "linear model, robust regression, and heated-only exact permutation sensitivity",
    "Cox time-to-event model",
    "Cox time-to-event model",
    "Cox time-to-event model",
    "bias-reduced discrete-time transition model with complementary log-log link",
    "descriptive bootstrap only; n = 3 source patches",
    "leave-growth-out heat-penalty index; descriptive bootstrap only; n = 3 source patches",
    "descriptive bootstrap only; n = 3 source patches",
    "do not test yet"
  ),
  diagnostic_plan = c(
    "residual plots, Q-Q plot, Shapiro-Wilk, Cook's distance, leave-one-out, robust-regression weights",
    "residual plots, Q-Q plot, Shapiro-Wilk, Cook's distance, leave-one-out, robust-regression weights",
    "Schoenfeld proportional hazards check, dfbeta influence",
    "Schoenfeld proportional hazards check, dfbeta influence",
    "Schoenfeld proportional hazards check, dfbeta influence",
    "convergence, DHARMa screens, calibration plot, residual-over-time plot, leverage",
    "no model diagnostics; uncertainty summarized by bootstrap sign",
    "no model diagnostics; uncertainty summarized by bootstrap sign",
    "no model diagnostics; uncertainty summarized by bootstrap sign",
    "joinability audit only"
  ),
  current_status = c(
    rep("tested here", 9),
    paste0(
      "blocked for growth/regeneration; direct SNP join currently only for ",
      paste(snp_joinability$response_variable[
        snp_joinability$direct_sample_level_join
      ], collapse = ", ")
    )
  )
)

write_csv(scope_inventory,
          file.path(TBL_DIR, "37_tradeoff_scope_inventory.csv"))

# ---- Endpoint growth/regeneration models ----------------------------------
endpoint_fits <- list(
  fit_lm_growth_cost(wounded, "new_corallites_on_tip",
                     "new corallites by Day 15",
                     "lm_growth_by_new_corallites_all_wounded",
                     include_treatment = TRUE),
  fit_lm_growth_cost(wounded_heat, "new_corallites_on_tip",
                     "new corallites by Day 15",
                     "lm_growth_by_new_corallites_heated_wounded",
                     include_treatment = FALSE),
  fit_lm_growth_cost(wounded, "tip_extension",
                     "tip extension by Day 15",
                     "lm_growth_by_tip_extension_all_wounded",
                     include_treatment = TRUE),
  fit_lm_growth_cost(wounded_heat, "tip_extension",
                     "tip extension by Day 15",
                     "lm_growth_by_tip_extension_heated_wounded",
                     include_treatment = FALSE),
  fit_brglm_regeneration(wounded, "new_corallites_on_tip",
                         "new corallites by Day 15",
                         "brglm_new_corallites_by_growth_all_wounded",
                         include_treatment = TRUE),
  fit_brglm_regeneration(wounded_heat, "new_corallites_on_tip",
                         "new corallites by Day 15",
                         "brglm_new_corallites_by_growth_heated_wounded",
                         include_treatment = FALSE),
  fit_brglm_regeneration(wounded, "tip_extension",
                         "tip extension by Day 15",
                         "brglm_tip_extension_by_growth_all_wounded",
                         include_treatment = TRUE),
  fit_brglm_regeneration(wounded_heat, "tip_extension",
                         "tip extension by Day 15",
                         "brglm_tip_extension_by_growth_heated_wounded",
                         include_treatment = FALSE)
)

endpoint_results <- map_dfr(endpoint_fits, "result")
endpoint_diagnostics <- map_dfr(endpoint_fits, "diagnostics")
endpoint_influence <- bind_rows(map(endpoint_fits, \(x) x$influence))
endpoint_leave_one_out <- bind_rows(map(endpoint_fits, \(x) x$leave_one_out))

robust_fits <- list(
  fit_robust_growth_cost(wounded, "new_corallites_on_tip",
                         "new corallites by Day 15",
                         "lmrob_growth_by_new_corallites_all_wounded",
                         include_treatment = TRUE),
  fit_robust_growth_cost(wounded_heat, "new_corallites_on_tip",
                         "new corallites by Day 15",
                         "lmrob_growth_by_new_corallites_heated_wounded",
                         include_treatment = FALSE),
  fit_robust_growth_cost(wounded, "tip_extension",
                         "tip extension by Day 15",
                         "lmrob_growth_by_tip_extension_all_wounded",
                         include_treatment = TRUE),
  fit_robust_growth_cost(wounded_heat, "tip_extension",
                         "tip extension by Day 15",
                         "lmrob_growth_by_tip_extension_heated_wounded",
                         include_treatment = FALSE)
)

robust_results <- map_dfr(robust_fits, "result")
robust_diagnostics <- map_dfr(robust_fits, "diagnostics")

heated_permutation <- bind_rows(
  exact_perm_mean_diff(wounded_heat$pct_growth,
                       wounded_heat$new_corallites_on_tip) |>
    mutate(event = "new corallites by Day 15"),
  exact_perm_mean_diff(wounded_heat$pct_growth,
                       wounded_heat$tip_extension) |>
    mutate(event = "tip extension by Day 15")
) |>
  mutate(
    model = paste0("exact_permutation_growth_difference_heated_",
                   str_replace_all(str_to_lower(event), "[^a-z0-9]+", "_")),
    analysis_type = "heated-only exact permutation",
    estimate = raw_yes_minus_no,
    lower_ci = NA_real_,
    upper_ci = NA_real_,
    p_value = permutation_p,
    estimate_scale = "raw percentage-point growth difference: event yes minus no",
    tradeoff_direction = "negative estimate would support a growth cost of regeneration",
    diagnostic_figure = NA_character_
  ) |>
  select(model, analysis_type, event, n_event, n_nonevent,
         estimate, lower_ci, upper_ci, p_value, estimate_scale,
         tradeoff_direction, diagnostic_figure)

endpoint_results <- bind_rows(endpoint_results, heated_permutation)

write_csv(endpoint_results,
          file.path(TBL_DIR, "37_tradeoff_endpoint_results.csv"))
write_csv(robust_results,
          file.path(TBL_DIR, "37_tradeoff_robust_growth_results.csv"))
write_csv(endpoint_influence,
          file.path(TBL_DIR, "37_tradeoff_endpoint_influence.csv"))
write_csv(endpoint_leave_one_out,
          file.path(TBL_DIR, "37_tradeoff_endpoint_leave_one_out.csv"))

# ---- Time-to-event models --------------------------------------------------
cox_fits <- list(
  fit_cox_timing(phys, bw, "wound_smoothed", "wound smoothing"),
  fit_cox_timing(phys, bw, "tip_extension", "tip extension"),
  fit_cox_timing(phys, bw, "new_corallites_on_tip", "new corallites")
)

cox_results <- map_dfr(cox_fits, "result")
cox_diagnostics <- map_dfr(cox_fits, "diagnostics")

write_csv(cox_results,
          file.path(TBL_DIR, "37_tradeoff_time_to_event_results.csv"))

transition_fits <- list(
  fit_discrete_transition(phys, bw, "wound_smoothed", "wound smoothing"),
  fit_discrete_transition(phys, bw, "tip_extension", "tip extension"),
  fit_discrete_transition(phys, bw, "new_corallites_on_tip", "new corallites")
)

transition_results <- map_dfr(transition_fits, "result")
transition_diagnostics <- map_dfr(transition_fits, "diagnostics")
transition_at_risk <- map_dfr(transition_fits, \(x) {
  x$data |>
    count(model = x$result$model[[1]], event_label = x$result$event[[1]],
          treatment, day_end, event, name = "n_intervals")
})

write_csv(transition_results,
          file.path(TBL_DIR, "37_tradeoff_discrete_transition_results.csv"))
write_csv(transition_at_risk,
          file.path(TBL_DIR, "37_tradeoff_discrete_transition_at_risk.csv"))

# ---- Source-level screens --------------------------------------------------
growth_means <- read_csv(file.path(TBL_DIR, "13_genet_emmeans.csv"),
                         show_col_types = FALSE) |>
  filter(response == "Growth (% mass change)") |>
  mutate(source_label = str_to_upper(thicket))

growth_tradeoff <- growth_means |>
  select(thicket, source_label, treatment, mean, se, n) |>
  pivot_wider(names_from = treatment,
              values_from = c(mean, se, n),
              names_sep = "_") |>
  mutate(
    ambient_growth_pct = mean_28C,
    heated_growth_pct = mean_31C,
    growth_heat_drop_pct = mean_28C - mean_31C,
    growth_heat_drop_se = sqrt(se_28C^2 + se_31C^2)
  )

heated_regen_by_source <- wounded_heat |>
  group_by(thicket, source_label) |>
  summarise(
    n_heated_wounded = n(),
    new_corallites_count = sum(new_corallites_on_tip == 1, na.rm = TRUE),
    tip_extension_count = sum(tip_extension == 1, na.rm = TRUE),
    new_corallites_pct = 100 * new_corallites_count / n_heated_wounded,
    tip_extension_pct = 100 * tip_extension_count / n_heated_wounded,
    heated_growth_wounded_pct = mean(pct_growth, na.rm = TRUE),
    .groups = "drop"
  )

heat_effects <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                         show_col_types = FALSE) |>
  filter(is.finite(estimate), is.finite(SE))

heat_index <- function(response_set, value_col = "estimate") {
  heat_effects |>
    filter(response %in% response_set) |>
    group_by(response, thicket) |>
    summarise(effect = mean(.data[[value_col]], na.rm = TRUE),
              .groups = "drop") |>
    group_by(response) |>
    mutate(scaled = effect / max(abs(effect), na.rm = TRUE)) |>
    ungroup() |>
    group_by(thicket) |>
    summarise(index = mean(scaled, na.rm = TRUE),
              n_responses = n(),
              .groups = "drop")
}

whole_fragment_responses <- c("pam_fvfm", "color_dscale", "growth_pct",
                              "log_zoox_density")
non_growth_responses <- setdiff(whole_fragment_responses, "growth_pct")
non_symbiont_responses <- setdiff(whole_fragment_responses,
                                  "log_zoox_density")

heat_penalty_sets <- tibble(
  index_name = c("all whole-fragment traits",
                 "leave growth out",
                 "leave photosynthesis score out",
                 "leave color out",
                 "leave symbiont density out",
                 "leave regeneration out"),
  omitted_response = c("none",
                       "skeletal growth",
                       "photosynthesis score",
                       "color",
                       "symbiont density",
                       "all regeneration traits"),
  response_set = list(
    whole_fragment_responses,
    non_growth_responses,
    setdiff(whole_fragment_responses, "pam_fvfm"),
    setdiff(whole_fragment_responses, "color_dscale"),
    non_symbiont_responses,
    whole_fragment_responses
  )
)

heat_penalty_leave_one_out <- heat_penalty_sets |>
  mutate(values = map2(index_name, response_set, \(nm, responses) {
    heat_index(responses) |>
      mutate(responses_used = paste(responses, collapse = "; "))
  })) |>
  select(index_name, omitted_response, values) |>
  unnest(values) |>
  rename(heat_penalty = index) |>
  select(index_name, omitted_response, thicket, heat_penalty,
         n_responses, responses_used)

write_csv(heat_penalty_leave_one_out,
          file.path(TBL_DIR,
                    "37_tradeoff_heat_penalty_leave_one_trait_out.csv"))

non_growth_heat_index <- heat_index(non_growth_responses) |>
  rename(non_growth_heat_penalty = index,
         non_growth_heat_n = n_responses)

non_regeneration_heat_index <- heat_index(whole_fragment_responses) |>
  rename(non_regeneration_heat_penalty = index,
         non_regeneration_heat_n = n_responses)

source_screen <- growth_tradeoff |>
  left_join(heated_regen_by_source, by = c("thicket", "source_label")) |>
  left_join(non_growth_heat_index, by = "thicket") |>
  left_join(non_regeneration_heat_index, by = "thicket")

B <- 5000

source_screen <- source_screen |>
  arrange(thicket)

draw_normal_matrix <- function(mean, se, B) {
  matrix(
    rnorm(B * length(mean),
          mean = rep(mean, each = B),
          sd = rep(se, each = B)),
    nrow = B,
    ncol = length(mean)
  )
}

draw_beta_matrix <- function(success, total, B) {
  matrix(
    rbeta(B * length(success),
          shape1 = rep(success + 0.5, each = B),
          shape2 = rep(total - success + 0.5, each = B)) * 100,
    nrow = B,
    ncol = length(success)
  )
}

matrix_slope <- function(x, y) {
  x_centered <- sweep(x, 1, rowMeans(x), "-")
  y_centered <- sweep(y, 1, rowMeans(y), "-")
  rowSums(x_centered * y_centered) / rowSums(x_centered^2)
}

draw_heat_penalty_matrix <- function(response_set) {
  base <- heat_effects |>
    filter(response %in% response_set) |>
    group_by(response, thicket) |>
    summarise(
      estimate = mean(estimate, na.rm = TRUE),
      se = sqrt(sum(SE^2, na.rm = TRUE)) / n(),
      .groups = "drop"
    ) |>
    arrange(response, thicket)

  draws <- draw_normal_matrix(base$estimate, base$se, B)
  n_response <- n_distinct(base$response)
  n_source <- n_distinct(base$thicket)
  draw_array <- array(draws, dim = c(B, n_source, n_response))
  scaled <- draw_array
  for (j in seq_len(n_response)) {
    denom <- apply(abs(draw_array[, , j, drop = FALSE]), 1, max)
    scaled[, , j] <- draw_array[, , j] / denom
  }
  apply(scaled, c(1, 2), mean)
}

non_growth_penalty_mat <- draw_heat_penalty_matrix(non_growth_responses)
non_regeneration_penalty_mat <- draw_heat_penalty_matrix(
  whole_fragment_responses
)

ambient_growth_mat <- draw_normal_matrix(source_screen$ambient_growth_pct,
                                         source_screen$se_28C, B)
growth_drop_mat <- draw_normal_matrix(source_screen$growth_heat_drop_pct,
                                      source_screen$growth_heat_drop_se, B)
new_corallites_mat <- draw_beta_matrix(source_screen$new_corallites_count,
                                       source_screen$n_heated_wounded, B)
tip_extension_mat <- draw_beta_matrix(source_screen$tip_extension_count,
                                      source_screen$n_heated_wounded, B)

summarise_source_screen <- function(screen, x_variable, y_variable,
                                    observed_x, observed_y,
                                    x_mat, y_mat,
                                    expected_tradeoff_direction) {
  slopes <- matrix_slope(x_mat, y_mat)
  prob_direction <- if (expected_tradeoff_direction == "positive") {
    mean(slopes > 0, na.rm = TRUE)
  } else {
    mean(slopes < 0, na.rm = TRUE)
  }

  tibble(
    screen = screen,
    x_variable = x_variable,
    y_variable = y_variable,
    n_source_patches = nrow(source_screen),
    observed_slope = source_slope(observed_x, observed_y),
    boot_slope_median = median(slopes, na.rm = TRUE),
    boot_slope_lo = quantile(slopes, 0.025, na.rm = TRUE),
    boot_slope_hi = quantile(slopes, 0.975, na.rm = TRUE),
    expected_tradeoff_direction = expected_tradeoff_direction,
    boot_prob_tradeoff_direction = prob_direction,
    interpretation = case_when(
      screen %in% c("ambient growth vs growth heat drop",
                    "ambient growth vs non-growth heat penalty") &
        prob_direction >= 0.8 ~
        "possible source-level growth/heat-tolerance trade-off",
      screen %in% c("ambient growth vs growth heat drop",
                    "ambient growth vs non-growth heat penalty") ~
        "weak source-level growth/heat-tolerance trade-off screen",
      prob_direction <= 0.2 ~
        "pattern is opposite the expected trade-off direction",
      TRUE ~
        "little support for expected trade-off direction"
    )
  )
}

source_results <- bind_rows(
  summarise_source_screen(
    "ambient growth vs growth heat drop",
    "growth at 28C", "growth lost under heat",
    source_screen$ambient_growth_pct, source_screen$growth_heat_drop_pct,
    ambient_growth_mat, growth_drop_mat, "positive"
  ),
  summarise_source_screen(
    "ambient growth vs non-growth heat penalty",
    "growth at 28C", "heat penalty excluding growth",
    source_screen$ambient_growth_pct, source_screen$non_growth_heat_penalty,
    ambient_growth_mat, non_growth_penalty_mat, "positive"
  ),
  summarise_source_screen(
    "ambient growth vs heated new corallites",
    "growth at 28C", "new corallites under heat",
    source_screen$ambient_growth_pct, source_screen$new_corallites_pct,
    ambient_growth_mat, new_corallites_mat, "negative"
  ),
  summarise_source_screen(
    "ambient growth vs heated tip extension",
    "growth at 28C", "tip extension under heat",
    source_screen$ambient_growth_pct, source_screen$tip_extension_pct,
    ambient_growth_mat, tip_extension_mat, "negative"
  ),
  summarise_source_screen(
    "growth heat drop vs heated new corallites",
    "growth lost under heat", "new corallites under heat",
    source_screen$growth_heat_drop_pct, source_screen$new_corallites_pct,
    growth_drop_mat, new_corallites_mat, "positive"
  ),
  summarise_source_screen(
    "non-regeneration heat penalty vs heated new corallites",
    "heat penalty excluding regeneration", "new corallites under heat",
    source_screen$non_regeneration_heat_penalty,
    source_screen$new_corallites_pct,
    non_regeneration_penalty_mat, new_corallites_mat, "positive"
  ),
  summarise_source_screen(
    "non-regeneration heat penalty vs heated tip extension",
    "heat penalty excluding regeneration", "tip extension under heat",
    source_screen$non_regeneration_heat_penalty,
    source_screen$tip_extension_pct,
    non_regeneration_penalty_mat, tip_extension_mat, "positive"
  )
)

write_csv(source_screen,
          file.path(TBL_DIR, "37_tradeoff_source_patch_values.csv"))
write_csv(source_results,
          file.path(TBL_DIR, "37_tradeoff_source_level_screen.csv"))

# ---- Diagnostics and result synthesis -------------------------------------
all_diagnostics <- bind_rows(endpoint_diagnostics, robust_diagnostics,
                             cox_diagnostics, transition_diagnostics) |>
  arrange(model, diagnostic)

write_csv(all_diagnostics,
          file.path(TBL_DIR, "37_tradeoff_model_diagnostics.csv"))

main_results <- bind_rows(
  endpoint_results |>
    filter(
      model %in% c(
        "lm_growth_by_new_corallites_heated_wounded",
        "lm_growth_by_tip_extension_heated_wounded",
        "brglm_new_corallites_by_growth_heated_wounded",
        "brglm_tip_extension_by_growth_heated_wounded"
      ) |
        analysis_type == "heated-only exact permutation"
    ) |>
    mutate(
      result_family = "same-fragment endpoint",
      plain_language = case_when(
        str_detect(model, "^lm_") & p_value >= 0.05 & estimate < 0 ~
          "adjusted estimate is slightly negative but too uncertain to support a growth cost",
        str_detect(model, "^lm_") & p_value >= 0.05 & estimate >= 0 ~
          "fragments with regeneration did not have lower growth",
        str_detect(model, "^lm_") & estimate < 0 ~
          "fragments with regeneration had lower growth",
        str_detect(model, "^lm_") ~
          "fragments with regeneration had higher growth",
        str_detect(model, "^brglm") & p_value >= 0.05 ~
          "no clear relationship between growth and regeneration odds",
        str_detect(model, "^brglm") & estimate < 1 ~
          "higher-growth fragments were less likely to regenerate",
        str_detect(model, "^brglm") ~
          "higher-growth fragments were at least as likely to regenerate",
        analysis_type == "heated-only exact permutation" & estimate < 0 ~
          "raw heated-only difference points toward a growth cost",
        analysis_type == "heated-only exact permutation" & estimate >= 0 ~
          "raw heated-only difference does not point toward a growth cost",
        TRUE ~ "see estimate"
      )
    ) |>
    select(result_family, model, event, n_event, n_nonevent, estimate,
           lower_ci, upper_ci, p_value, estimate_scale, tradeoff_direction,
           plain_language),
  robust_results |>
    filter(
      model %in% c(
        "lmrob_growth_by_new_corallites_heated_wounded",
        "lmrob_growth_by_tip_extension_heated_wounded"
      )
    ) |>
    mutate(
      result_family = "same-fragment robust endpoint",
      plain_language = case_when(
        p_value >= 0.05 & estimate < 0 ~
          "robust estimate is slightly negative but too uncertain to support a growth cost",
        p_value >= 0.05 & estimate >= 0 ~
          "robust model does not show lower growth in regenerating fragments",
        estimate < 0 ~
          "robust model points toward a growth cost",
        TRUE ~
          "robust model points toward higher growth in regenerating fragments"
      )
    ) |>
    select(result_family, model, event, n_event, n_nonevent, estimate,
           lower_ci, upper_ci, p_value, estimate_scale, tradeoff_direction,
           plain_language),
  cox_results |>
    mutate(
      result_family = "same-fragment timing",
      plain_language = if_else(
        estimate < 1,
        "higher-growth fragments tended to reach the stage more slowly",
        "higher-growth fragments did not reach the stage more slowly"
      )
    ) |>
    select(result_family, model, event, n_event, n_nonevent, estimate,
           lower_ci, upper_ci, p_value, estimate_scale, tradeoff_direction,
           plain_language),
  transition_results |>
    mutate(
      result_family = "same-fragment daily transition",
      plain_language = case_when(
        p_value >= 0.05 & estimate < 1 ~
          "daily transition model is slightly below 1 but too uncertain to support slower regeneration in faster-growing fragments",
        p_value >= 0.05 & estimate >= 1 ~
          "daily transition model does not show slower regeneration in faster-growing fragments",
        estimate < 1 ~
          "daily transition model supports slower regeneration in faster-growing fragments",
        TRUE ~
          "daily transition model supports faster regeneration in faster-growing fragments"
      )
    ) |>
    select(result_family, model, event, n_event, n_nonevent, estimate,
           lower_ci, upper_ci, p_value, estimate_scale, tradeoff_direction,
           plain_language),
  source_results |>
    transmute(
      result_family = "source-patch screen",
      model = screen,
      event = y_variable,
      n_event = NA_integer_,
      n_nonevent = NA_integer_,
      estimate = observed_slope,
      lower_ci = boot_slope_lo,
      upper_ci = boot_slope_hi,
      p_value = NA_real_,
      estimate_scale = "source-level slope; descriptive only",
      tradeoff_direction = paste("expected slope direction:",
                                 expected_tradeoff_direction),
      plain_language = interpretation
    )
)

write_csv(main_results,
          file.path(TBL_DIR, "37_tradeoff_main_results.csv"))

cat("\n=== Formal trade-off main results ===\n")
print(main_results |>
        mutate(across(where(is.numeric), \(x) round(x, 3))) |>
        select(result_family, event, estimate, lower_ci, upper_ci,
               p_value, plain_language),
      n = Inf)

cat("\n=== Diagnostic status counts ===\n")
print(all_diagnostics |> count(status))
