# =============================================================================
# Purpose: Add Trinity Conn / Ross Cunning Hauru colony metadata as external
#          context for future genetic matching.
#
#          Molly forwarded Trinity Conn's 2026-08-10 email on 2026-09-03.
#          Trinity attached Hauru colony coordinates, acute CBASS ED50 values,
#          clonal groups, and dominant symbiont calls. Molly noted that our
#          source patch C coordinate may match their Apul-115, but Trinity was
#          clear that the exact test requires calling SNPs from our RNA-seq data
#          and comparing them to their whole-genome data.
#
# What & why: this script keeps that preliminary external context reproducible without
#   treating it as LTH phenotype data. It standardizes Trinity's spreadsheet,
#   converts our source-patch coordinates from south/west magnitudes to signed
#   decimal degrees, finds the nearest Hauru colonies to source patches A/C/D,
#   and makes a simple map. Geographic proximity is planning evidence only; it
#   is not a genotype match.
#
# Input:   data/external/trinity_conn_hauru_colonies_2026-08-10.xlsx
#          data/raw/metadata/metadata.csv
# Output:  data/external/trinity_conn_hauru_colonies_2026-08-10.csv
#          output/tables/26b_trinity_hauru_table_summary.csv
#          output/tables/26b_trinity_hauru_nearest_source_patch_matches.csv
#          figures/26b_trinity_hauru_context.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

external_dir <- here::here("data", "external")
hauru_xlsx <- file.path(external_dir, "trinity_conn_hauru_colonies_2026-08-10.xlsx")
hauru_csv <- file.path(external_dir, "trinity_conn_hauru_colonies_2026-08-10.csv")
parse_number_or_na <- function(x) {
  suppressWarnings(parse_double(na_if(as.character(x), "NA")))
}

if (!file.exists(hauru_xlsx)) {
  cat("\n=== Trinity Hauru context skipped ===\n")
  cat("Missing data/external/trinity_conn_hauru_colonies_2026-08-10.xlsx\n")
} else {
  hauru <- read_excel(hauru_xlsx) |>
    clean_names() |>
    transmute(
      site = str_to_lower(site),
      colony_id,
      lat = parse_number_or_na(lat),
      lon = parse_number_or_na(lon),
      ed50 = as.numeric(ed50),
      clonal_group = as.character(clonal_group),
      dominant_sym = as.character(dominant_sym)
    )

  write_csv(hauru, hauru_csv)

  source_patch_coords <- suppressWarnings(
    read_csv(file.path(DATA_RAW, "metadata", "metadata.csv"),
             col_types = cols(.default = col_character()))
  ) |>
    clean_names() |>
    filter(thicket %in% c("a", "c", "d")) |>
    mutate(
      coord_lat = parse_number_or_na(coord_lat),
      coord_long = parse_number_or_na(coord_long)
    ) |>
    group_by(thicket) |>
    summarise(
      n_fragments = n(),
      lth_lat_raw = first(coord_lat),
      lth_lon_raw = first(coord_long),
      .groups = "drop"
    ) |>
    mutate(
      source_patch = str_to_upper(thicket),
      # LTH metadata store Moorea south/west coordinates as positive magnitudes.
      lth_lat = if_else(lth_lat_raw > 0, -lth_lat_raw, lth_lat_raw),
      lth_lon = if_else(lth_lon_raw > 0, -lth_lon_raw, lth_lon_raw)
    ) |>
    select(source_patch, thicket, n_fragments, lth_lat_raw, lth_lon_raw,
           lth_lat, lth_lon)

  haversine_m <- function(lat1, lon1, lat2, lon2) {
    radius_m <- 6371000
    to_rad <- pi / 180
    phi1 <- lat1 * to_rad
    phi2 <- lat2 * to_rad
    dphi <- (lat2 - lat1) * to_rad
    dlambda <- (lon2 - lon1) * to_rad
    a <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
    2 * radius_m * asin(pmin(1, sqrt(a)))
  }

  hauru_with_coords <- hauru |>
    filter(!is.na(lat), !is.na(lon))

  nearest <- tidyr::crossing(source_patch_coords, hauru_with_coords) |>
    mutate(distance_m = haversine_m(lth_lat, lth_lon, lat, lon)) |>
    group_by(source_patch) |>
    arrange(distance_m, .by_group = TRUE) |>
    mutate(match_rank = row_number()) |>
    filter(match_rank <= 3) |>
    ungroup() |>
    transmute(
      source_patch,
      lth_lat,
      lth_lon,
      nearest_rank = match_rank,
      hauru_colony_id = colony_id,
      hauru_lat = lat,
      hauru_lon = lon,
      distance_m = round(distance_m, 1),
      hauru_ed50 = round(ed50, 3),
      hauru_clonal_group = clonal_group,
      hauru_dominant_sym = dominant_sym,
      interpretation = if_else(
        nearest_rank == 1,
        "closest coordinate match; requires SNP comparison for exact identity",
        "nearby external colony; context only"
      )
    )

  write_csv(nearest, file.path(TBL_DIR, "26b_trinity_hauru_nearest_source_patch_matches.csv"))

  summary_tbl <- hauru |>
    summarise(
      source = "Trinity Conn 2026-08-10 Hauru colony table",
      n_colonies = n(),
      n_with_coordinates = sum(!is.na(lat) & !is.na(lon)),
      n_clonal_groups = n_distinct(clonal_group),
      ed50_min = round(min(ed50, na.rm = TRUE), 3),
      ed50_mean = round(mean(ed50, na.rm = TRUE), 3),
      ed50_max = round(max(ed50, na.rm = TRUE), 3)
    )
  write_csv(summary_tbl, file.path(TBL_DIR, "26b_trinity_hauru_table_summary.csv"))

  nearest_labels <- nearest |>
    filter(nearest_rank == 1) |>
    distinct(hauru_colony_id, hauru_lon, hauru_lat) |>
    mutate(
      label_x = case_when(
        hauru_colony_id == "Apul-115" ~ hauru_lon - 0.00010,
        TRUE ~ hauru_lon + 0.00010
      ),
      label_y = case_when(
        hauru_colony_id == "Apul-115" ~ hauru_lat + 0.00010,
        TRUE ~ hauru_lat + 0.00012
      )
    )

  source_patch_labels <- source_patch_coords |>
    mutate(
      label_x = case_when(
        source_patch == "A" ~ lth_lon - 0.00005,
        source_patch == "D" ~ lth_lon - 0.00015,
        TRUE ~ lth_lon - 0.00012
      ),
      label_y = case_when(
        source_patch == "A" ~ lth_lat - 0.00017,
        source_patch == "D" ~ lth_lat - 0.00012,
        TRUE ~ lth_lat - 0.00013
      )
    )

  x_vals <- c(hauru_with_coords$lon, source_patch_coords$lth_lon)
  y_vals <- c(hauru_with_coords$lat, source_patch_coords$lth_lat)
  x_pad <- max(diff(range(x_vals, na.rm = TRUE)) * 0.10, 0.00022)
  y_pad <- max(diff(range(y_vals, na.rm = TRUE)) * 0.10, 0.00022)

  p <- ggplot() +
    geom_point(data = hauru_with_coords,
               aes(lon, lat, fill = ed50),
               shape = 21, size = 2.8, alpha = 0.85, colour = "grey20") +
    geom_text(data = nearest_labels,
              aes(label_x, label_y, label = hauru_colony_id),
              size = 2.4, colour = "grey15") +
    geom_point(data = source_patch_coords,
               aes(lth_lon, lth_lat),
               shape = 8, size = 4.2, stroke = 1.1, colour = "#1F3864") +
    geom_text(data = source_patch_labels,
              aes(label_x, label_y, label = source_patch),
              size = 2.6, fontface = "bold", colour = "#1F3864") +
    scale_fill_viridis_c(name = "Acute ED50 (°C)", option = "C") +
    scale_x_continuous(labels = label_number(accuracy = 0.0005), breaks = pretty_breaks(3)) +
    scale_y_continuous(labels = label_number(accuracy = 0.0005), breaks = pretty_breaks(4)) +
    coord_quickmap(xlim = range(x_vals, na.rm = TRUE) + c(-x_pad, x_pad),
                   ylim = range(y_vals, na.rm = TRUE) + c(-y_pad, y_pad),
                   expand = FALSE) +
    labs(
      title = "Hauru colonies near LTH patches",
      subtitle = paste0(
        "Preliminary coordinate matches only.\n",
        "Genetic identity still requires DNA-marker comparison."
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    guides(fill = guide_colorbar(barwidth = unit(35, "mm"), barheight = unit(3.5, "mm"))) +
    theme_pub(8) +
    theme(
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey35"),
      plot.margin = margin(6, 6, 6, 6, "pt")
    )

  # Match the canvas to the north-south extent instead of padding a landscape
  # export with blank space. The geographic aspect ratio remains unchanged.
  save_fig(p, "26b_trinity_hauru_context", width = 140, height = 190)

  cat("\n=== Trinity Hauru context ===\n")
  print(as.data.frame(summary_tbl))
  print(as.data.frame(nearest))
  cat("Wrote 26b_trinity_hauru_* tables, cleaned external CSV, and 26b_trinity_hauru_context.{pdf,png}\n")
}
