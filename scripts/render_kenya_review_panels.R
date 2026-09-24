#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
  library(terra)
})

source("R/functions/project_logging.R")
source("R/functions/ndvi_dates.R")

config <- read_project_config()
run <- new_run_context("kenya-review-panels", config$country$iso3)
log_event(run, event = "run_started", message = "Rendering Kenya calibration review panels")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
report_dir <- file.path(data_root, "models", "validation_reports", "kenya_detectability")
panel_dir <- file.path(report_dir, "review_panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

sample_file <- file.path(report_dir, "review_sample.csv")
detectability_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")
monthly_rain_file <- file.path(base_dir, "KEN_baseline_monthly_rainfall.parquet")
phenology_file <- file.path(base_dir, "countries", "KEN_seasonal-phenology.parquet")
windows_file <- file.path(base_dir, "KEN_candidate_season_windows.parquet")
ndvi_dir <- file.path(data_root, "climate_raw", "glass_ndvi_tif")

required <- c(sample_file, detectability_file, monthly_rain_file, phenology_file, windows_file)
if (!all(file.exists(required))) stop("Missing review-panel inputs.")

review <- fread(sample_file)
review_pixels <- sort(unique(review$pixel))
review_key <- unique(review[, .(pixel, admin1_name, review_stratum)])
review_key[, facet_label := paste(admin1_name, pixel, sep = " | ")]

ndvi_files <- sort(list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE))
ndvi_dates <- parse_glass_ndvi_dates(ndvi_files)
keep <- as.integer(format(ndvi_dates, "%Y")) >= config$baseline$start_year &
  as.integer(format(ndvi_dates, "%Y")) <= config$baseline$end_year
ndvi_files <- ndvi_files[keep]
ndvi_dates <- ndvi_dates[keep]

ndvi_stack <- rast(ndvi_files)
ndvi_matrix <- as.matrix(terra::extract(ndvi_stack, review_pixels))
ndvi_long <- data.table(
  pixel = rep(review_pixels, each = length(ndvi_dates)),
  date = rep(ndvi_dates, times = length(review_pixels)),
  ndvi = as.vector(t(ndvi_matrix))
)
ndvi_long[, `:=`(
  year = as.integer(format(date, "%Y")),
  doy = as.integer(format(date, "%j"))
)]
ndvi_climatology <- ndvi_long[, .(
  ndvi_p10 = quantile(ndvi, 0.1, na.rm = TRUE),
  ndvi_median = median(ndvi, na.rm = TRUE),
  ndvi_p90 = quantile(ndvi, 0.9, na.rm = TRUE)
), by = .(pixel, doy)]
ndvi_climatology[, `:=`(
  signal_min = min(ndvi_p10, na.rm = TRUE),
  signal_max = max(ndvi_p90, na.rm = TRUE)
), by = pixel]
ndvi_climatology[, `:=`(
  ndvi_scaled = (ndvi_median - signal_min) / (signal_max - signal_min),
  ndvi_p10_scaled = (ndvi_p10 - signal_min) / (signal_max - signal_min),
  ndvi_p90_scaled = (ndvi_p90 - signal_min) / (signal_max - signal_min)
)]

rain <- as.data.table(read_parquet(monthly_rain_file))[
  pixel %in% review_pixels
]
rain[, doy := round((month - 0.5) * 365 / 12)]
rain[, rain_scaled := {
  span <- max(median_monthly_rain, na.rm = TRUE) - min(median_monthly_rain, na.rm = TRUE)
  if (is.finite(span) && span > 0) {
    (median_monthly_rain - min(median_monthly_rain, na.rm = TRUE)) / span
  } else rep(0, .N)
}, by = pixel]

pheno <- as.data.table(read_parquet(
  phenology_file,
  col_select = c("pixel", "flag", "Greenup", "Senescence")
))[pixel %in% review_pixels]
parts <- tstrsplit(pheno$flag, "_")
pheno[, `:=`(
  year = as.integer(parts[[1]]),
  raw_season = as.integer(parts[[2]])
)]
event_climatology <- pheno[
  !is.na(Greenup) & raw_season <= config$phenology$max_seasons,
  .(greenup_median_doy = as.numeric(median(as.integer(format(Greenup, "%j")), na.rm = TRUE))),
  by = .(pixel, raw_season)
]

windows <- as.data.table(read_parquet(windows_file))[pixel %in% review_pixels]
windows[, `:=`(
  window_start_doy = round((window_start_month - 1) * 365 / 12 + 1),
  window_end_doy = round(window_end_month * 365 / 12)
)]
windows_nowrap <- windows[window_start_doy <= window_end_doy]
windows_wrap <- windows[window_start_doy > window_end_doy]
windows <- rbindlist(list(
  windows_nowrap,
  windows_wrap[, .SD][, window_end_doy := 365],
  windows_wrap[, .SD][, window_start_doy := 1]
), use.names = TRUE)

annual <- as.data.table(read_parquet(detectability_file))[
  pixel %in% review_pixels
]
annual <- annual[review[, .(pixel, review_stratum)], on = .(pixel), allow.cartesian = TRUE]
annual <- annual[raw_season == fifelse(
  review_stratum == "wet_missing_season2", 2L, 1L
)]

ndvi_climatology <- review_key[ndvi_climatology, on = .(pixel)]
rain <- review_key[rain, on = .(pixel)]
event_climatology <- review_key[event_climatology, on = .(pixel)]
windows <- review_key[windows, on = .(pixel)]
annual <- review_key[annual, on = .(pixel, review_stratum)]

strata <- sort(unique(review$review_stratum))
for (stratum in strata) {
  ndvi_s <- ndvi_climatology[review_stratum == stratum]
  rain_s <- rain[review_stratum == stratum]
  event_s <- event_climatology[review_stratum == stratum]
  windows_s <- windows[review_stratum == stratum]
  annual_s <- annual[review_stratum == stratum]

  p_climatology <- ggplot() +
    geom_rect(
      data = windows_s,
      aes(
        xmin = window_start_doy, xmax = window_end_doy,
        ymin = -Inf, ymax = Inf, fill = factor(season_id)
      ),
      alpha = 0.08,
      inherit.aes = FALSE
    ) +
    geom_ribbon(
      data = ndvi_s,
      aes(doy, ymin = ndvi_p10_scaled, ymax = ndvi_p90_scaled),
      fill = "#3B8D5A", alpha = 0.18
    ) +
    geom_line(data = ndvi_s, aes(doy, ndvi_scaled), colour = "#176B3A", linewidth = 0.55) +
    geom_line(
      data = rain_s,
      aes(doy, rain_scaled),
      colour = "#2C7FB8", linewidth = 0.65
    ) +
    geom_vline(
      data = event_s,
      aes(xintercept = greenup_median_doy, colour = factor(raw_season)),
      linewidth = 0.55, linetype = "dashed"
    ) +
    facet_wrap(~facet_label, ncol = 4) +
    scale_x_continuous(
      breaks = c(15, 74, 135, 196, 258, 319),
      labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
      limits = c(1, 365)
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = c("1" = "#66C2A5", "2" = "#FC8D62"), guide = "none") +
    scale_colour_manual(values = c("1" = "#1B9E77", "2" = "#D95F02"), name = "Event season") +
    labs(
      title = paste("Seasonal climatology:", stratum),
      subtitle = "Green = scaled NDVI median (band: p10–p90); blue = scaled monthly rainfall; dashed = median detected greenup",
      x = NULL, y = "Within-pixel scaled signal"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(panel_dir, paste0("climatology_", stratum, ".png")),
    p_climatology, width = 16, height = 11, dpi = 180, bg = "white"
  )

  p_annual <- ggplot(annual_s, aes(year, wet_anomaly)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_line(colour = "grey45", linewidth = 0.4) +
    geom_point(aes(fill = season_detected), shape = 21, size = 1.8, colour = "grey20") +
    facet_wrap(~facet_label, ncol = 4) +
    scale_fill_manual(values = c("TRUE" = "#2CA25F", "FALSE" = "#DE2D26")) +
    labs(
      title = paste("Annual wet anomaly and event detection:", stratum),
      x = NULL, y = "Rainfall anomaly (baseline MAD units)", fill = "Event detected"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(panel_dir, paste0("annual_", stratum, ".png")),
    p_annual, width = 16, height = 11, dpi = 180, bg = "white"
  )
}

review_sheet <- unique(review[, .(
  review_stratum, pixel, admin1_name, x, y, raw_season,
  ndvi_amplitude_median, event_coverage, timing_concentration,
  baseline_rainfall_signal, valley_strength,
  suggested_label = NA_character_,
  reviewer_label = NA_character_,
  reviewer_confidence = NA_character_,
  reviewer_notes = NA_character_,
  reviewer = NA_character_,
  review_date = as.IDate(NA)
)])
fwrite(review_sheet, file.path(report_dir, "review_sheet.csv"), na = "")

failed <- FALSE
finish_run(
  run,
  "success",
  list(report_dir = report_dir, pixels = length(review_pixels), panels = length(strata) * 2L)
)
cat("Wrote review panels under", panel_dir, "\n")
