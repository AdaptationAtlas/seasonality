#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")
source("R/functions/season_windows.R")

config <- read_project_config()
run <- new_run_context("kenya-season-windows", config$country$iso3)
log_event(run, event = "run_started", message = "Building candidate baseline season windows")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
rainfall_file <- file.path(base_dir, "KEN_baseline_monthly_rainfall.parquet")
annual_rain_file <- file.path(base_dir, "KEN_annual_rainfall_metrics.parquet")
phenology_file <- file.path(base_dir, "countries", "KEN_seasonal-phenology.parquet")
environment_file <- file.path(base_dir, "pixel_index_aridity_elevation.parquet")
landcover_file <- file.path(base_dir, "pixel_index_lulc-commodities.parquet")
output_file <- file.path(base_dir, "KEN_candidate_season_windows.parquet")

if (!file.exists(rainfall_file)) stop("Missing rainfall baseline: ", rainfall_file)
required <- c(annual_rain_file, phenology_file, environment_file, landcover_file)
if (!all(file.exists(required))) stop("Missing season-window metadata inputs.")

rainfall <- as.data.table(read_parquet(rainfall_file))
setorder(rainfall, pixel, month)

windows <- rainfall[, {
  monthly <- median_monthly_rain[match(1:12, month)]
  derive_two_season_windows(
    monthly,
    min_peak_separation = config$season_windows$min_peak_separation_months
  )
}, by = pixel]

annual_rain <- as.data.table(read_parquet(annual_rain_file))
annual_rain[, rainfall_signal_strength := pmax(rain_h1, rain_h2, na.rm = TRUE)]
annual_rain[!is.finite(rainfall_signal_strength), rainfall_signal_strength := NA_real_]
pixel_rain <- annual_rain[, .(
  baseline_rainfall_si = median(rainfall_si, na.rm = TRUE),
  baseline_rain_h1 = median(rain_h1, na.rm = TRUE),
  baseline_rain_h2 = median(rain_h2, na.rm = TRUE),
  baseline_rainfall_signal = median(rainfall_signal_strength, na.rm = TRUE)
), by = pixel]

admin1 <- unique(as.data.table(read_parquet(
  phenology_file,
  col_select = c("pixel", "admin1_name")
)))
environment <- as.data.table(read_parquet(environment_file))[, .(
  pixel, aridity_bin, elev_bin
)]
landcover <- as.data.table(read_parquet(landcover_file))
livestock_columns <- c("buffalo", "chickens", "cattle", "goats", "pigs", "sheep")
noncrop_columns <- c("code", "pixel", "short_name", "description", livestock_columns)
crop_columns <- setdiff(names(landcover), noncrop_columns)
landcover[, crop_area_total := rowSums(.SD, na.rm = TRUE), .SDcols = crop_columns]
landcover[, has_mapped_crop := crop_area_total > 0]
landcover <- landcover[, .(
  pixel, landcover_short = short_name, crop_area_total, has_mapped_crop
)]

metadata <- pixel_rain[admin1, on = "pixel"]
metadata <- environment[metadata, on = "pixel"]
metadata <- landcover[metadata, on = "pixel"]
windows <- metadata[windows, on = .(pixel), nomatch = 0L]
setcolorder(windows, c(
  "pixel", "admin1_name", "season_id", "peak_month",
  "window_start_month", "window_end_month"
))
setorder(windows, pixel, season_id)
write_parquet(windows, output_file)

failed <- FALSE
finish_run(
  run,
  "success",
  list(output = output_file, rows = nrow(windows), pixels = uniqueN(windows$pixel))
)
cat("Wrote", output_file, "\n")
