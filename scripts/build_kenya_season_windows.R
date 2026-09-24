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
detectability_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")
output_file <- file.path(base_dir, "KEN_candidate_season_windows.parquet")

if (!file.exists(rainfall_file)) stop("Missing rainfall baseline: ", rainfall_file)
if (!file.exists(detectability_file)) stop("Missing detectability inputs: ", detectability_file)

rainfall <- as.data.table(read_parquet(rainfall_file))
setorder(rainfall, pixel, month)

windows <- rainfall[, {
  monthly <- median_monthly_rain[match(1:12, month)]
  derive_two_season_windows(
    monthly,
    min_peak_separation = config$season_windows$min_peak_separation_months
  )
}, by = pixel]

metadata <- unique(as.data.table(read_parquet(
  detectability_file,
  col_select = c(
    "pixel", "admin1_name", "aridity_bin", "elev_bin",
    "landcover_short", "crop_area_total", "has_mapped_crop",
    "baseline_rainfall_si", "baseline_rain_h1", "baseline_rain_h2",
    "baseline_rainfall_signal"
  )
)))
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
