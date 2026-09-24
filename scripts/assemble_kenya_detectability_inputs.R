#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("kenya-detectability-inputs", config$country$iso3)
log_event(run, event = "run_started", message = "Assembling Kenya detectability calibration table")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")

ndvi_file <- file.path(base_dir, "KEN_annual_signal_metrics.parquet")
rain_file <- file.path(base_dir, "KEN_annual_rainfall_metrics.parquet")
baseline_season_file <- file.path(base_dir, "KEN_baseline_season_signal_metrics.parquet")
pheno_file <- file.path(base_dir, "countries", "KEN_seasonal-phenology.parquet")
environment_file <- file.path(base_dir, "pixel_index_aridity_elevation.parquet")
landcover_file <- file.path(base_dir, "pixel_index_lulc-commodities.parquet")
output_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")

required_files <- c(
  ndvi_file, rain_file, baseline_season_file,
  pheno_file, environment_file, landcover_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = ", "))

ndvi <- as.data.table(read_parquet(ndvi_file))
rain <- as.data.table(read_parquet(rain_file))
rain[, admin1_name := NULL]
annual <- merge(ndvi, rain, by = c("pixel", "year"), all.x = TRUE)
annual[, rainfall_signal_strength := pmax(rain_h1, rain_h2, na.rm = TRUE)]
annual[!is.finite(rainfall_signal_strength), rainfall_signal_strength := NA_real_]

pixel_rain <- annual[, .(
  baseline_annual_rain = median(annual_rain, na.rm = TRUE),
  baseline_rainfall_si = median(rainfall_si, na.rm = TRUE),
  baseline_rain_h1 = median(rain_h1, na.rm = TRUE),
  baseline_rain_h2 = median(rain_h2, na.rm = TRUE),
  baseline_rainfall_signal = median(rainfall_signal_strength, na.rm = TRUE)
), by = pixel]

environment <- as.data.table(read_parquet(environment_file))
if (!"aridity_bin" %in% names(environment)) {
  environment[, aridity_bin := cut(
    aridity,
    breaks = c(-Inf, 0.05, 0.20, 0.50, 0.65, Inf),
    labels = c("hyper-arid", "arid", "semi-arid", "sub-humid", "humid")
  )]
}
environment <- environment[, .(pixel, aridity, aridity_bin, elevation, elev_bin)]

landcover <- as.data.table(read_parquet(landcover_file))
livestock_columns <- c("buffalo", "chickens", "cattle", "goats", "pigs", "sheep")
landcover_columns <- c("code", "pixel", "short_name", "description")
crop_columns <- setdiff(names(landcover), c(landcover_columns, livestock_columns))
landcover[, crop_area_total := rowSums(.SD, na.rm = TRUE), .SDcols = crop_columns]
landcover[, has_mapped_crop := crop_area_total > 0]
landcover <- landcover[, .(
  pixel, code, short_name, description, crop_area_total, has_mapped_crop
)]
setnames(
  landcover,
  c("code", "short_name", "description"),
  c("landcover_code", "landcover_short", "landcover_description")
)

baseline_season <- as.data.table(read_parquet(baseline_season_file))
baseline_season <- baseline_season[raw_season <= config$phenology$max_seasons]
baseline_season <- pixel_rain[baseline_season, on = .(pixel)]
baseline_season <- environment[baseline_season, on = .(pixel)]
baseline_season <- landcover[baseline_season, on = .(pixel)]

pheno <- as.data.table(read_parquet(
  pheno_file,
  col_select = c(
    "pixel", "flag", "Greenup", "Senescence", "R2", "NSE", "RMSE"
  )
))
flag_parts <- tstrsplit(pheno$flag, "_")
pheno[, `:=`(
  year = as.integer(flag_parts[[1]]),
  raw_season = as.integer(flag_parts[[2]]),
  complete_event = !is.na(Greenup) & !is.na(Senescence)
)]

event_year <- pheno[
  year >= config$baseline$start_year &
    year <= config$baseline$end_year &
    raw_season <= config$phenology$max_seasons,
  .(
    season_detected = any(complete_event),
    greenup_doy = if (any(complete_event)) {
      as.numeric(median(as.integer(format(Greenup[complete_event], "%j")), na.rm = TRUE))
    } else NA_real_,
    senescence_doy = if (any(complete_event)) {
      as.numeric(median(as.integer(format(Senescence[complete_event], "%j")), na.rm = TRUE))
    } else NA_real_,
    event_R2 = median(R2[complete_event], na.rm = TRUE),
    event_NSE = median(NSE[complete_event], na.rm = TRUE),
    event_RMSE = median(RMSE[complete_event], na.rm = TRUE)
  ),
  by = .(pixel, year, raw_season)
]

inputs <- baseline_season[
  annual,
  on = .(pixel),
  allow.cartesian = TRUE,
  nomatch = 0L
]
inputs <- event_year[
  inputs,
  on = .(pixel, year, raw_season)
]
inputs[is.na(season_detected), season_detected := FALSE]

leading <- c(
  "pixel", "admin1_name", "year", "raw_season",
  "aridity_bin", "elev_bin", "landcover_short", "crop_area_total", "has_mapped_crop",
  "season_detected", "greenup_doy", "senescence_doy"
)
setcolorder(inputs, c(leading, setdiff(names(inputs), leading)))
setorder(inputs, pixel, raw_season, year)
write_parquet(inputs, output_file)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    output = output_file,
    rows = nrow(inputs),
    pixels = uniqueN(inputs$pixel),
    admin1_units = uniqueN(inputs$admin1_name)
  )
)
cat("Wrote", output_file, "\n")
