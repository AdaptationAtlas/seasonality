#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")
source("R/functions/circular_utils.R")
source("R/functions/detectability.R")

config <- read_project_config()
run <- new_run_context("kenya-detectability-inputs", config$country$iso3)
log_event(run, event = "run_started", message = "Assembling stable-season Kenya detectability table")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
ndvi_file <- file.path(base_dir, "KEN_annual_signal_metrics.parquet")
rain_file <- file.path(base_dir, "KEN_annual_rainfall_metrics.parquet")
stable_file <- file.path(base_dir, "KEN_stable_phenology_events.parquet")
output_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")

required_files <- c(ndvi_file, rain_file, stable_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = ", "))

ndvi <- as.data.table(read_parquet(ndvi_file))
setnames(
  ndvi,
  c("candidate_events", "complete_events", "detected_seasons", "median_R2", "median_NSE"),
  c("raw_candidate_events", "raw_complete_events", "raw_detected_seasons", "raw_median_R2", "raw_median_NSE")
)
rain <- as.data.table(read_parquet(rain_file))
rain[, admin1_name := NULL]
annual <- merge(ndvi, rain, by = c("pixel", "year"), all.x = TRUE)
annual[, admin1_name := NULL]
annual[, rainfall_signal_strength := pmax(rain_h1, rain_h2, na.rm = TRUE)]
annual[!is.finite(rainfall_signal_strength), rainfall_signal_strength := NA_real_]

pixel_ndvi <- annual[, .(
  years_with_ndvi = sum(sufficient_data, na.rm = TRUE),
  ndvi_amplitude_median = median(ndvi_amplitude[sufficient_data], na.rm = TRUE),
  ndvi_amplitude_p10 = as.numeric(quantile(ndvi_amplitude[sufficient_data], 0.1, na.rm = TRUE)),
  ndvi_amplitude_p90 = as.numeric(quantile(ndvi_amplitude[sufficient_data], 0.9, na.rm = TRUE))
), by = pixel]
pixel_rain <- annual[, .(
  baseline_annual_rain = median(annual_rain, na.rm = TRUE),
  baseline_rainfall_si = median(rainfall_si, na.rm = TRUE),
  baseline_rain_h1 = median(rain_h1, na.rm = TRUE),
  baseline_rain_h2 = median(rain_h2, na.rm = TRUE),
  baseline_rainfall_signal = median(rainfall_signal_strength, na.rm = TRUE)
), by = pixel]

stable <- as.data.table(read_parquet(stable_file))
stable[, `:=`(
  quality_event = season_detected & !is.na(fit_pass) & fit_pass,
  complete_event = season_detected & !is.na(complete_pair) & complete_pair,
  greenup_doy = as.integer(format(Greenup, "%j")),
  senescence_doy = as.integer(format(Senescence, "%j"))
)]

stable_baseline <- stable[, {
  quality_doy <- greenup_doy[quality_event]
  q <- quantile_circular_safe(quality_doy, 0.5)
  list(
    observed_years = .N,
    detected_event_coverage = mean(season_detected),
    complete_event_coverage = mean(complete_event),
    quality_event_coverage = mean(quality_event),
    event_coverage = mean(quality_event),
    greenup_median_doy = as.numeric(q),
    timing_concentration = circular_concentration(quality_doy),
    median_event_R2 = median(R2[quality_event], na.rm = TRUE),
    median_event_NSE = median(NSE[quality_event], na.rm = TRUE)
  )
}, by = .(pixel, season_id)]
stable_baseline[!is.finite(median_event_R2), median_event_R2 := NA_real_]
stable_baseline[!is.finite(median_event_NSE), median_event_NSE := NA_real_]

annual_seasons <- stable[, .(
  detected_seasons = sum(quality_event)
), by = .(pixel, year)]
pixel_events <- annual_seasons[, .(
  any_event_year_coverage = mean(detected_seasons >= 1L),
  bimodal_year_coverage = mean(detected_seasons >= 2L)
), by = pixel]

event_fields <- stable[, .(
  pixel, year, season_id, admin1_name,
  aridity_bin, elev_bin, landcover_short, crop_area_total, has_mapped_crop,
  peak_month, window_start_month, window_end_month, valley_strength,
  source_raw_season = raw_season, event_method = meth,
  Greenup, Maturity, Senescence, Dormancy,
  season_detected, complete_event, quality_event,
  greenup_doy, senescence_doy, season_length_days,
  event_R2 = R2, event_NSE = NSE, event_RMSE = RMSE,
  peak_distance_months, event_candidate_count
)]

inputs <- merge(event_fields, annual, by = c("pixel", "year"), all.x = TRUE)
inputs <- stable_baseline[inputs, on = .(pixel, season_id)]
inputs <- pixel_ndvi[inputs, on = "pixel"]
inputs <- pixel_rain[inputs, on = "pixel"]
inputs <- pixel_events[inputs, on = "pixel"]
inputs <- annual_seasons[inputs, on = .(pixel, year)]

leading <- c(
  "pixel", "admin1_name", "year", "season_id",
  "aridity_bin", "elev_bin", "landcover_short", "crop_area_total", "has_mapped_crop",
  "season_detected", "complete_event", "quality_event", "greenup_doy", "senescence_doy",
  "season_length_days", "source_raw_season"
)
setcolorder(inputs, c(leading, setdiff(names(inputs), leading)))
setorder(inputs, pixel, season_id, year)
write_parquet(inputs, output_file)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    output = output_file,
    rows = nrow(inputs),
    pixels = uniqueN(inputs$pixel),
    admin1_units = uniqueN(inputs$admin1_name),
    stable_season_keys = uniqueN(inputs, by = c("pixel", "year", "season_id"))
  )
)
cat("Wrote", output_file, "\n")
