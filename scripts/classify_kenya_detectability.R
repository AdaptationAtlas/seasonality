#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")
source("R/functions/detectability.R")

config <- read_project_config()
run <- new_run_context("classify-kenya-detectability", config$country$iso3)
log_event(run, event = "run_started", message = "Applying provisional Kenya detectability rules")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
report_dir <- file.path(data_root, "models", "validation_reports", "kenya_detectability")
input_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")
annual_output <- file.path(base_dir, "KEN_detectability_classes.parquet")
baseline_output <- file.path(base_dir, "KEN_baseline_detectability.parquet")
if (!file.exists(input_file)) stop("Missing detectability inputs: ", input_file)

dat <- as.data.table(read_parquet(input_file))
baseline_fields <- c(
  "pixel", "admin1_name", "season_id", "aridity_bin", "elev_bin",
  "landcover_short", "crop_area_total", "has_mapped_crop",
  "years_with_ndvi", "ndvi_amplitude_median", "event_coverage",
  "timing_concentration", "baseline_rainfall_signal", "valley_strength",
  "bimodal_year_coverage"
)
baseline <- unique(dat[, ..baseline_fields])

classified <- as.data.table(classify_baseline_signal(
  season_id = baseline$season_id,
  aridity_bin = baseline$aridity_bin,
  years_with_ndvi = baseline$years_with_ndvi,
  ndvi_amplitude = baseline$ndvi_amplitude_median,
  event_coverage = baseline$event_coverage,
  timing_concentration = baseline$timing_concentration,
  rainfall_signal = baseline$baseline_rainfall_signal,
  valley_strength = baseline$valley_strength,
  bimodal_year_coverage = baseline$bimodal_year_coverage,
  thresholds = config$detectability
))
baseline <- cbind(baseline, classified)
baseline[, threshold_status := config$detectability$status]
setorder(baseline, pixel, season_id)
write_parquet(baseline, baseline_output)

annual <- baseline[, .(
  pixel, season_id, baseline_pathway, ndvi_evidence,
  rainfall_evidence, bimodal_supported, season_supported, threshold_status
)][dat, on = .(pixel, season_id)]
annual[, annual_pathway := classify_annual_pathway(
  baseline_pathway = baseline_pathway,
  season_id = season_id,
  quality_event = quality_event,
  bimodal_supported = bimodal_supported,
  detected_seasons = detected_seasons,
  wet_anomaly = wet_anomaly,
  wet_anomaly_min = config$detectability$wet_anomaly_min
)]
setorder(annual, pixel, year, season_id)
write_parquet(annual, annual_output)

summary <- baseline[, .(
  pixels = uniqueN(pixel),
  crop_pixels = uniqueN(pixel[has_mapped_crop == TRUE])
), by = .(season_id, aridity_bin, baseline_pathway)]
fwrite(summary, file.path(report_dir, "detectability_class_summary.csv"))

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    annual_output = annual_output,
    baseline_output = baseline_output,
    baseline_rows = nrow(baseline),
    annual_rows = nrow(annual),
    not_identifiable = baseline[baseline_pathway == "not_identifiable", .N]
  )
)
cat("Wrote detectability classes under", base_dir, "\n")
