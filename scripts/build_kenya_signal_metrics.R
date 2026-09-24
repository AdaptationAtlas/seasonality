#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(terra)
})

source("R/functions/project_logging.R")
source("R/functions/circular_utils.R")
source("R/functions/detectability.R")

config <- read_project_config()
run <- new_run_context("kenya-signal-metrics", config$country$iso3)
log_event(run, event = "run_started", message = "Building annual Kenya NDVI signal metrics")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
ndvi_dir <- file.path(data_root, "climate_raw", "glass_ndvi_tif")
pheno_file <- file.path(
  data_root, "climate_derived", "glass_phenology", "countries",
  "KEN_seasonal-phenology.parquet"
)
output_file <- file.path(
  data_root, "climate_derived", "glass_phenology",
  "KEN_annual_signal_metrics.parquet"
)
baseline_output_file <- file.path(
  data_root, "climate_derived", "glass_phenology",
  "KEN_baseline_season_signal_metrics.parquet"
)

ndvi_files <- sort(list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE))
if (!length(ndvi_files)) stop("No GLASS NDVI GeoTIFF files found in ", ndvi_dir)
if (!file.exists(pheno_file)) stop("Kenya phenology file not found: ", pheno_file)

tokens <- sub(".*\\.A([0-9]{7})\\..*", "\\1", basename(ndvi_files))
ndvi_dates <- as.Date(substr(tokens, 1, 4), format = "%Y") + as.integer(substr(tokens, 5, 7)) - 1L
keep <- as.integer(format(ndvi_dates, "%Y")) >= config$baseline$start_year &
  as.integer(format(ndvi_dates, "%Y")) <= config$baseline$end_year
ndvi_files <- ndvi_files[keep]
ndvi_dates <- ndvi_dates[keep]

pheno <- as.data.table(read_parquet(pheno_file))
kenya_pixels <- sort(unique(pheno$pixel))
pixel_admin <- unique(pheno[, .(pixel, admin1_name)])
flag_parts <- tstrsplit(pheno$flag, "_")
pheno[, `:=`(
  year = as.integer(flag_parts[[1]]),
  raw_season = as.integer(flag_parts[[2]]),
  complete_event = !is.na(Greenup) & !is.na(Senescence)
)]

event_metrics <- pheno[
  year >= config$baseline$start_year & year <= config$baseline$end_year,
  .(
    candidate_events = .N,
    complete_events = sum(complete_event),
    detected_seasons = uniqueN(raw_season[complete_event]),
    median_R2 = median(R2, na.rm = TRUE),
    median_NSE = median(NSE, na.rm = TRUE)
  ),
  by = .(pixel, year)
]

r <- rast(ndvi_files)
years <- sort(unique(as.integer(format(ndvi_dates, "%Y"))))

annual <- rbindlist(lapply(years, function(year) {
  layer_index <- which(as.integer(format(ndvi_dates, "%Y")) == year)
  values <- as.matrix(terra::extract(r[[layer_index]], kenya_pixels))
  valid_n <- rowSums(is.finite(values))
  valid_fraction <- valid_n / ncol(values)

  row_quantile <- function(x, probability) {
    if (!any(is.finite(x))) return(NA_real_)
    as.numeric(stats::quantile(x, probability, na.rm = TRUE, names = FALSE))
  }
  q10 <- apply(values, 1, row_quantile, probability = 0.1)
  q50 <- apply(values, 1, row_quantile, probability = 0.5)
  q90 <- apply(values, 1, row_quantile, probability = 0.9)

  log_event(
    run,
    event = "year_processed",
    message = as.character(year),
    data = list(pixels = length(kenya_pixels), layers = length(layer_index))
  )

  data.table(
    pixel = kenya_pixels,
    year = year,
    valid_fraction = valid_fraction,
    ndvi_q10 = q10,
    ndvi_q50 = q50,
    ndvi_q90 = q90,
    ndvi_amplitude = q90 - q10,
    ndvi_mean = rowMeans(values, na.rm = TRUE),
    ndvi_sd = apply(values, 1, stats::sd, na.rm = TRUE),
    sufficient_data = valid_fraction >= 0.75
  )
}), use.names = TRUE)

annual <- event_metrics[annual, on = .(pixel, year)]
annual <- pixel_admin[annual, on = .(pixel)]
setcolorder(annual, c("pixel", "admin1_name", "year"))
write_parquet(annual, output_file)

baseline_years <- config$baseline$end_year - config$baseline$start_year + 1L
pixel_signal <- annual[, .(
  years_with_ndvi = sum(sufficient_data),
  ndvi_amplitude_median = median(ndvi_amplitude, na.rm = TRUE),
  ndvi_amplitude_p10 = quantile(ndvi_amplitude, 0.1, na.rm = TRUE),
  ndvi_amplitude_p90 = quantile(ndvi_amplitude, 0.9, na.rm = TRUE),
  any_event_year_coverage = mean(!is.na(complete_events) & complete_events > 0),
  bimodal_year_coverage = mean(!is.na(detected_seasons) & detected_seasons >= 2)
), by = .(pixel, admin1_name)]

pheno_complete <- pheno[
  complete_event &
    year >= config$baseline$start_year & year <= config$baseline$end_year
]
pheno_complete[, greenup_doy := as.integer(format(Greenup, "%j"))]

baseline_season <- pheno_complete[, .(
  observed_years = uniqueN(year),
  event_coverage = uniqueN(year) / baseline_years,
  greenup_median_doy = as.numeric(quantile_circular_safe(greenup_doy, 0.5)),
  timing_concentration = circular_concentration(greenup_doy),
  median_R2 = median(R2, na.rm = TRUE),
  median_NSE = median(NSE, na.rm = TRUE)
), by = .(pixel, raw_season)]

baseline_season <- pixel_signal[baseline_season, on = .(pixel)]
setcolorder(baseline_season, c("pixel", "admin1_name", "raw_season"))
write_parquet(baseline_season, baseline_output_file)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    annual_output = output_file,
    baseline_output = baseline_output_file,
    annual_rows = nrow(annual),
    baseline_rows = nrow(baseline_season),
    pixels = uniqueN(annual$pixel)
  )
)
cat("Wrote", output_file, "\n")
cat("Wrote", baseline_output_file, "\n")
