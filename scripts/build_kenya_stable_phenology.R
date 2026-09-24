#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")
source("R/functions/season_windows.R")

config <- read_project_config()
run <- new_run_context("kenya-stable-phenology", config$country$iso3)
log_event(run, event = "run_started", message = "Assigning fitted events to stable rainfall windows")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
phenology_file <- file.path(base_dir, "countries", "KEN_seasonal-phenology.parquet")
windows_file <- file.path(base_dir, "KEN_candidate_season_windows.parquet")
output_file <- file.path(base_dir, "KEN_stable_phenology_events.parquet")
if (!all(file.exists(c(phenology_file, windows_file)))) stop("Missing stable-phenology inputs.")

events <- as.data.table(read_parquet(
  phenology_file,
  col_select = c(
    "pixel", "flag", "origin", "Greenup", "Maturity", "Senescence", "Dormancy",
    "meth", "R2", "NSE", "RMSE", "pvalue"
  )
))
flag_parts <- tstrsplit(events$flag, "_", fixed = TRUE)
events[, `:=`(
  year = as.integer(flag_parts[[1]]),
  raw_season = as.integer(flag_parts[[2]]),
  greenup_month = as.integer(format(Greenup, "%m"))
)]
events <- events[
  year >= config$baseline$start_year & year <= config$baseline$end_year & !is.na(Greenup)
]

windows <- as.data.table(read_parquet(windows_file))
window_keys <- c(
  "pixel", "admin1_name", "season_id", "peak_month", "window_start_month",
  "window_end_month", "aridity_bin", "elev_bin", "landcover_short",
  "crop_area_total", "has_mapped_crop", "valley_strength"
)
windows <- windows[, ..window_keys]

candidates <- windows[events, on = "pixel", allow.cartesian = TRUE, nomatch = 0L]
candidates <- candidates[in_circular_month_window(
  greenup_month, window_start_month, window_end_month
)]
candidates[, `:=`(
  complete_pair = !is.na(Greenup) & !is.na(Senescence),
  fit_pass = !is.na(R2) & !is.na(NSE) &
    R2 >= config$phenology$ndvi_r2_min & NSE >= config$phenology$ndvi_nse_min,
  peak_distance_months = circular_month_distance(greenup_month, peak_month)
)]
candidates[, event_candidate_count := .N, by = .(pixel, year, season_id)]
setorder(
  candidates,
  pixel, year, season_id,
  -complete_pair, -fit_pass, peak_distance_months, -NSE, -R2
)
selected <- unique(candidates, by = c("pixel", "year", "season_id"))

grid <- windows[, .(
  year = config$baseline$start_year:config$baseline$end_year
), by = window_keys]
stable <- merge(
  grid,
  selected[, .(
    pixel, year, season_id, raw_season, flag, origin,
    Greenup, Maturity, Senescence, Dormancy, meth, R2, NSE, RMSE, pvalue,
    complete_pair, fit_pass, peak_distance_months, event_candidate_count
  )],
  by = c("pixel", "year", "season_id"),
  all.x = TRUE,
  sort = FALSE
)
stable[, `:=`(
  season_detected = !is.na(Greenup),
  event_candidate_count = fifelse(is.na(event_candidate_count), 0L, event_candidate_count),
  season_length_days = as.integer(Senescence - Greenup)
)]
stable[season_length_days < 0L | season_length_days > 365L, season_length_days := NA_integer_]
setorder(stable, pixel, year, season_id)
write_parquet(stable, output_file)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    output = output_file,
    rows = nrow(stable),
    pixels = uniqueN(stable$pixel),
    detected = sum(stable$season_detected),
    duplicate_groups = stable[event_candidate_count > 1L, .N]
  )
)
cat("Wrote", output_file, "\n")
