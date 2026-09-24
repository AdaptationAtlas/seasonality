#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

source("R/functions/project_logging.R")
source("R/functions/circular_utils.R")

config <- read_project_config()
run <- new_run_context("score-kenya-ceepa-validation", config$country$iso3)
log_event(run, event = "run_started", message = "Comparing CEEPA planting dates with stable GLASS events")
ceepa_config <- config$validation$ceepa

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
report_dir <- file.path(
  data_root, "models", "validation_reports", "kenya_detectability", "independent_validation"
)
observations_file <- file.path(report_dir, "KEN_ceepa_crop_observations.parquet")
stable_file <- file.path(base_dir, "KEN_stable_phenology_events.parquet")
if (!all(file.exists(c(observations_file, stable_file)))) stop("Missing CEEPA scoring inputs.")

observed <- as.data.table(read_parquet(observations_file))[
  planting_validation_eligible == TRUE & !is.na(survey_season_name)
]
observed[, season_id := fifelse(survey_season_name == "long rains", 1L, 2L)]
observed_summary <- observed[, {
  planting_q <- quantile_circular_safe(planting_doy, c(0.1, 0.5, 0.9))
  length_values <- season_length_days[!is.na(season_length_days)]
  list(
    observed_records = .N,
    observed_households = uniqueN(hhcode),
    observed_planting_q10 = as.numeric(planting_q[1]),
    observed_planting_median = as.numeric(planting_q[2]),
    observed_planting_q90 = as.numeric(planting_q[3]),
    observed_exact_share = mean(planting_precision == "day"),
    observed_length_records = length(length_values),
    observed_length_median = if (length(length_values)) as.numeric(median(length_values)) else NA_real_
  )
}, by = .(admin1_name, season_id, survey_season_name)]

stable <- as.data.table(read_parquet(stable_file))[
  year == ceepa_config$comparison_year & has_mapped_crop == TRUE
]
remote_summary <- stable[, {
  detected_doy <- as.integer(format(Greenup[season_detected], "%j"))
  quality_doy <- as.integer(format(Greenup[season_detected & fit_pass], "%j"))
  length_values <- season_length_days[season_detected & fit_pass & !is.na(season_length_days)]
  q <- quantile_circular_safe(quality_doy, c(0.1, 0.5, 0.9))
  list(
    remote_candidate_pixels = .N,
    remote_detected_pixels = sum(season_detected),
    remote_quality_pixels = sum(season_detected & fit_pass, na.rm = TRUE),
    remote_detection_coverage = mean(season_detected),
    remote_quality_coverage = mean(season_detected & fit_pass, na.rm = TRUE),
    remote_greenup_q10 = as.numeric(q[1]),
    remote_greenup_median = as.numeric(q[2]),
    remote_greenup_q90 = as.numeric(q[3]),
    remote_length_pixels = length(length_values),
    remote_length_median = if (length(length_values)) as.numeric(median(length_values)) else NA_real_,
    remote_detected_doy_median = if (length(detected_doy)) as.numeric(quantile_circular_safe(detected_doy, 0.5)) else NA_real_
  )
}, by = .(admin1_name, season_id)]

comparison <- merge(
  observed_summary,
  remote_summary,
  by = c("admin1_name", "season_id"),
  all.x = TRUE
)
comparison[, `:=`(
  greenup_minus_planting_days = signed_circ_diff(
    remote_greenup_median, observed_planting_median
  ),
  absolute_timing_difference_days = circ_dist(
    remote_greenup_median, observed_planting_median
  ),
  remote_minus_observed_length_days = remote_length_median - observed_length_median
)]
setorder(comparison, admin1_name, season_id)
fwrite(comparison, file.path(report_dir, "KEN_ceepa_remote_validation.csv"), na = "")
write_parquet(comparison, file.path(report_dir, "KEN_ceepa_remote_validation.parquet"))

score_summary <- comparison[
  observed_records >= ceepa_config$min_observed_records &
    remote_quality_pixels >= ceepa_config$min_remote_pixels,
  .(
    county_seasons = .N,
    median_signed_greenup_lag_days = median(greenup_minus_planting_days, na.rm = TRUE),
    median_absolute_timing_difference_days = median(absolute_timing_difference_days, na.rm = TRUE),
    median_remote_quality_coverage = median(remote_quality_coverage, na.rm = TRUE),
    median_length_difference_days = median(remote_minus_observed_length_days, na.rm = TRUE)
  )
]
fwrite(score_summary, file.path(report_dir, "KEN_ceepa_validation_score_summary.csv"), na = "")

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    comparison_rows = nrow(comparison),
    scored_county_seasons = score_summary$county_seasons,
    report_dir = report_dir
  )
)
cat("Wrote CEEPA comparison under", report_dir, "\n")
