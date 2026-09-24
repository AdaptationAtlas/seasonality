#!/usr/bin/env Rscript

source("R/functions/circular_utils.R")
source("R/functions/season_assignment.R")
source("R/functions/project_logging.R")
source("R/functions/detectability.R")

assert_equal <- function(actual, expected, tolerance = 1e-8) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop(
      "Expected ", paste(expected, collapse = ", "),
      "; got ", paste(actual, collapse = ", "),
      call. = FALSE
    )
  }
}

assert_equal(season_length_doy(350, 20), 35)
assert_equal(circ_dist(360, 5), 10)
assert_equal(forward_circ_dist(360, 5), 10)
assert_equal(in_circular_window(c(355, 5, 100), 350, 20), c(TRUE, TRUE, FALSE))

q <- quantile_circular_safe(c(350, 355, 5, 10), probs = 0.5)
if (circ_dist(q, 365) > 10) stop("Circular median failed near year boundary.")

tmp <- tempfile("seasonality-test-logs-")
run <- new_run_context("unit-test", log_dir = tmp)
log_event(run, event = "test", message = "logging works", data = list(value = 1))
finish_run(run)
lines <- readLines(run$log_file)
if (length(lines) != 2L) stop("Expected two JSONL log records.")
invisible(lapply(lines, jsonlite::fromJSON))

dates <- seq(as.Date("2001-01-01"), as.Date("2001-12-31"), by = "8 days")
seasonal <- 0.5 + 0.2 * sin(2 * pi * as.numeric(format(dates, "%j")) / 365)
metrics <- annual_ndvi_metrics(seasonal, dates)
if (!metrics$sufficient_data || metrics$ndvi_amplitude < 0.3) {
  stop("Seasonal NDVI metrics failed.")
}

thresholds <- list(
  ndvi_amplitude_min = 0.1,
  timing_concentration_min = 0.7,
  event_coverage_min = 0.6,
  rainfall_seasonality_min = 0.5,
  wet_anomaly_min = 1.5
)
paths <- classify_signal_path(
  sufficient_data = c(TRUE, TRUE, TRUE, FALSE),
  ndvi_amplitude = c(0.2, 0.02, 0.2, 0.2),
  timing_concentration = c(0.9, 0.2, 0.9, 0.9),
  event_coverage = c(0.8, 0.2, 0.8, 0.8),
  rainfall_seasonality = c(0.8, 0.8, 0.8, 0.8),
  expected_seasons = c(2, 2, 2, 2),
  detected_seasons = c(2, 2, 1, 2),
  wet_anomaly = c(0, 0, 2, 0),
  thresholds = thresholds
)
assert_equal(
  as.character(paths),
  c("ndvi_seasonal", "rainfall_proxy", "wet_merged", "insufficient_data")
)

cat("All tests passed.\n")
