#!/usr/bin/env Rscript

source("R/functions/circular_utils.R")
source("R/functions/season_assignment.R")
source("R/functions/project_logging.R")
source("R/functions/detectability.R")
source("R/functions/rainfall_seasonality.R")
source("R/functions/season_windows.R")
source("R/functions/ndvi_dates.R")
source("R/functions/ceepa_validation.R")

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

uniform_rain <- rep(100, 12)
uniform_metrics <- rainfall_seasonality_metrics(uniform_rain)
assert_equal(uniform_metrics$rainfall_si, 0)
assert_equal(uniform_metrics$rain_h1, 0, tolerance = 1e-12)
assert_equal(uniform_metrics$rain_h2, 0, tolerance = 1e-12)

unimodal_rain <- rep(0, 12)
unimodal_rain[4] <- 100
unimodal_metrics <- rainfall_seasonality_metrics(unimodal_rain)
assert_equal(unimodal_metrics$rain_h1, 1)
assert_equal(unimodal_metrics$rain_h2, 1)
assert_equal(unimodal_metrics$wettest_month, 4)

bimodal_rain <- rep(0, 12)
bimodal_rain[c(4, 10)] <- 100
bimodal_metrics <- rainfall_seasonality_metrics(bimodal_rain)
assert_equal(bimodal_metrics$rain_h1, 0, tolerance = 1e-12)
assert_equal(bimodal_metrics$rain_h2, 1)

window_rain <- c(5, 10, 40, 100, 60, 15, 5, 10, 30, 80, 50, 10)
windows <- derive_two_season_windows(window_rain, min_peak_separation = 3L)
assert_equal(windows$peak_month, c(4, 10))
if (nrow(windows) != 2L || any(windows$valley_strength <= 0)) {
  stop("Two-season rainfall windows failed.")
}

flat_windows <- derive_two_season_windows(rep(100, 12))
if (nrow(flat_windows) != 0L) stop("Flat rainfall should not produce two peaks.")

glass_dates <- parse_glass_ndvi_dates(c(
  "GLASS13B01.V10.A2000001.2023068.tif",
  "GLASS13B01.V10.A2000365.2023068.tif"
))
assert_equal(
  as.character(glass_dates),
  c("2000-01-01", "2000-12-30")
)

ceepa_dates <- parse_ceepa_calendar_dates(c(
  "15mar03", "ddapr03", "w2oct03", "ddsept03", "continous", "ddmon90", NA
))
assert_equal(
  ceepa_dates$doy,
  c(74L, 105L, 284L, 258L, NA_integer_, NA_integer_, NA_integer_)
)
assert_equal(
  ceepa_dates$precision,
  c("day", "month", "week", "month", NA_character_, NA_character_, NA_character_)
)
assert_equal(
  ceepa_dates$parse_status,
  c("parsed", "parsed", "parsed", "parsed", "continuous", "unparsed", "missing")
)
assert_equal(
  ceepa_current_county(c("TRNASNZOIA", "meru north", "NITHI")),
  c("Trans Nzoia", "Meru", "Tharaka-Nithi")
)

cat("All tests passed.\n")
