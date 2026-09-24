#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
})

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("kenya-calibration-report", config$country$iso3)
log_event(run, event = "run_started", message = "Building Kenya detectability calibration report")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
report_dir <- file.path(data_root, "models", "validation_reports", "kenya_detectability")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(base_dir, "KEN_detectability_inputs.parquet")
windows_file <- file.path(base_dir, "KEN_candidate_season_windows.parquet")
coords_file <- file.path(base_dir, "pixel_index.parquet")
if (!all(file.exists(c(input_file, windows_file, coords_file)))) {
  stop("Run detectability and season-window builders before calibration report.")
}

dat <- as.data.table(read_parquet(input_file))
windows <- as.data.table(read_parquet(windows_file))
coords <- as.data.table(read_parquet(coords_file))[, .(pixel, x, y)]
valley_by_pixel <- unique(windows[, .(pixel, valley_strength)])
dat <- valley_by_pixel[dat, on = .(pixel)]

baseline_fields <- c(
  "pixel", "admin1_name", "raw_season", "aridity_bin", "elev_bin",
  "landcover_short", "crop_area_total", "has_mapped_crop",
  "ndvi_amplitude_median", "event_coverage", "timing_concentration",
  "baseline_rainfall_si", "baseline_rain_h1", "baseline_rain_h2",
  "baseline_rainfall_signal"
)
baseline <- unique(dat[, ..baseline_fields])
window_pixel <- windows[, .(
  valley_strength = valley_strength[1],
  peak_1 = peak_month[season_id == 1L][1],
  peak_2 = peak_month[season_id == 2L][1]
), by = pixel]
baseline <- window_pixel[baseline, on = .(pixel)]
baseline <- coords[baseline, on = .(pixel)]

pixel_summary <- dcast(
  baseline,
  pixel + admin1_name + aridity_bin + elev_bin + landcover_short +
    crop_area_total + has_mapped_crop + x + y + valley_strength + peak_1 + peak_2 ~ raw_season,
  value.var = c("event_coverage", "timing_concentration"),
  fill = NA_real_
)
pixel_signal <- unique(baseline[, .(
  pixel,
  ndvi_amplitude_median,
  baseline_rainfall_si,
  baseline_rain_h1,
  baseline_rain_h2,
  baseline_rainfall_signal
)])
pixel_summary <- pixel_signal[pixel_summary, on = .(pixel)]
write_parquet(pixel_summary, file.path(report_dir, "pixel_signal_summary.parquet"))

quantile_probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
quantiles <- baseline[, {
  values <- list(
    ndvi_amplitude = ndvi_amplitude_median,
    event_coverage = event_coverage,
    timing_concentration = timing_concentration,
    rainfall_signal = baseline_rainfall_signal,
    valley_strength = valley_strength
  )
  rbindlist(lapply(names(values), function(metric_name) {
    data.table(
      metric = metric_name,
      probability = quantile_probs,
      value = as.numeric(quantile(values[[metric_name]], quantile_probs, na.rm = TRUE))
    )
  }))
}, by = .(has_mapped_crop, aridity_bin, raw_season)]
fwrite(quantiles, file.path(report_dir, "calibration_quantiles.csv"))

crop_s1 <- baseline[has_mapped_crop == TRUE & raw_season == 1L]
thresholds <- list(
  amplitude_low = as.numeric(quantile(crop_s1$ndvi_amplitude_median, 0.25, na.rm = TRUE)),
  amplitude_high = as.numeric(quantile(crop_s1$ndvi_amplitude_median, 0.75, na.rm = TRUE)),
  coverage_high = as.numeric(quantile(crop_s1$event_coverage, 0.75, na.rm = TRUE)),
  concentration_high = as.numeric(quantile(crop_s1$timing_concentration, 0.75, na.rm = TRUE)),
  rainfall_signal_low = as.numeric(quantile(crop_s1$baseline_rainfall_signal, 0.25, na.rm = TRUE)),
  rainfall_signal_high = as.numeric(quantile(crop_s1$baseline_rainfall_signal, 0.75, na.rm = TRUE)),
  valley_high = as.numeric(quantile(crop_s1$valley_strength, 0.75, na.rm = TRUE))
)

crop_s2_year <- dat[has_mapped_crop == TRUE & raw_season == 2L]
thresholds$wet_anomaly_high <- as.numeric(quantile(crop_s2_year$wet_anomaly, 0.9, na.rm = TRUE))
jsonlite::write_json(
  thresholds,
  file.path(report_dir, "review_sampling_thresholds.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

wet_merge_candidates <- crop_s2_year[
  !season_detected &
    wet_anomaly >= thresholds$wet_anomaly_high &
    valley_strength >= thresholds$valley_high
]
wet_merge_candidates <- coords[wet_merge_candidates, on = .(pixel)]
write_parquet(
  wet_merge_candidates,
  file.path(report_dir, "wet_merge_candidates.parquet")
)

candidate_sets <- list(
  strong_ndvi = baseline[
    has_mapped_crop == TRUE & raw_season == 1L &
      ndvi_amplitude_median >= thresholds$amplitude_high &
      event_coverage >= thresholds$coverage_high &
      timing_concentration >= thresholds$concentration_high
  ],
  weak_ndvi_rain_seasonal = baseline[
    has_mapped_crop == TRUE & raw_season == 1L &
      ndvi_amplitude_median <= thresholds$amplitude_low &
      baseline_rainfall_signal >= thresholds$rainfall_signal_high
  ],
  humid_weak = baseline[
    has_mapped_crop == TRUE & raw_season == 1L & aridity_bin == "humid" &
      ndvi_amplitude_median <= thresholds$amplitude_low &
      baseline_rainfall_signal <= thresholds$rainfall_signal_low
  ],
  wet_missing_season2 = wet_merge_candidates
)

set.seed(20260924)
review_sample <- rbindlist(lapply(names(candidate_sets), function(stratum) {
  candidate <- copy(candidate_sets[[stratum]])
  if (!nrow(candidate)) return(NULL)
  candidate <- candidate[sample(.N)]
  candidate <- candidate[!duplicated(pixel)]
  candidate <- candidate[seq_len(min(12L, .N))]
  candidate[, review_stratum := stratum]
  candidate
}), fill = TRUE)
fwrite(review_sample, file.path(report_dir, "review_sample.csv"))

plot_base <- baseline[has_mapped_crop == TRUE]
if (nrow(plot_base) > 60000L) plot_base <- plot_base[sample(.N, 60000L)]

p_signal <- ggplot(
  plot_base,
  aes(ndvi_amplitude_median, event_coverage, colour = baseline_rainfall_signal)
) +
  geom_point(alpha = 0.18, size = 0.6) +
  facet_grid(raw_season ~ aridity_bin) +
  scale_colour_viridis_c(option = "C", na.value = "grey70") +
  labs(
    title = "Kenya crop-pixel seasonal signal diagnostics",
    x = "Median annual NDVI amplitude",
    y = "Phenology event-year coverage",
    colour = "Rainfall signal"
  ) +
  theme_bw(base_size = 10)
ggsave(
  file.path(report_dir, "ndvi_rainfall_signal.png"),
  p_signal, width = 13, height = 7, dpi = 180, bg = "white"
)

p_windows <- ggplot(
  pixel_summary[has_mapped_crop == TRUE],
  aes(valley_strength, event_coverage_2, colour = aridity_bin)
) +
  geom_point(alpha = 0.22, size = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
  labs(
    title = "Second-season evidence versus rainfall dry-valley strength",
    x = "Baseline dry-valley strength",
    y = "Season-2 event-year coverage",
    colour = "Aridity class"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(report_dir, "season2_valley_evidence.png"),
  p_windows, width = 10, height = 7, dpi = 180, bg = "white"
)

p_sample <- ggplot(
  review_sample,
  aes(x, y, colour = review_stratum)
) +
  geom_point(size = 2, alpha = 0.85) +
  coord_equal() +
  labs(
    title = "Stratified Kenya calibration sample",
    x = "Longitude", y = "Latitude", colour = "Review stratum"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(report_dir, "review_sample_map.png"),
  p_sample, width = 8, height = 8, dpi = 180, bg = "white"
)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    report_dir = report_dir,
    review_rows = nrow(review_sample),
    wet_merge_candidates = nrow(wet_merge_candidates)
  )
)
cat("Wrote calibration report assets under", report_dir, "\n")
