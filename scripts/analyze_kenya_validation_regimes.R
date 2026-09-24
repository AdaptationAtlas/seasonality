#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(ggplot2)
})

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("kenya-validation-regimes", config$country$iso3)
log_event(run, event = "run_started", message = "Profiling CEEPA timing differences by ecological regime")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
base_dir <- file.path(data_root, "climate_derived", "glass_phenology")
report_dir <- file.path(
  data_root, "models", "validation_reports", "kenya_detectability", "independent_validation"
)
comparison_file <- file.path(report_dir, "KEN_ceepa_remote_validation.csv")
baseline_file <- file.path(base_dir, "KEN_baseline_detectability.parquet")
if (!all(file.exists(c(comparison_file, baseline_file)))) stop("Missing ecological-validation inputs.")

comparison <- fread(comparison_file)
baseline <- as.data.table(read_parquet(baseline_file))[has_mapped_crop == TRUE]

ecology <- baseline[, .(
  mapped_crop_pixels = .N,
  arid_share = mean(aridity_bin == "arid", na.rm = TRUE),
  semi_arid_share = mean(aridity_bin == "semi-arid", na.rm = TRUE),
  sub_humid_share = mean(aridity_bin == "sub-humid", na.rm = TRUE),
  humid_share = mean(aridity_bin == "humid", na.rm = TRUE),
  median_ndvi_amplitude = median(ndvi_amplitude_median, na.rm = TRUE),
  median_event_coverage = median(event_coverage, na.rm = TRUE),
  median_rainfall_signal = median(baseline_rainfall_signal, na.rm = TRUE),
  median_valley_strength = median(valley_strength, na.rm = TRUE),
  ndvi_pathway_share = mean(baseline_pathway == "ndvi_seasonal"),
  rainfall_proxy_share = mean(baseline_pathway == "rainfall_proxy"),
  not_identifiable_share = mean(baseline_pathway == "not_identifiable")
), by = .(admin1_name, season_id)]

share_columns <- c("arid_share", "semi_arid_share", "sub_humid_share", "humid_share")
share_names <- c("arid", "semi-arid", "sub-humid", "humid")
ecology[, dominant_aridity := share_names[max.col(.SD, ties.method = "first")], .SDcols = share_columns]

diagnostics <- ecology[comparison, on = .(admin1_name, season_id)]
diagnostics[, sufficient_validation_sample :=
  observed_records >= config$validation$ceepa$min_observed_records &
  remote_quality_pixels >= config$validation$ceepa$min_remote_pixels
]
diagnostics[, timing_class := fifelse(
  !sufficient_validation_sample, "insufficient_sample",
  fifelse(
    absolute_timing_difference_days <= 30, "aligned",
    fifelse(absolute_timing_difference_days < 45, "review", "discordant")
  )
)]
setorder(diagnostics, -absolute_timing_difference_days, admin1_name, season_id)
fwrite(diagnostics, file.path(report_dir, "KEN_ceepa_ecological_diagnostics.csv"), na = "")
write_parquet(diagnostics, file.path(report_dir, "KEN_ceepa_ecological_diagnostics.parquet"))

summary <- diagnostics[sufficient_validation_sample == TRUE, .(
  county_seasons = .N,
  median_absolute_difference = median(absolute_timing_difference_days, na.rm = TRUE),
  median_humid_share = median(humid_share, na.rm = TRUE),
  median_not_identifiable_share = median(not_identifiable_share, na.rm = TRUE)
), by = .(timing_class, dominant_aridity)]
setorder(summary, timing_class, dominant_aridity)
fwrite(summary, file.path(report_dir, "KEN_ceepa_ecological_summary.csv"))

correlations <- rbindlist(lapply(c(
  "humid_share", "median_ndvi_amplitude", "median_event_coverage",
  "median_rainfall_signal", "median_valley_strength", "not_identifiable_share"
), function(metric) {
  data.table(
    metric = metric,
    spearman_rho = cor(
      diagnostics[sufficient_validation_sample == TRUE][[metric]],
      diagnostics[sufficient_validation_sample == TRUE]$absolute_timing_difference_days,
      method = "spearman",
      use = "complete.obs"
    )
  )
}))
fwrite(correlations, file.path(report_dir, "KEN_ceepa_ecological_correlations.csv"))

plot_data <- diagnostics[sufficient_validation_sample == TRUE]
p <- ggplot(
  plot_data,
  aes(humid_share, absolute_timing_difference_days, colour = survey_season_name)
) +
  geom_hline(yintercept = c(30, 45), linetype = c("dashed", "solid"), colour = "grey55") +
  geom_point(aes(size = observed_records, shape = timing_class), alpha = 0.8) +
  geom_text(
    data = plot_data[timing_class == "discordant"],
    aes(label = admin1_name),
    size = 2.6, check_overlap = TRUE, vjust = -0.8, show.legend = FALSE
  ) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    title = "CEEPA–GLASS timing differences concentrate in humid counties",
    subtitle = "County-season medians; GLASS greenup and farmer-reported planting are related but non-equivalent events",
    x = "Humid share of mapped-crop pixels",
    y = "Absolute timing difference (days)",
    colour = "Survey season", size = "CEEPA records", shape = "Diagnostic"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(
  file.path(report_dir, "KEN_ceepa_humid_timing_diagnostics.png"),
  p, width = 10, height = 7, dpi = 180, bg = "white"
)

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    county_seasons = nrow(diagnostics),
    sufficiently_sampled = diagnostics[sufficient_validation_sample == TRUE, .N],
    discordant = diagnostics[timing_class == "discordant", .N],
    report_dir = report_dir
  )
)
cat("Wrote ecological validation diagnostics under", report_dir, "\n")
