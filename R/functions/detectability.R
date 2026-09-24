# Metrics and explicit pathway assignment for seasonal-signal detectability.

circular_concentration <- function(doy, period = 365) {
  doy <- doy[is.finite(doy)]
  if (!length(doy)) return(NA_real_)
  angle <- 2 * pi * doy / period
  sqrt(mean(cos(angle))^2 + mean(sin(angle))^2)
}

annual_ndvi_metrics <- function(ndvi, dates, min_valid_fraction = 0.75) {
  stopifnot(length(ndvi) == length(dates))
  dates <- as.Date(dates)
  years <- as.integer(format(dates, "%Y"))

  pieces <- lapply(split(seq_along(dates), years), function(i) {
    x <- as.numeric(ndvi[i])
    valid <- is.finite(x)
    valid_fraction <- mean(valid)

    if (!any(valid)) {
      return(data.frame(
        year = years[i[1]], valid_fraction = 0,
        ndvi_q10 = NA_real_, ndvi_q50 = NA_real_, ndvi_q90 = NA_real_,
        ndvi_amplitude = NA_real_, ndvi_mean = NA_real_, ndvi_sd = NA_real_,
        sufficient_data = FALSE
      ))
    }

    q <- stats::quantile(x[valid], c(0.1, 0.5, 0.9), na.rm = TRUE, names = FALSE)
    data.frame(
      year = years[i[1]],
      valid_fraction = valid_fraction,
      ndvi_q10 = q[1],
      ndvi_q50 = q[2],
      ndvi_q90 = q[3],
      ndvi_amplitude = q[3] - q[1],
      ndvi_mean = mean(x[valid]),
      ndvi_sd = stats::sd(x[valid]),
      sufficient_data = valid_fraction >= min_valid_fraction
    )
  })

  do.call(rbind, pieces)
}

classify_signal_path <- function(
    sufficient_data,
    ndvi_amplitude,
    timing_concentration,
    event_coverage,
    rainfall_seasonality = NA_real_,
    expected_seasons = NA_integer_,
    detected_seasons = NA_integer_,
    wet_anomaly = NA_real_,
    thresholds) {
  required <- c(
    "ndvi_amplitude_min", "timing_concentration_min", "event_coverage_min",
    "rainfall_seasonality_min", "wet_anomaly_min"
  )
  missing_thresholds <- setdiff(required, names(thresholds))
  if (length(missing_thresholds)) {
    stop("Missing thresholds: ", paste(missing_thresholds, collapse = ", "))
  }

  n <- max(
    length(sufficient_data), length(ndvi_amplitude),
    length(timing_concentration), length(event_coverage),
    length(rainfall_seasonality), length(expected_seasons),
    length(detected_seasons), length(wet_anomaly)
  )
  recycle <- function(x) rep_len(x, n)
  sufficient_data <- recycle(sufficient_data)
  ndvi_amplitude <- recycle(ndvi_amplitude)
  timing_concentration <- recycle(timing_concentration)
  event_coverage <- recycle(event_coverage)
  rainfall_seasonality <- recycle(rainfall_seasonality)
  expected_seasons <- recycle(expected_seasons)
  detected_seasons <- recycle(detected_seasons)
  wet_anomaly <- recycle(wet_anomaly)

  out <- rep("not_identifiable", n)
  out[!sufficient_data | is.na(sufficient_data)] <- "insufficient_data"

  merged <- sufficient_data &
    !is.na(expected_seasons) & expected_seasons >= 2L &
    !is.na(detected_seasons) & detected_seasons < expected_seasons &
    !is.na(wet_anomaly) & wet_anomaly >= thresholds$wet_anomaly_min
  out[merged] <- "wet_merged"

  ndvi_ok <- sufficient_data & !merged &
    !is.na(ndvi_amplitude) & ndvi_amplitude >= thresholds$ndvi_amplitude_min &
    !is.na(timing_concentration) & timing_concentration >= thresholds$timing_concentration_min &
    !is.na(event_coverage) & event_coverage >= thresholds$event_coverage_min
  out[ndvi_ok] <- "ndvi_seasonal"

  rain_ok <- sufficient_data & !merged & !ndvi_ok &
    !is.na(rainfall_seasonality) &
    rainfall_seasonality >= thresholds$rainfall_seasonality_min
  out[rain_ok] <- "rainfall_proxy"

  factor(
    out,
    levels = c(
      "ndvi_seasonal", "rainfall_proxy", "wet_merged",
      "not_identifiable", "insufficient_data"
    )
  )
}

classify_baseline_signal <- function(
    season_id,
    aridity_bin,
    years_with_ndvi,
    ndvi_amplitude,
    event_coverage,
    timing_concentration,
    rainfall_signal,
    valley_strength,
    bimodal_year_coverage,
    thresholds) {
  required <- c(
    "min_years_with_ndvi", "ndvi_amplitude_min", "humid_ndvi_amplitude_min",
    "event_coverage_min", "timing_concentration_min", "rainfall_signal_min",
    "rainfall_signal_floor", "valley_strength_min", "bimodal_year_coverage_min"
  )
  missing_thresholds <- setdiff(required, names(thresholds))
  if (length(missing_thresholds)) {
    stop("Missing thresholds: ", paste(missing_thresholds, collapse = ", "))
  }

  n <- max(
    length(season_id), length(aridity_bin), length(years_with_ndvi),
    length(ndvi_amplitude), length(event_coverage), length(timing_concentration),
    length(rainfall_signal), length(valley_strength), length(bimodal_year_coverage)
  )
  recycle <- function(x) rep_len(x, n)
  season_id <- recycle(season_id)
  aridity_bin <- as.character(recycle(aridity_bin))
  years_with_ndvi <- recycle(years_with_ndvi)
  ndvi_amplitude <- recycle(ndvi_amplitude)
  event_coverage <- recycle(event_coverage)
  timing_concentration <- recycle(timing_concentration)
  rainfall_signal <- recycle(rainfall_signal)
  valley_strength <- recycle(valley_strength)
  bimodal_year_coverage <- recycle(bimodal_year_coverage)

  enough_years <- !is.na(years_with_ndvi) &
    years_with_ndvi >= thresholds$min_years_with_ndvi
  amplitude_min <- ifelse(
    aridity_bin == "humid",
    thresholds$humid_ndvi_amplitude_min,
    thresholds$ndvi_amplitude_min
  )
  ndvi_evidence <- enough_years &
    !is.na(ndvi_amplitude) & ndvi_amplitude >= amplitude_min &
    !is.na(event_coverage) & event_coverage >= thresholds$event_coverage_min &
    !is.na(timing_concentration) & timing_concentration >= thresholds$timing_concentration_min

  rain_evidence <- !is.na(rainfall_signal) &
    rainfall_signal >= thresholds$rainfall_signal_min
  rain_bimodal <- rain_evidence & !is.na(valley_strength) &
    valley_strength >= thresholds$valley_strength_min
  ndvi_bimodal <- ndvi_evidence &
    !is.na(bimodal_year_coverage) &
    bimodal_year_coverage >= thresholds$bimodal_year_coverage_min &
    !is.na(rainfall_signal) & rainfall_signal >= thresholds$rainfall_signal_floor
  bimodal_supported <- rain_bimodal | ndvi_bimodal
  season_supported <- season_id == 1L | (season_id == 2L & bimodal_supported)

  pathway <- rep("not_identifiable", n)
  pathway[!enough_years] <- "insufficient_data"
  pathway[enough_years & season_supported & rain_evidence] <- "rainfall_proxy"
  pathway[enough_years & season_supported & ndvi_evidence] <- "ndvi_seasonal"

  data.frame(
    baseline_pathway = factor(
      pathway,
      levels = c("ndvi_seasonal", "rainfall_proxy", "not_identifiable", "insufficient_data")
    ),
    ndvi_evidence = ndvi_evidence,
    rainfall_evidence = rain_evidence,
    bimodal_supported = bimodal_supported,
    season_supported = season_supported
  )
}

classify_annual_pathway <- function(
    baseline_pathway,
    season_id,
    quality_event,
    bimodal_supported,
    detected_seasons,
    wet_anomaly,
    wet_anomaly_min) {
  baseline_pathway <- as.character(baseline_pathway)
  out <- baseline_pathway
  out[baseline_pathway == "ndvi_seasonal" & quality_event] <- "ndvi_event"
  out[baseline_pathway == "ndvi_seasonal" & !quality_event] <- "missing_ndvi_event"

  wet_merged <- season_id == 2L & bimodal_supported & !quality_event &
    !is.na(detected_seasons) & detected_seasons < 2L &
    !is.na(wet_anomaly) & wet_anomaly >= wet_anomaly_min
  out[wet_merged] <- "wet_merged"

  factor(
    out,
    levels = c(
      "ndvi_event", "rainfall_proxy", "wet_merged", "missing_ndvi_event",
      "not_identifiable", "insufficient_data"
    )
  )
}
