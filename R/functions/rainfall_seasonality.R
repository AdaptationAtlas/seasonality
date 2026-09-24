# Rainfall seasonality metrics from 12 monthly totals.

rainfall_harmonic <- function(monthly_rain, harmonic = 1L) {
  stopifnot(length(monthly_rain) == 12L, harmonic %in% c(1L, 2L))
  rain <- as.numeric(monthly_rain)
  total <- sum(rain, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NA_real_)

  theta <- 2 * pi * ((seq_len(12) - 0.5) / 12)
  x <- sum(rain * cos(harmonic * theta), na.rm = TRUE)
  y <- sum(rain * sin(harmonic * theta), na.rm = TRUE)
  sqrt(x^2 + y^2) / total
}

rainfall_seasonality_metrics <- function(monthly_rain) {
  stopifnot(length(monthly_rain) == 12L)
  rain <- as.numeric(monthly_rain)
  valid <- is.finite(rain)
  total <- sum(rain[valid])

  if (!all(valid) || total <= 0) {
    return(list(
      months_observed = sum(valid),
      annual_rain = if (any(valid)) total else NA_real_,
      rainfall_si = NA_real_,
      rain_h1 = NA_real_,
      rain_h2 = NA_real_,
      wettest_month = NA_integer_
    ))
  }

  list(
    months_observed = 12L,
    annual_rain = total,
    rainfall_si = sum(abs(rain - total / 12)) / total,
    rain_h1 = rainfall_harmonic(rain, 1L),
    rain_h2 = rainfall_harmonic(rain, 2L),
    wettest_month = which.max(rain)
  )
}
