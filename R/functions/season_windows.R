# Stable candidate season windows from monthly rainfall climatology.

circular_month_distance <- function(a, b) {
  distance <- abs(a - b)
  pmin(distance, 12 - distance)
}

in_circular_month_window <- function(month, start, end) {
  n <- max(length(month), length(start), length(end))
  month <- rep_len(month, n)
  start <- rep_len(start, n)
  end <- rep_len(end, n)
  out <- rep(NA, n)
  ok <- !is.na(month) & !is.na(start) & !is.na(end)
  out[ok & start <= end] <- month[ok & start <= end] >= start[ok & start <= end] &
    month[ok & start <= end] <= end[ok & start <= end]
  out[ok & start > end] <- month[ok & start > end] >= start[ok & start > end] |
    month[ok & start > end] <= end[ok & start > end]
  out
}

circular_month_sequence <- function(start, end, include_start = TRUE, include_end = TRUE) {
  sequence <- ((start - 1L + 0:11) %% 12L) + 1L
  end_position <- match(end, sequence)
  sequence <- sequence[seq_len(end_position)]
  if (!include_start) sequence <- sequence[-1L]
  if (!include_end && length(sequence)) sequence <- sequence[-length(sequence)]
  sequence
}

smooth_monthly_circular <- function(monthly_rain) {
  stopifnot(length(monthly_rain) == 12L)
  x <- as.numeric(monthly_rain)
  previous <- x[c(12L, 1:11)]
  following <- x[c(2:12, 1L)]
  (previous + 2 * x + following) / 4
}

find_rainfall_peaks <- function(monthly_rain, max_peaks = 2L, min_separation = 3L) {
  stopifnot(length(monthly_rain) == 12L, max_peaks >= 1L, min_separation >= 1L)
  smoothed <- smooth_monthly_circular(monthly_rain)
  if (!any(is.finite(smoothed)) || sum(smoothed, na.rm = TRUE) <= 0) {
    return(integer())
  }

  previous <- smoothed[c(12L, 1:11)]
  following <- smoothed[c(2:12, 1L)]
  candidates <- which(
    is.finite(smoothed) & smoothed >= previous & smoothed > following
  )
  if (!length(candidates)) candidates <- which.max(smoothed)
  candidates <- candidates[order(smoothed[candidates], decreasing = TRUE)]

  selected <- integer()
  for (candidate in candidates) {
    separated <- !length(selected) ||
      all(circular_month_distance(candidate, selected) >= min_separation)
    if (separated) selected <- c(selected, candidate)
    if (length(selected) >= max_peaks) break
  }
  sort(selected)
}

derive_two_season_windows <- function(monthly_rain, min_peak_separation = 3L) {
  stopifnot(length(monthly_rain) == 12L)
  smoothed <- smooth_monthly_circular(monthly_rain)
  peaks <- find_rainfall_peaks(
    monthly_rain,
    max_peaks = 2L,
    min_separation = min_peak_separation
  )

  if (length(peaks) < 2L) {
    return(data.frame(
      season_id = integer(), peak_month = integer(),
      window_start_month = integer(), window_end_month = integer(),
      peak_rain = numeric(), peak_share = numeric(),
      valley_strength = numeric()
    ))
  }

  peak_1 <- peaks[1]
  peak_2 <- peaks[2]
  arc_1_to_2 <- circular_month_sequence(
    peak_1, peak_2, include_start = FALSE, include_end = FALSE
  )
  arc_2_to_1 <- circular_month_sequence(
    peak_2, peak_1, include_start = FALSE, include_end = FALSE
  )
  trough_1_to_2 <- arc_1_to_2[which.min(smoothed[arc_1_to_2])]
  trough_2_to_1 <- arc_2_to_1[which.min(smoothed[arc_2_to_1])]

  peak_values <- smoothed[peaks]
  trough_values <- smoothed[c(trough_1_to_2, trough_2_to_1)]
  valley_strength <- 1 - max(trough_values) / min(peak_values)
  total <- sum(smoothed, na.rm = TRUE)

  data.frame(
    season_id = 1:2,
    peak_month = peaks,
    window_start_month = c(
      trough_2_to_1 %% 12L + 1L,
      trough_1_to_2 %% 12L + 1L
    ),
    window_end_month = c(trough_1_to_2, trough_2_to_1),
    peak_rain = peak_values,
    peak_share = peak_values / total,
    valley_strength = rep(valley_strength, 2L)
  )
}
