#' @title ERA Circular Utilities
#' @description
#' Utility functions for working with circular variables (day-of-year) in ERA phenology workflows.
#' Includes circular quantiles, distances, and interval logic.

NULL

#' Circular quantiles for day-of-year data
#'
#' Computes quantiles for circular variables (e.g. DOY), avoiding edge artefacts.
#'
#' @param doy Numeric vector (1–365)
#' @param probs Quantiles
#' @param n_days Cycle length
#'
#' @return Named numeric vector
#' @export
quantile_circular_safe <- function(doy, probs = c(0.1, 0.5, 0.9), n_days = 365) {

  doy <- doy[!is.na(doy)]

  if (!length(doy)) {
    out <- rep(NA_real_, length(probs))
    names(out) <- paste0(probs * 100, "%")
    return(out)
  }

  ang <- 2 * pi * (doy / n_days)

  mu <- atan2(mean(sin(ang)), mean(cos(ang)))
  if (mu < 0) mu <- mu + 2 * pi

  centre <- mu * n_days / (2 * pi)

  x <- doy
  x[x < centre - n_days / 2] <- x[x < centre - n_days / 2] + n_days
  x[x > centre + n_days / 2] <- x[x > centre + n_days / 2] - n_days

  qs <- quantile(x, probs = probs, na.rm = TRUE)

  ((qs - 1) %% n_days) + 1
}

#' Season length on circular calendar
#' @export
season_length_doy <- function(sos, eos, n_days = 365) {
  (eos - sos + n_days) %% n_days
}

#' Circular distance (shortest arc)
#' @export
circ_dist <- function(x, y, n_days = 365) {
  d <- abs(x - y)
  pmin(d, n_days - d)
}

#' Signed shortest circular difference (x minus y)
#' @export
signed_circ_diff <- function(x, y, n_days = 365) {
  ((x - y + n_days / 2) %% n_days) - n_days / 2
}

#' Forward circular distance
#' @export
forward_circ_dist <- function(start, end, n_days = 365) {
  start <- rep_len(start, max(length(start), length(end)))
  end   <- rep_len(end,   max(length(start), length(end)))
  (end - start + n_days) %% n_days
}

#' Membership in circular window
#' @export
in_circular_window <- function(x, start, end, n_days = 365) {

  n <- max(length(x), length(start), length(end))

  x     <- rep_len(x, n)
  start <- rep_len(start, n)
  end   <- rep_len(end, n)

  out <- rep(NA, n)

  ok <- !is.na(x) & !is.na(start) & !is.na(end)

  nowrap <- ok & (start <= end)
  wrap   <- ok & (start > end)

  out[nowrap] <- x[nowrap] >= start[nowrap] & x[nowrap] <= end[nowrap]
  out[wrap]   <- x[wrap] >= start[wrap] | x[wrap] <= end[wrap]

  out
}
