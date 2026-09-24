#' @title ERA Season Assignment
#' @description
#' Functions to derive seasonal structure from SOS/EOS using circular clustering.

NULL

#' Convert date to DOY (365-day)
#' @export
doy365 <- function(x) {
  d <- lubridate::yday(x)
  leap <- lubridate::leap_year(x)
  d[leap & d > 59] <- d[leap & d > 59] - 1L
  d
}

#' Midpoint of circular interval
#' @export
circular_mid_doy <- function(sos, eos) {

  s <- doy365(sos)
  e <- doy365(eos)

  dur <- e - s
  dur[dur < 0] <- dur[dur < 0] + 365L

  mid <- s + dur / 2

  ((mid - 1) %% 365) + 1
}

#' Circular mean (angle)
#' @export
circ_mean_angle <- function(theta) {
  atan2(mean(sin(theta)), mean(cos(theta))) %% (2 * pi)
}

#' Fit circular clusters (season detection)
#' @export
fit_circular_seasons <- function(doy,
                                 max_seasons = 3L,
                                 min_prop = 0.05,
                                 nstart = 50L) {

  doy <- as.numeric(doy)
  ok <- !is.na(doy)

  out <- list(
    k = NA_integer_,
    cluster = rep(NA_integer_, length(doy))
  )

  if (sum(ok) < 3L) {
    out$k <- 1L
    out$cluster[ok] <- 1L
    return(out)
  }

  theta <- 2 * pi * (doy[ok] - 1) / 365
  xy <- cbind(cos(theta), sin(theta))

  n <- nrow(xy)
  p <- ncol(xy)
  k_max <- min(max_seasons, n - 1L)

  best_bic <- Inf
  best_km <- NULL
  best_k <- 1L

  for (k in 1L:k_max) {

    km <- kmeans(xy, centers = k, nstart = nstart, iter.max = 100)

    props <- tabulate(km$cluster, nbins = k) / n
    if (any(props < min_prop)) next

    bic <- n * log(km$tot.withinss / n) + log(n) * k * p

    if (bic < best_bic) {
      best_bic <- bic
      best_km <- km
      best_k <- k
    }
  }

  if (is.null(best_km)) {
    best_km <- kmeans(xy, centers = 1L, nstart = 1L)
    best_k <- 1L
  }

  cl <- best_km$cluster

  mu <- tapply(theta, cl, circ_mean_angle)
  ord <- order(mu)

  relab <- setNames(seq_along(ord), names(mu)[ord])
  cl2 <- unname(relab[as.character(cl)])

  out$k <- best_k
  out$cluster[ok] <- cl2

  out
}

#' Assign seasons within admin unit
#' @export
assign_admin_seasons <- function(sos, eos, max_seasons = 3L) {

  mid_doy <- circular_mid_doy(sos, eos)

  fit <- fit_circular_seasons(
    doy = mid_doy,
    max_seasons = max_seasons
  )

  list(
    mid_doy = mid_doy,
    season = fit$cluster,
    n_seasons = rep(fit$k, length(mid_doy))
  )
}
