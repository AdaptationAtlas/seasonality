library(data.table)
library(lubridate)

# Collapse leap years onto a 365-day circle
doy365 <- function(x) {
  d <- yday(x)
  leap <- leap_year(x)
  # shift days after Feb 28 back by 1 in leap years
  d[leap & d > 59] <- d[leap & d > 59] - 1L
  d
}

# Circular midpoint between sos and eos on a 365-day year
# Works even when eos is in the next calendar year relative to sos
circular_mid_doy <- function(sos, eos) {
  s <- doy365(sos)
  e <- doy365(eos)

  dur <- e - s
  dur[dur < 0] <- dur[dur < 0] + 365L

  mid <- s + dur / 2
  ((mid - 1) %% 365) + 1
}

# Circular mean angle for a cluster
circ_mean_angle <- function(theta) {
  atan2(mean(sin(theta)), mean(cos(theta))) %% (2 * pi)
}

# Fit circular clusters using k-means on the unit circle
# and choose k with a simple BIC-like criterion
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

    # avoid silly tiny clusters
    if (any(props < min_prop)) next

    # BIC-like penalty
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

  # relabel clusters in circular calendar order
  cl <- best_km$cluster
  mu <- tapply(theta, cl, circ_mean_angle)
  ord <- order(mu)
  relab <- setNames(seq_along(ord), names(mu)[ord])
  cl2 <- unname(relab[as.character(cl)])

  out$k <- best_k
  out$cluster[ok] <- cl2
  out
}

# Wrapper for one admin unit
assign_admin_seasons <- function(sos, eos, max_seasons = 3L) {
  # choose one timing variable for clustering
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
