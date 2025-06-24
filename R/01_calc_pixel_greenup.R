#' Analyse one NDVI pixel time-series with {phenofit}
#'
#' Fits seasonal curves to an NDVI time-series from **one raster pixel**
#'    and returns the dates of phenological events (SOS, EOS, UD …),
#'    together with the goodness-of-fit (GOF) of the *best* curve
#'    chosen per season.
#'
#' @param ndvi_ts Numeric vector – raw NDVI values (may contain NA & water code)
#' @param dates   Date vector   – same length as ndvi_ts
#' @param nptperyear Integer – points per year (46 for 8-day, 23 for 16-day, …)
#' @param w_low    Numeric – weight to assign to missing observations (default 0.1)
#' @param min_frac    Numeric – minimum fraction of valid observations (default 0.20)
#' @param min_years   Integer – minimum number of distinct calendar years (default 3)
#'
#' @details
#' * Seasons are first segmented with `phenofit::season_mov()` using a
#'   Whittaker smoother (`smooth_wWHIT` + adaptive weights `wTSM`).
#' * Each season is then fitted with every model in `methods_to_fit`.
#'   The model giving the **lowest RMSE** is taken as *best* for that
#'   season (ties broken by NSE).
#' * Output combines `get_pheno()` dates **and** the GOF statistics
#'   (`R2`, `NSE`, `RMSE`, `R`) for the winning model.
#'   get_pheno uses method ="AG" with TRS = c(0.2, 0.5)
#'
#' @return  A `data.table` of phenology dates (one row per detected season) **or**
#'          `NULL` if the time series fails QC or no season is found
#'
#' @import data.table
#' @import phenofit
#'
#'  @export
calc_pixel_greenup <- function(ndvi_ts,
                             dates,
                             nptperyear  = 46,
                             w_low       = 0.1,
                             methods_to_fit= c("Beck", "Elmore","Zhang"),
                             min_frac    = 0.20,
                             min_years   = 3) {

  ##───────────────────────────────────────────────
  ## 0. QC & data-sufficiency filter             ░
  ##───────────────────────────────────────────────
  stopifnot(length(ndvi_ts) == length(dates))

  ts_tbl <- data.frame(t = dates, y = as.numeric(ndvi_ts))
  ts_tbl <- ts_tbl[!is.na(ts_tbl$y), ]                 # drop NA

  # adaptive test: require ≥ min_frac valid obs AND cover ≥ min_years
  obs_ratio     <- nrow(ts_tbl) / length(ndvi_ts)
  unique_years  <- length(unique(format(ts_tbl$t, "%Y")))
  if (obs_ratio < min_frac || unique_years < min_years) return(NULL)

  ##───────────────────────────────────────────────
  ## 1. Build weights – missing → low weight      ░
  ##───────────────────────────────────────────────
  w <- rep(1, length(ndvi_ts))
  w[is.na(ndvi_ts)] <- w_low

  ##───────────────────────────────────────────────
  ## 2. Build INPUT object (no ‘south’ flag)      ░
  ##───────────────────────────────────────────────
  INPUT <- phenofit::check_input(
    t            = dates,
    y            = as.numeric(ndvi_ts),
    w            = w,
    nptperyear   = nptperyear,
    maxgap       = floor(nptperyear / 4),   # tolerate gaps ≤ ¼ yr
    wmin         = w_low
  )

  ##───────────────────────────────────────────────
  ## 3. Season segmentation & curve fitting       ░
  ##───────────────────────────────────────────────
  brks <-suppressWarnings(phenofit::season_mov(
    INPUT,
    list(
      FUN             = "smooth_wWHIT",      # Whittaker smoother
      wFUN            = phenofit::wTSM,      # adaptive weights
      maxExtendMonth  = 3,                   # allow ±3-mo extension
      wmin            = w_low,
      r_min           = 0.1                  # min relative amplitude
    )
  ))

 if(!is.null(brks)){
    fit <- curvefits(
      INPUT, brks,
      list(methods = methods_to_fit,
           wFUN    = wTSM,
           iters   = 2)
    )

    gof <- data.table(phenofit::get_GOF(fit))
    best <- gof[ , .SD[which.min(RMSE)], by = flag][,n_sim:=NULL]

    # Extract phenology metrics
    pheno <- phenofit::get_pheno(fit, method = "AG", TRS = c(0.2, 0.5))$date$AG

    pheno <- merge(pheno,                 # phenology dates
                       best,
                       by = "flag",
                       all.x = TRUE)

    # If no season detected → NULL
    if (is.null(pheno) || nrow(pheno) == 0) {
      return(NULL)
    }else{
      return(pheno)
    }
 }else{
   warning("Can't find any growing seasons!")
   return(NULL)
 }
}
