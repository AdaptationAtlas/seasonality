#' Download AgERA5 v2.0 Climate Data via CDS API
#'
#' Downloads daily AgERA5 reanalysis climate data using the CDS API for specified years,
#' months, and a defined bounding box. Supports batch downloads with a parallel `furrr` map
#' and retry logic.
#'
#' **Pre‐requisite:** you must first configure your Copernicus CDS API key.
#' For example:
#' ```r
#' # once-only, in your R startup:
#' ecmwfr::wf_set_key(
#'   user = "<YOUR_CDS_EMAIL>@domain.org",
#'   key  = "<YOUR_CDSAPI_KEY>"
#' )
#' ```
#' You can obtain your CDS API key by logging in at
#' <https://cds.climate.copernicus.eu/api-how-to>.
#'
#' The number of concurrent workers should reflect your system resources and the Copernicus CDS API’s fair-use policy.
#’ Up to 4–6 workers is typically safe. Using more may result in throttling or queued requests by the CDS server.
#'
#’ For more details on available parameters and structure of the requests, refer to the documentation for
#’ ecmwfr::wf_request() and ecmwfr::wf_request_batch() functions in the ecmwfr package.
#'
#' @param variable Character. The climate variable to download (e.g., "solar_radiation_flux").
#' @param stat Character or NULL. Statistic to request (e.g., "24_hour_mean").
#' @param time Character or NULL. Optional time string for certain variables (e.g., "12_00").
#' @param years Character vector of years (e.g., "1995" to "2014").
#' @param months Character vector of 2-digit months (e.g., "01" to "12").
#' @param datadir Output directory where files are saved.
#' @param version Character dataset version. Default is "2_0".
#' @param workers Integer. Number of concurrent futures/workers.
#' @param bbox Numeric vector [xmin, xmax, ymin, ymax]. Default is Africa.
#' @param format Character file format. Default is "zip".
#' @param timeout Numeric seconds before timing out. Default 3600.
#' @param retry Integer number of download retries. Default 5.
#' @import data.table
#' @importFrom lubridate days_in_month
#' @importFrom furrr future_map
#' @importFrom future plan multisession
#' @importFrom ecmwfr wf_set_key
#' @export
getAgERA5_v2 <- function(variable,
                         stat = NULL,
                         time = NULL,
                         years,
                         months,
                         datadir,
                         version = "2_0",
                         workers = 4,
                         user = "ecmwfr",
                         bbox    = c(-26, 64, -35, 38),
                         format  = "zip",
                         timeout = 3600,
                         parallel_retries = 3,
                         retry   = 5) {


    start_time <- Sys.time()

  # dependencies
  # Vector of required packages
  pkgs <- c("data.table", "lubridate", "furrr", "future", "ecmwfr", "progressr")

  # Install missing ones
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    library(p, character.only = TRUE)
  }

  dt <- data.table::as.data.table(expand.grid(year=years, month=months, stringsAsFactors=FALSE))
  dt[, path := file.path(datadir, variable,
                     paste0(paste(c(variable, na.omit(stat), na.omit(time), year, month), collapse = "-"), ".", format)),by=.I
     ][, file := basename(path),
       ][,ok := file.exists(path)]

  pending <- dt[ok == FALSE]
  total <- nrow(dt)
  already_complete<-total-nrow(pending)
  to_do <- nrow(pending)
  if (to_do == 0) {
    message("✔ Already have all ", total, " files for ", variable)
    return(invisible(NULL))
  }

  # create directory
  dir.create(dirname(pending$path[1]), recursive = TRUE, showWarnings = FALSE)

  # area ordering N/W/S/E
  area_vec <- c(bbox[4], bbox[1], bbox[3], bbox[2])

  # prepare requests
  requests <- lapply(1:nrow(pending), function(i) {
    req<-list(
      dataset_short_name = "sis-agrometeorological-indicators",
      version            = version,
      variable           = variable,
      statistic          = stat,
      year               = pending$year[i],
      month              = pending$month[i],
      day                = sprintf("%02d", seq_len(lubridate::days_in_month(
      as.Date(paste0(pending$year[i], "-", pending$month[i], "-01"))))),
      area               = area_vec,
      time               = time,
      format             = format,
      target             = pending$file[i]
    )
    req <- Filter(Negate(is.null), req)
    req[!sapply(req, function(x) all(is.na(x)))]
  })

  # parallel download with retries
  future::plan(future::multisession, workers = workers)

  handlers(global = TRUE); handlers("progress")

  results <- progressr::with_progress({
    p <- progressr::progressor(length(requests))

  results <- furrr::future_map(seq_len(length(requests)), function(i) {
    attempt <- 1L
    success <- FALSE
    while(attempt <= parallel_retries && !success) {
      try({
        ecmwfr::wf_request(
          request      = requests[[i]],
          transfer     = TRUE,
          path         = dirname(pending$path[1]),
          time_out     = timeout,
          retry        = retry,
          user         = user,
          verbose      = FALSE
        )
        success <- TRUE
      }, silent = TRUE)
      attempt <- attempt + 1L
    }
    p(sprintf("→ %s [%d/%d] %s", requests[[i]]$variable, i, length(requests),
              if (success) "✓" else "✗"))

    data.table::data.table(
      seq = i,
      ok  = success
    )
  }, .options = furrr::furrr_options(seed = TRUE))
  })

  future::plan(sequential)
  res_dt <- data.table::rbindlist(results)


  # build and write run‐level log
  end_time   <- Sys.time()
  duration_s <- as.numeric(difftime(end_time, start_time, units = "secs"))

  n_ok <- sum(res_dt$ok) + already_complete
  files_failed <- total - n_ok

  log_dt <- data.table::data.table(
    timestamp    = format(end_time, "%Y-%m-%dT%H:%M:%SZ"),
    variable     = variable,
    statistic    = ifelse(is.null(stat), NA_character_, stat),
    time         = ifelse(is.null(time),  NA_character_, time),
    years        = paste(min(years),  max(years),  sep = "-"),
    months       = paste(min(months), max(months), sep = "-"),
    bbox         = paste(bbox, collapse = " "),
    workers      = workers,
    files_total  = total,
    files_ok     = n_ok,
    files_failed = total - n_ok,
    out_dir      = normalizePath(datadir),
    dataset_id   = "sis-agrometeorological-indicators",
    base_url     = "https://cds.climate.copernicus.eu/api/v2",
    timeout      = timeout,
    duration_s   = duration_s,
    user         = Sys.getenv("USER"),
    host         = Sys.info()[["nodename"]]
  )

  logf <- file.path(datadir, "agera5_v2_log.csv")
  if (file.exists(logf)) {
    data.table::fwrite(log_dt, logf, append = TRUE)
  } else {
    data.table::fwrite(log_dt, logf)
  }

  # success message at end
  message(sprintf("✔ sis-agrometeorological-indicators: %d/%d OK → out_dir @ %s",
                  n_ok, total, normalizePath(datadir)))
  invisible(res_dt)
}
