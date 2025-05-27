#' Download AgERA5 v2.0 Climate Data via CDS API
#'
#' Downloads daily AgERA5 reanalysis climate data using the CDS API for specified years,
#' months, and a defined bounding box. Supports batch downloads using `wf_request_batch`.
#'
#' @param variable Character. The climate variable to download (e.g., "solar_radiation_flux").
#' @param stat Character or NULL. The statistic to request (e.g., "24_hour_mean").
#' @param time Character or NULL. Optional time string for certain variables (e.g., "12_00").
#' @param years Character vector of years (e.g., "1995":"2014").
#' @param months Character vector of 2-digit months (e.g., "01" to "12").
#' @param datadir Output directory where files are saved.
#' @param version Dataset version. Default is "2_0".
#' @param workers Integer. Number of concurrent download workers.
#' @param bbox Numeric vector of bounding box [xmin, xmax, ymin, ymax]. Default is Africa.
#' @param format File format for download. Default is "zip".
#' @param timeout Request timeout in seconds.
#' @param retry Number of times to retry failed requests.
#'
#' @return NULL. Writes files to disk.
#'
#' @import data.table
#' @importFrom lubridate days_in_month
#' @importFrom ecmwfr wf_request_batch
#' @importFrom utils file.exists
#' @export
#'
#' @examples
#' \dontrun{
#' # Set UID and key beforehand using ecmwfr::wf_set_key()
#'
#' years  <- as.character(1995:2014)
#' months <- sprintf("%02d", 1:12)
#' datadir <- "~/agera5_downloads"
#'
#' vars_to_download <- data.frame(
#'   variable = c("solar_radiation_flux", "2m_relative_humidity", "10m_wind_speed"),
#'   statistics = c(NA, NA, "24_hour_mean"),
#'   time = c(NA, "12_00", NA),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (i in seq_len(nrow(vars_to_download))) {
#'   getAgERA5_v2(
#'     variable = vars_to_download$variable[i],
#'     stat     = vars_to_download$statistics[i],
#'     time     = vars_to_download$time[i],
#'     years    = years,
#'     months   = months,
#'     datadir  = datadir,
#'     workers  = 4
#'   )
#' }
#' }
getAgERA5_v2 <- function(variable,
                         stat,
                         time,
                         years,
                         months,
                         datadir,
                         version = "2_0",
                         workers = 2,
                         bbox    = c(-26, 64, -35, 38),
                         format  = "zip",
                         time_out = 3600,
                         retry   = 5) {

  # install / load deps
  pkgs <- c("data.table","ecmwfr")
  for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  lapply(pkgs, library, character.only=TRUE)


  start_time <- Sys.time()
  inputs <- data.table::data.table(expand.grid(year = years, month = months, stringsAsFactors = FALSE))

  inputs[, path := file.path(datadir, variable,
                             paste0(paste(c(variable, na.omit(stat), na.omit(time), year, month),
                                          collapse = "-"), ".", format)), by = .I
  ][, file := basename(path)
  ][, exists := file.exists(path)]

  inputs <- inputs[exists == FALSE]
  files_total <- nrow(inputs)

  if (files_total > 0) {
    dir_name <- dirname(inputs$path[1])
    if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)

    message("→ Downloading ", variable, " [", stat, "][", time, "] for years: ",
            min(years), "-", max(years), " & months: ", min(months), "-", max(months))

    area_str <- c(bbox[4], bbox[1], bbox[3], bbox[2]) # N/W/S/E

    request_list <- lapply(1:nrow(inputs), function(i) {
      request <- list(
        dataset_short_name = "sis-agrometeorological-indicators",
        version            = version,
        variable           = variable,
        statistic          = stat,
        year               = inputs$year[i],
        month              = inputs$month[i],
        day                = sprintf("%02d", seq_len(lubridate::days_in_month(
         as.Date(paste0(inputs$year[i], "-", inputs$month[i], "-01"))))),
        area               = area_str,
        time               = time,
        format             = format,
        target             = inputs$file[i]
      )
      request <- Filter(Negate(is.null), request)
      request[!sapply(request, function(x) all(is.na(x)))]
    })

    ecmwfr::wf_request_batch(
      user         = UID,
      request_list = request_list,
      path         = dir_name,
      workers      = workers,
      retry        = retry,
      time_out      = timeout
    )

    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

    log_line <- paste0(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%OSZ"), "\t",
      paste0(min(years), "-", max(years)), "\t",
      paste(bbox, collapse = " "), "\t",
      1, "\t",
      workers, "\t",
      files_total, "\t",
      sum(file.exists(inputs$path)), "\t",
      files_total - sum(file.exists(inputs$path)), "\t",
      datadir, "\t",
      "sis-agrometeorological-indicators", "\t",
      "https://cds.climate.copernicus.eu/api/v2", "\t",
      timeout, "\t",
      round(duration, 2), "\t",
      Sys.getenv("USER"), "\t",
      Sys.info()[["nodename"]]
    )
    cat(log_line, "\n")

  } else {
    message("✔ Already have all files.")
  }

  invisible(NULL)
}
