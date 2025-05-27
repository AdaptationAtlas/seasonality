# CHC DATA MODULE ####

#’ @importFrom future.apply future_lapply
#’ @importFrom progressr with_progress progressor handlers
#’ @importFrom data.table rbindlist data.table fread fwrite
#’ @importFrom fs dir_create file_exists file_info path path_file
#’ @importFrom glue glue
#’ @importFrom arrow write_parquet
#’
#’ @title Core downloader for daily CHC GeoTIFFs → COGs
#’
#’ @description
#’ Downloads daily CHC raster (e.g. CHIRPS-v3, CHIRTS-ERA5 tmax/tmin/vpd),
#’ optionally crops to a bounding box, rounds values, converts to
#’ an internal-tiled ZSTD-compressed COG, and logs progress.
#’
#’ @param start_year Integer. First year to download.
#’ @param end_year   Integer. Last year to download.
#’ @param out_dir    Character. Folder to write COGs (subdirs per year).
#’ @param varname    Character. File prefix, e.g. `"chirps-v3.0"` or `"CHIRTS-ERA5.daily_Tmax"`.
#’ @param base_url   Character. CHC URL root *without* trailing slash.
#’ @param dataset_id Character. Short ID for logs/index, e.g. `"chirps_v3"`.
#’ @param bbox       Optional numeric(4): c(xmin, ymin, xmax, ymax) for `-projwin`.
#’ @param round_digits Optional integer: how many decimals to keep (e.g. 1).
#’ @param workers    Integer. Parallel workers (default = half of cores).
#’
#’ @return Invisibly, a data.table of `{date, cog, url, ok}`; also writes:
#’ * `out_dir/{dataset_id}_index.parquet`
#’ * `out_dir/{dataset_id}_log.csv`
#’
#’ @export
#’ @param timeout Integer. Seconds to wait for each download.file attempt (default: 60).
# Core downloader for daily CHC GeoTIFFs → COGs (with timeout in system2) ####
#’ @title download_chc_daily
#’ @description
#’ Downloads a daily CHC raster series (COG-converts, crops, rounds) in parallel,
#’ with per-call timeouts.
#’ @param start_year,end_year,out_dir,varname,base_url,dataset_id,bbox,
#’   round_digits,workers as in get_chirps_v3()
#’ @param timeout Integer. Seconds to allow each GDAL or gdal_calc call.
#’ @return Invisible data.table of \{date, cog, url, ok\}; writes index + log.
download_chc_daily <- function(start_year, end_year, out_dir,
                               varname, base_url, dataset_id,
                               bbox = NULL, round_digits = NULL,
                               workers = parallel::detectCores() %/% 2,
                               timeout = 60) {
  message("→ [download_chc_daily] dataset=", dataset_id,
          " years=", start_year, "-", end_year,
          " out_dir=", out_dir)

  # install / load deps
  pkgs <- c("data.table","fs","glue","future.apply","progressr","arrow")
  for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  lapply(pkgs, library, character.only=TRUE)

  stopifnot(start_year <= end_year)
  fs::dir_create(out_dir, recurse=TRUE)
  dates <- seq(as.Date(sprintf("%d-01-01", start_year)),
               as.Date(sprintf("%d-12-31", end_year)), by="day")
  start_time <- Sys.time()

  one_day <- function(d) {
    y  <- format(d, "%Y"); m <- format(d, "%m"); dd <- format(d, "%d")
    fname   <- glue::glue("{varname}.{y}.{m}.{dd}.tif")
    src     <- glue::glue("{base_url}/{y}/{fname}")
    src_vsi <- paste0("/vsicurl/", src)
    dest_dir<- fs::path(out_dir, y); fs::dir_create(dest_dir)
    dest    <- fs::path(dest_dir, fname)

    # skip existing
    if (fs::file_exists(dest) && fs::file_info(dest)$size > 1e4) {
      return(list(date=d, cog=dest, url=src, ok=TRUE))
    }

    # build bbox option
    bb_opt <- if (is.null(bbox)) character() else {
      c("-projwin",
        as.character(bbox[1]), as.character(bbox[4]),
        as.character(bbox[2]), as.character(bbox[3]))
    }

    tmp1 <- tempfile(fileext = ".tif")
    res <- tryCatch({
      # 1) streaming crop
      system2("gdal_translate",
              args = c("-q", bb_opt, shQuote(src_vsi), shQuote(tmp1)),
              stdout = NULL, stderr = NULL, timeout = timeout)

      # fallback if too small
      if (!fs::file_exists(tmp1) || fs::file_info(tmp1)$size < 1e3) {
        old_to <- getOption("timeout"); options(timeout = timeout)
        ok_dl <- FALSE
        for (i in 1:3) {
          try(suppressWarnings(download.file(src, tmp1, mode="wb", quiet=TRUE)), silent=TRUE)
          if (!is.null(bbox)) {
            td <- tempfile(fileext=".tif")
            system2("gdal_translate",
                    args = c("-q", bb_opt, shQuote(tmp1), shQuote(td)),
                    stdout = NULL, stderr = NULL, timeout = timeout)
            file.rename(td, tmp1)
          }
          ok_dl <- fs::file_exists(tmp1) && fs::file_info(tmp1)$size > 1e3
          if (ok_dl) break
          Sys.sleep(1.5 * i)
        }
        options(timeout = old_to)
        if (!ok_dl) stop("download failed after 3 attempts")
      }

      # 2) optional rounding
      tmp2 <- tmp1
      if (!is.null(round_digits)) {
        tmp2 <- tempfile(fileext = ".tif")
        scale <- 10^round_digits
        calc  <- sprintf("round(A*%d)/%d", scale, scale)
        system2("gdal_calc.py",
                args = c("-A", shQuote(tmp1),
                         "--outfile", shQuote(tmp2),
                         "--calc", shQuote(calc),
                         "--NoDataValue=-9999", "--quiet"),
                stdout = NULL, stderr = NULL, timeout = timeout)
      }

      # 3) COG conversion
      # 3) convert to COG
      cmd_args <- c("-q", "-a_nodata", "-9999",
                    "-of", "COG", "-co", "COMPRESS=ZSTD",
                    "-co", "BLOCKSIZE=512", "-co", "NUM_THREADS=1",
                    shQuote(tmp2), shQuote(dest))

      ok <- FALSE
      for (k in 1:3) {
        ret <- system2("gdal_translate", args = cmd_args,
                       stdout = NULL, stderr = NULL, timeout = timeout)
        if (ret == 0 &&
            fs::file_exists(dest) &&
            fs::file_info(dest)$size > 1e4) {
          ok <- TRUE ; break
        }
        Sys.sleep(1.5 * k)
      }
      if (!ok) stop("COG conversion failed (exit code ", ret, ")")
      list(date=d, cog=dest, url=src, ok=ok)

    }, error = function(e) {
      message("✖ [", dataset_id, "] ", fname, ": ", e$message)
      list(date=d, cog=dest, url=src, ok=FALSE)
    })

    res
  }

  # parallel + progressr
  future::plan(multisession, workers = workers)
  handlers(global = TRUE); handlers("progress")
  results <- progressr::with_progress({
    p <- progressr::progressor(length(dates))
    future.apply::future_lapply(dates, function(dd) { r <- one_day(dd); p(); r },
                                future.seed = TRUE)
  })
  future::plan(sequential)

  # write index
  dt <- data.table::rbindlist(results)
  idx <- fs::path(out_dir, paste0(dataset_id, "_index.parquet"))
  arrow::write_parquet(dt, idx, compression = "zstd")

  # write log
  log <- data.table::data.table(
    timestamp = Sys.time(),
    years     = paste(start_year, end_year, sep="-"),
    bbox      = if (is.null(bbox)) NA_character_ else paste(bbox, collapse=" "),
    round_digits = if (is.null(round_digits)) NA_integer_ else round_digits,
    workers   = workers,
    files_total = nrow(dt),
    files_ok    = sum(dt$ok),
    files_failed= sum(!dt$ok),
    out_dir   = normalizePath(out_dir),
    dataset_id= dataset_id,
    base_url  = base_url,
    timeout   = timeout,
    duration  = as.numeric(difftime(Sys.time(), start_time, units="secs")),
    user      = Sys.info()[["user"]],
    host      = Sys.info()[["nodename"]]
  )
  logf <- fs::path(out_dir, paste0(dataset_id, "_log.csv"))
  if (fs::file_exists(logf)) {
    old <- data.table::fread(logf)
    data.table::fwrite(rbind(old, log, fill=TRUE), logf)
  } else {
    data.table::fwrite(log, logf)
  }

  message(sprintf("✔ %s: %d/%d OK → index @ %s",
                  dataset_id, sum(dt$ok), nrow(dt), idx))
  invisible(dt)
}

# Wrappers for daily datasets ####

## chirps v3 ####

#’ @rdname download_chc_daily
#’ @export
get_chirps_v3 <- function(start_year, end_year, out_dir,
                          bbox=NULL, round_digits=1, workers=NULL,timeout=60) {
  message("→ [get_chirps_v3] called for years ", start_year, "–", end_year,
          "; out_dir = ", out_dir,
          "; bbox = ", if (is.null(bbox)) "NULL" else paste(bbox, collapse = " "),
          "; workers = ", if (is.null(workers)) "default" else workers,
          "; timeout = ", timeout, "\n")
  download_chc_daily(start_year, end_year, out_dir,
                     varname     = "chirps-v3.0",
                     base_url    = "https://data.chc.ucsb.edu/experimental/CHIRPS/v3.0/daily/final/ERA5",
                     dataset_id  = "chirps_v3",
                     bbox        = bbox,
                     round_digits= round_digits,
                     workers     = workers %||% parallel::detectCores()%/%2,
                     timeout = timeout)
}

## chirts era5 tmax ####
#’ @rdname download_chc_daily
#’ @export
get_chirts_era5_tmax <- function(start_year, end_year, out_dir,
                                 bbox=NULL, round_digits=NULL, workers=NULL,timeout=60) {
  download_chc_daily(start_year,
                     end_year,
                     out_dir,
                     varname     = "CHIRTS-ERA5.daily_Tmax",
                     base_url    = "https://data.chc.ucsb.edu/experimental/CHIRTS-ERA5/tmax/tifs/daily",
                     dataset_id  = "chirts_era5_tmax",
                     bbox        = NULL,
                     round_digits= NULL,
                     workers     = NULL,
                     timeout = 60)
}

## chirts era5 tmin ####
#’ @rdname download_chc_daily
#’ @export
get_chirts_era5_tmin <- function(start_year, end_year, out_dir,
                                 bbox=NULL, round_digits=NULL, workers=NULL,timeout=60) {
  download_chc_daily(start_year, end_year, out_dir,
                     varname     = "CHIRTS-ERA5.daily_Tmin",
                     base_url    = "https://data.chc.ucsb.edu/experimental/CHIRTS-ERA5/tmin/tifs/daily",
                     dataset_id  = "chirts_era5_tmin",
                     bbox        = bbox,
                     round_digits= round_digits,
                     workers     = workers %||% parallel::detectCores()%/%2,
                     timeout = timeout)
}

## chirts era5 vpd ####
#’ @rdname download_chc_daily - files do not exist at present,
#’ @export
get_chirts_era5_vpd <- function(start_year, end_year, out_dir,
                                bbox=NULL, round_digits=NULL, workers=NULL,timeout=60) {
  download_chc_daily(start_year, end_year, out_dir,
                     varname     = "CHIRTS-ERA5.daily_vpd",
                     base_url    = "https://data.chc.ucsb.edu/experimental/CHC_CMIP6/observations/vpd",
                     dataset_id  = "chc_cmip6_observations",
                     bbox        = NULL,
                     round_digits= 1,
                     workers     = workers %||% parallel::detectCores()%/%2,
                     timeout = timeout)
}

## hobbins ref et ####
#' @title Download Hobbins Reference ET (dekadal) → COGs (Corrected)
#'
#' @description
#' Downloads dekadal Hobbins reference evapotranspiration (ETos_p05)
#' from CHC, uncompresses .tif.gz correctly, optionally crops, COG-converts,
#' and writes an index + log.
#'
#' @param start_year Integer. First year to download.
#' @param end_year   Integer. Last year to download.
#' @param out_dir    Character. Folder to write COGs (subdirs per year).
#' @param bbox       Optional numeric(4): c(xmin, ymin, xmax, ymax) for `-projwin`.
#' @param round_digits Optional integer: how many decimals to round (e.g., 1).
#' @param workers    Integer. Number of parallel workers.
#'
get_hobbins_refet <- function(start_year,
                              end_year,
                              out_dir,
                              dataset_id="CHC-products-Hobbins_RefET",
                              base_url="https://data.chc.ucsb.edu/products/Hobbins_RefET/ETos_p05_dekad_global/tifs",
                              workers,
                              bbox = NULL,
                              round_digits= 1,
                              timeout = 60) {
  message("→ [download_hobbins_refet] dataset=", dataset_id,
          " years=", start_year, "-", end_year,
          " out_dir=", out_dir)

  pkgs <- c("data.table","fs","glue","future.apply","progressr","arrow")
  for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  invisible(lapply(pkgs, library, character.only=TRUE))

  fs::dir_create(out_dir, recurse=TRUE)
  dekads <- expand.grid(year = start_year:end_year, dekad = 1:36)
  start_time <- Sys.time()

  one_dekad <- function(yr, dk, base_url,timeout) {
    fname <- sprintf("ETos_p05.%d%02d.tif", yr, dk)
    url   <- glue::glue("{base_url}/{fname}")
    vsi   <- paste0("/vsicurl/", url)
    dest  <- fs::path(out_dir, as.character(yr), fname)
    fs::dir_create(fs::path_dir(dest))

    if (fs::file_exists(dest) && fs::file_info(dest)$size > 1e4) {
      return(list(year=yr, dekad=dk, cog=dest, url=url, ok=TRUE))
    }

    bb_opt <- if (is.null(bbox)) character() else {
      c("-projwin",
        as.character(bbox[1]), as.character(bbox[4]),
        as.character(bbox[3]), as.character(bbox[2]))
    }

    tmp1 <- tempfile(fileext = ".tif")
    res <- tryCatch({
      rc <- system2("gdal_translate",
                    args = c("-q", bb_opt, shQuote(vsi), shQuote(tmp1)),
                    stdout = NULL, stderr = NULL, timeout = timeout)

      good_stream <- rc == 0 &&
        file.exists(tmp1) &&
        fs::file_info(tmp1)$size > 1024

      if (!good_stream) {

        ok_dl <- FALSE
        for (i in 1:3) {                       # ➊  three attempts
          # -- download the full TIFF -------------------------------------------
          try(suppressWarnings(
            download.file(url, tmp1, mode = "wb", quiet = TRUE)
          ), silent = TRUE)

          # -- crop locally if bbox requested ------------------------------------
          if (!is.null(bbox) && file.exists(tmp1)) {
            td  <- tempfile(fileext = ".tif")
            rc2 <- system2("gdal_translate",
                           args = c("-q", bb_opt, shQuote(tmp1), shQuote(td)),
                           stdout = NULL, stderr = NULL, timeout = timeout)

            if (rc2 == 0 && file.exists(td) && fs::file_info(td)$size > 1024) {
              fs::file_copy(td, tmp1, overwrite = TRUE)
              fs::file_delete(td)
            }
          }

          ok_dl <- file.exists(tmp1) && fs::file_info(tmp1)$size > 1024
          if (ok_dl) break                    # ➋  success: leave the loop
          Sys.sleep(1.5 * i)                  # ➌  back-off before next try
        }

        if (!ok_dl)
          stop("download failed after 3 attempts")
      }

      # optional rounding
      tmp2 <- tmp1
      if (!is.null(round_digits)) {
        tmp2 <- tempfile(fileext = ".tif")
        scale <- 10^round_digits
        calc  <- sprintf("round(A*%d)/%d", scale, scale)
        system2("gdal_calc.py",
                args = c("-A", shQuote(tmp1),
                         "--outfile", shQuote(tmp2),
                         "--calc", shQuote(calc),
                         "--NoDataValue=-9999", "--quiet"),
                stdout = NULL, stderr = NULL, timeout = timeout)
      }

      cmd_args <- c("-q", "-a_nodata", "-9999",
                    "-of", "COG", "-co", "COMPRESS=ZSTD",
                    "-co", "BLOCKSIZE=512", "-co", "NUM_THREADS=1",
                    shQuote(tmp2), shQuote(dest))

      ok <- FALSE
      for (k in 1:3) {
        ret <- system2("gdal_translate", args = cmd_args,
                       stdout = NULL, stderr = NULL, timeout = timeout)
        if (ret == 0 &&
            fs::file_exists(dest) &&
            fs::file_info(dest)$size > 1e4) {
          ok <- TRUE ; break
        }
        Sys.sleep(1.5 * k)
      }
      if (!ok) stop("COG conversion failed")
      list(year=yr, dekad=dk, cog=dest, url=url, ok=TRUE)

    }, error = function(e) {
      message("✖ [", dataset_id, "] ", fname, ": ", e$message)
      list(year=yr, dekad=dk, cog=NA_character_, url=url, ok=FALSE)
    })

    res
  }

  future::plan(multisession, workers = future::availableCores() %/% 2)
  handlers(global = TRUE); handlers("progress")
  results <- progressr::with_progress({
    p <- progressr::progressor(nrow(dekads))
    future.apply::future_lapply(1:nrow(dekads), function(i) {
      row <- dekads[i, ]
      r <- one_dekad(yr=row$year, dk=row$dekad,base_url=base_url,timeout=timeout); p(); r
    }, future.seed = TRUE)
  })
  future::plan(sequential)

  dt <- data.table::rbindlist(results, ignore.attr=TRUE)
  idx <- fs::path(out_dir, paste0(dataset_id, "_index.parquet"))
  arrow::write_parquet(dt, idx, compression = "zstd")

  log <- data.table::data.table(
    timestamp = Sys.time(),
    years     = paste(start_year, end_year, sep="-"),
    bbox      = if (is.null(bbox)) NA_character_ else paste(bbox, collapse=" "),
    workers   = future::availableCores() %/% 2,
    files_total = nrow(dt),
    files_ok    = sum(dt$ok),
    files_failed= sum(!dt$ok),
    out_dir   = normalizePath(out_dir),
    dataset_id= dataset_id,
    base_url  = base_url,
    timeout   = timeout,
    duration  = as.numeric(difftime(Sys.time(), start_time, units="secs")),
    user      = Sys.info()[["user"]],
    host      = Sys.info()[["nodename"]]
  )
  logf <- fs::path(out_dir, paste0(dataset_id, "_log.csv"))
  if (fs::file_exists(logf)) {
    old <- data.table::fread(logf)
    data.table::fwrite(rbind(old, log, fill=TRUE), logf)
  } else {
    data.table::fwrite(log, logf)
  }

  message(sprintf("✔ %s: %d/%d OK → index @ %s",
                  dataset_id, sum(dt$ok), nrow(dt), idx))
  invisible(dt)
}
