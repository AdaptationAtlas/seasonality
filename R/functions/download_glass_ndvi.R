#' @title Parallel downloader for GLASS NDVI 0.05 degree products
#' @description Download .hdf files from GLASS HKU server in parallel
#' @param year Year to download (e.g. 2000)
#' @param out_dir Output directory
#' @param overwrite Logical, whether to overwrite existing files
#' @param workers Number of parallel workers
#' @return data.table with columns: file, url, dest, ok, error
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor handlers
#' @importFrom data.table data.table rbindlist fwrite
#' @importFrom fs dir_create file_exists
#' @export
download_glass_ndvi <- function(year = 2000,
                                         out_dir = "data/glass_005_ndvi",
                                         overwrite = FALSE,
                                         workers = future::availableCores() %/% 2) {

  # Load required packages
  pkgs <- c("fs", "httr", "rvest", "xml2", "progressr", "future.apply", "data.table")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  invisible(lapply(pkgs, library, character.only = TRUE))

  base_url <- sprintf("http://glass.hku.hk/archive/NDVI/MODIS/0.05/%d/", year)
  fs::dir_create(out_dir, recurse = TRUE)

  # Scrape file names
  page <- xml2::read_html(base_url)
  files <- page %>%
    rvest::html_elements("a") %>%
    rvest::html_attr("href") %>%
    grep("\\.hdf$", ., value = TRUE)

  urls <- paste0(base_url, files)
  dests <- file.path(out_dir, files)

  # Define individual download task
  download_one <- function(i) {
    url <- urls[i]
    dest <- dests[i]
    file <- files[i]
    attempts <- 0
    ok <- FALSE
    err <- NA_character_

    if (!fs::file_exists(dest) || overwrite) {
      while (attempts < 3 && !ok) {
        attempts <- attempts + 1
        tryCatch({
          httr::GET(url, httr::write_disk(dest, overwrite = TRUE), timeout(60))
          ok <- fs::file_exists(dest) && fs::file_info(dest)$size > 1024
        }, error = function(e) {
          err <<- e$message
          Sys.sleep(1.5 * attempts)
        })
      }
    } else {
      ok <- TRUE
    }
    data.table::data.table(file = file, url = url, dest = dest, ok = ok, error = ifelse(ok, NA, err))
  }

  # Run in parallel with progress
  future::plan(multisession, workers = workers)
  progressr::handlers(global = TRUE)

  results <- progressr::with_progress({
    p <- progressr::progressor(along = urls)
    future.apply::future_lapply(seq_along(urls), function(i) {
      res <- download_one(i); p(); res
    })
  })

  future::plan(sequential)

  log_file<-file.path(out_dir, "glass_download_log.csv")
  log <- data.table::rbindlist(results)
  log$dat_year<-year

  message(sprintf("✔ GLASS NDVI %d: %d/%d files OK.", year, sum(log$ok), nrow(log)))

  if(file.exists(log_file)){
    log_ex<-data.table::fread(log_file)
    log_ex<-log_ex[dat_year!=year]
    log<-rbind(log_ex,log)
  }

  data.table::fwrite(log, log_file)

  return(invisible(log))
}
