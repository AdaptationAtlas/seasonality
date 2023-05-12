#' Downloads and unzips a file
#'
#' This function downloads a file from a specified URL and optionally unzips it to a destination folder.
#'
#' @param url character string specifying the URL of the file to be downloaded
#' @param filename character string specifying the name of the file to be saved.
#' @param destname character string specifying the name of the output file after unzipping.
#' @param quiet logical indicating whether to display download progress messages.
#' @param unzip logical indicating whether to unzip the downloaded file, default is TRUE.
#' @param remove logical indicating whether to remove the .gz file after unzipping, default is TRUE.
#'
#' @return NULL
#'
#' @examples
#' # Download and unzip a file
#' url <- "https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_0.8.0.tar.gz"
#' download_and_unzip(url, "dplyr_0.8.0.tar.gz", "dplyr")
#'
#' @import R.utils
#' @importFrom utils download.file
#' @export

download_and_unzip <- function(url, filename, destname, quiet = TRUE, unzip = T, remove = T) {

  # Download the .gz file
  download.file(url, destfile = filename, mode = "wb", quiet = quiet)

  # Unzip the .gz file
  if (unzip) {
    R.utils::gunzip(filename = filename, destname = destname, overwrite = TRUE, remove = remove)
  }

}
