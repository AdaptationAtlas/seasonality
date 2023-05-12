#' Download and Unzip Function
#'
#' This function downloads a .gz file from an input URL and then unzips the file to a destination folder.
#' The compressed file is then removed from the system.
#'
#' @param url A character vector specifying the url address of the .gz file to download.
#' @param destination A character vector specifying the location and filename where the downloaded and unzipped file should be saved.
#' @param quiet A logical parameter to suppress messages or not during download (default=FALSE)
#'
#' @return Returns the unzipped file.
#'
#' @examples
#' download_and_unzip("https://example.com/data.gz", "downloads/data.csv")
#'
#' @importFrom tools gunzip
#' @importFrom base download.file
#'
#' @export
download_and_unzip <- function(url, destination,quiet=T) {
  # Download the .gz file
  download.file(url, destfile = destination, mode = "wb",quiet=quiet)

  # Unzip the .gz file
  unzipped_file <- tools::gunzip(destination, overwrite = TRUE)

  # Remove the compressed file
  file.remove(destination)

  return(unzipped_file)
}
