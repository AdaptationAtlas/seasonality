# Parse GLASS filenames containing AYYYYDDD acquisition tokens.

parse_glass_ndvi_dates <- function(files) {
  names_only <- basename(files)
  tokens <- sub(".*\\.A([0-9]{7})\\..*", "\\1", names_only)
  valid <- grepl("^[0-9]{7}$", tokens)
  if (!all(valid)) {
    stop("Could not parse GLASS AYYYYDDD date from: ", paste(names_only[!valid], collapse = ", "))
  }

  years <- as.integer(substr(tokens, 1, 4))
  days <- as.integer(substr(tokens, 5, 7))
  dates <- as.Date(paste0(years, "-01-01")) + days - 1L

  if (any(as.integer(format(dates, "%Y")) != years)) {
    stop("GLASS day-of-year falls outside filename year.")
  }
  dates
}
