# Helpers for CEEPA household-survey crop calendar validation.

ceepa_crop_lookup <- function() {
  crop_name <- c(
    "alfalfa", "banana", "barley", "beans", "cashew", "cassava",
    "citrus fruit", "chickpeas", "clover", "cocoa", "cocoyam", "cowpea",
    "coffee", "cotton", "cucumber", "enset", "field pea", "flax",
    "garden-eggs", "garlic", "grape", "groundnut", "kola", "lentil",
    "mango", "maize", "millet", "oil palm", "okra", "onion",
    "palm dates", "paprika", "peanuts", "pepper", "pigeon pea", "pineapple",
    "plantain", "potato", "rice", "safflower", "sesame", "shallots",
    "sheanut", "sorghum", "soybean", "spinach", "squash", "sugarcane",
    "sunflower", "tea", "tef", "tobacco", "tomato", "wheat", "yam",
    "other"
  )
  perennial <- c(
    "alfalfa", "banana", "cashew", "citrus fruit", "cocoa", "coffee",
    "enset", "grape", "kola", "mango", "oil palm", "palm dates",
    "pineapple", "plantain", "sheanut", "sugarcane", "tea"
  )
  uncertain <- c("cassava", "clover", "cocoyam", "pepper", "pigeon pea", "yam", "other")
  data.frame(
    crop_code = seq_along(crop_name),
    crop_name = crop_name,
    crop_cycle = ifelse(
      crop_name %in% perennial, "perennial",
      ifelse(crop_name %in% uncertain, "uncertain", "annual")
    ),
    stringsAsFactors = FALSE
  )
}

ceepa_current_county <- function(district) {
  cleaned <- toupper(trimws(district))
  cleaned <- gsub("[^A-Z]", "", cleaned)
  lookup <- c(
    BARINGO = "Baringo", BOMET = "Bomet", BUNGOMA = "Bungoma",
    BUSIA = "Busia", ELGEYOMARAKWET = "Elgeyo-Marakwet", EMBU = "Embu",
    HOMABAY = "Homa Bay", KAJIADO = "Kajiado", KAKAMEGA = "Kakamega",
    KERICHO = "Kericho", KIAMBU = "Kiambu", KILIFI = "Kilifi",
    KIRINYAGA = "Kirinyaga", KISII = "Kisii", KISUMU = "Kisumu",
    KITUI = "Kitui", KWALE = "Kwale", LAIKIPIA = "Laikipia",
    MACHAKOS = "Machakos", MAKUENI = "Makueni", MERU = "Meru",
    MERUNORTH = "Meru", NYAMBENE = "Meru", MIGORI = "Migori",
    MURANGA = "Murang'A", NAKURU = "Nakuru", NANDI = "Nandi",
    NAROK = "Narok", NITHI = "Tharaka-Nithi", NYAMIRA = "Nyamira",
    NYANDARUA = "Nyandarua", NYERI = "Nyeri", SIAYA = "Siaya",
    TAITATAVETA = "Taita Taveta", TANARIVER = "Tana River",
    TRNASNZOIA = "Trans Nzoia", TRANSNZOIA = "Trans Nzoia",
    UASINGISHU = "Uasin Gishu", VIHIGA = "Vihiga", WESTPOKOT = "West Pokot"
  )
  unname(lookup[cleaned])
}

parse_ceepa_calendar_dates <- function(x) {
  raw <- as.character(x)
  normalized <- tolower(trimws(raw))
  normalized <- gsub("\\s+", "", normalized)
  normalized <- gsub("sept", "sep", normalized, fixed = TRUE)
  normalized <- gsub("frb", "feb", normalized, fixed = TRUE)
  normalized <- gsub("mat", "mar", normalized, fixed = TRUE)

  missing_tokens <- c("", ".", "0", "-98", "-99", "-100", "-999")
  continuous <- grepl("contin|perennial", normalized)
  missing <- is.na(normalized) | normalized %in% missing_tokens

  n <- length(normalized)
  year <- rep(NA_integer_, n)
  month <- rep(NA_integer_, n)
  day <- rep(NA_integer_, n)
  precision <- rep(NA_character_, n)

  month_names <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
  month_pattern <- paste(month_names, collapse = "|")

  parse_year <- function(token) {
    value <- suppressWarnings(as.integer(token))
    recent <- !is.na(value) & nchar(token) == 2L & value <= 29L
    historic <- !is.na(value) & nchar(token) == 2L & value >= 30L
    value[recent] <- value[recent] + 2000L
    value[historic] <- value[historic] + 1900L
    value[token %in% c("yr", "yy", "")] <- NA_integer_
    value
  }

  exact_pattern <- paste0("^([0-9]{1,2})(", month_pattern, ")([0-9]{2,4}|yr|yy)$")
  exact <- !missing & !continuous & grepl(exact_pattern, normalized)
  if (any(exact)) {
    parts <- regmatches(normalized[exact], regexec(exact_pattern, normalized[exact]))
    day[exact] <- as.integer(vapply(parts, `[`, character(1), 2L))
    month[exact] <- match(vapply(parts, `[`, character(1), 3L), month_names)
    year[exact] <- parse_year(vapply(parts, `[`, character(1), 4L))
    precision[exact] <- "day"
  }

  hyphen_pattern <- paste0("^([0-9]{1,2})-?(", month_pattern, ")-?([0-9]{2,4})$")
  hyphen <- !missing & !continuous & is.na(precision) & grepl(hyphen_pattern, normalized)
  if (any(hyphen)) {
    parts <- regmatches(normalized[hyphen], regexec(hyphen_pattern, normalized[hyphen]))
    day[hyphen] <- as.integer(vapply(parts, `[`, character(1), 2L))
    month[hyphen] <- match(vapply(parts, `[`, character(1), 3L), month_names)
    year[hyphen] <- parse_year(vapply(parts, `[`, character(1), 4L))
    precision[hyphen] <- "day"
  }

  week_pattern <- paste0("^w([1-5])(", month_pattern, ")([0-9]{2,4}|yr|yy)?$")
  week <- !missing & !continuous & is.na(precision) & grepl(week_pattern, normalized)
  if (any(week)) {
    parts <- regmatches(normalized[week], regexec(week_pattern, normalized[week]))
    week_number <- as.integer(vapply(parts, `[`, character(1), 2L))
    day[week] <- c(4L, 11L, 18L, 25L, 28L)[week_number]
    month[week] <- match(vapply(parts, `[`, character(1), 3L), month_names)
    year[week] <- parse_year(vapply(parts, function(z) if (length(z) >= 4L) z[4L] else "", character(1)))
    precision[week] <- "week"
  }

  month_only_pattern <- paste0("^(?:dd|mon)?(", month_pattern, ")([0-9]{2,4}|yr|yy)?$")
  month_only <- !missing & !continuous & is.na(precision) & grepl(month_only_pattern, normalized)
  if (any(month_only)) {
    parts <- regmatches(normalized[month_only], regexec(month_only_pattern, normalized[month_only]))
    day[month_only] <- 15L
    month[month_only] <- match(vapply(parts, `[`, character(1), 2L), month_names)
    year[month_only] <- parse_year(vapply(parts, function(z) if (length(z) >= 3L) z[3L] else "", character(1)))
    precision[month_only] <- "month"
  }

  valid_calendar <- !is.na(month) & !is.na(day)
  reference_date <- rep(as.Date(NA), n)
  reference_date[valid_calendar] <- as.Date(sprintf(
    "2001-%02d-%02d", month[valid_calendar], day[valid_calendar]
  ))
  invalid_day <- valid_calendar & is.na(reference_date)
  precision[invalid_day] <- NA_character_

  actual_date <- rep(as.Date(NA), n)
  dated <- !is.na(year) & !is.na(month) & !is.na(day) & !invalid_day
  actual_date[dated] <- as.Date(sprintf("%04d-%02d-%02d", year[dated], month[dated], day[dated]))

  status <- ifelse(
    missing, "missing",
    ifelse(continuous, "continuous", ifelse(is.na(precision), "unparsed", "parsed"))
  )
  uncertainty <- ifelse(precision == "day", 0L, ifelse(precision == "week", 3L, ifelse(precision == "month", 15L, NA_integer_)))

  data.frame(
    raw = raw,
    normalized = normalized,
    date = actual_date,
    year = year,
    month = month,
    day = day,
    doy = as.integer(format(reference_date, "%j")),
    precision = precision,
    uncertainty_days = uncertainty,
    parse_status = status,
    stringsAsFactors = FALSE
  )
}
