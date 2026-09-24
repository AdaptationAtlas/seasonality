#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(haven)
})

source("R/functions/project_logging.R")
source("R/functions/ceepa_validation.R")

config <- read_project_config()
run <- new_run_context("kenya-ceepa-validation", config$country$iso3)
log_event(run, event = "run_started", message = "Building Kenya CEEPA crop-calendar validation data")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
source_file <- file.path(data_root, "static_raw", "validation", "ceepa", "CEEPASurvey.dta")
output_dir <- file.path(
  data_root, "models", "validation_reports", "kenya_detectability", "independent_validation"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(source_file)) stop("Missing CEEPA source: run scripts/download_ceepa_validation.R")

survey <- as.data.table(read_dta(source_file))
survey <- survey[as.numeric(adm0) == 6]
if (nrow(survey) != 816L) stop("Expected 816 Kenya CEEPA households; found ", nrow(survey))

decode <- function(x) as.character(haven::as_factor(x, levels = "labels"))
survey[, `:=`(
  admin1_legacy = decode(adm1),
  admin2_legacy = decode(adm2),
  subdivision = decode(admindiv),
  village = decode(adminvil),
  season_1_name = decode(seas1nam),
  season_2_name = decode(seas2nam)
)]
survey[, admin1_name := ceepa_current_county(admin2_legacy)]
if (anyNA(survey$admin1_name)) {
  stop("Unmapped CEEPA districts: ", paste(unique(survey[is.na(admin1_name), admin2_legacy]), collapse = ", "))
}

slot_names <- grep("^s[123]p[123]c[1-6]$", names(survey), value = TRUE)
observations <- rbindlist(lapply(slot_names, function(slot) {
  match_parts <- regmatches(slot, regexec("^s([123])p([123])c([1-6])$", slot))[[1]]
  plant_name <- paste0(slot, "plant")
  harvest_name <- paste0(slot, "harv")
  area_name <- paste0(slot, "area")
  data.table(
    hhcode = survey$hhcode,
    admin1_name = survey$admin1_name,
    admin1_legacy = survey$admin1_legacy,
    admin2_legacy = survey$admin2_legacy,
    subdivision = survey$subdivision,
    village = survey$village,
    survey_season = as.integer(match_parts[2]),
    survey_season_name = if (match_parts[2] == "1") survey$season_1_name else if (match_parts[2] == "2") survey$season_2_name else NA_character_,
    plot = as.integer(match_parts[3]),
    crop_position = as.integer(match_parts[4]),
    crop_code = as.integer(survey[[slot]]),
    planting_raw = decode(survey[[plant_name]]),
    harvest_raw = decode(survey[[harvest_name]]),
    plot_area_percent = as.numeric(survey[[area_name]])
  )
}), use.names = TRUE)

crop_lookup <- as.data.table(ceepa_crop_lookup())
observations <- crop_lookup[observations, on = "crop_code"]
observations[, survey_season_name := trimws(tolower(survey_season_name))]
observations[grepl("^long rain", survey_season_name), survey_season_name := "long rains"]
observations[grepl("^short rain", survey_season_name), survey_season_name := "short rains"]
observations[survey_season_name == "", survey_season_name := NA_character_]

planting <- as.data.table(parse_ceepa_calendar_dates(observations$planting_raw))
setnames(planting, names(planting), paste0("planting_", names(planting)))
harvest <- as.data.table(parse_ceepa_calendar_dates(observations$harvest_raw))
setnames(harvest, names(harvest), paste0("harvest_", names(harvest)))
observations <- cbind(observations, planting, harvest)

observations <- observations[!is.na(crop_code) | planting_parse_status != "missing" | harvest_parse_status != "missing"]
observations[, planting_validation_eligible :=
  !is.na(crop_cycle) & crop_cycle == "annual" &
  planting_parse_status == "parsed" &
  (is.na(planting_year) | planting_year %in% 2003:2004)
]
observations[, season_length_days := fifelse(
  planting_validation_eligible & harvest_parse_status == "parsed",
  ((harvest_doy - planting_doy) %% 365L),
  NA_integer_
)]
observations[season_length_days < 20L | season_length_days > 365L, season_length_days := NA_integer_]

summary <- observations[
  planting_validation_eligible == TRUE,
  .(
    observations = .N,
    households = uniqueN(hhcode),
    planting_doy_q10 = as.numeric(quantile(planting_doy, 0.1, na.rm = TRUE)),
    planting_doy_median = as.numeric(median(planting_doy, na.rm = TRUE)),
    planting_doy_q90 = as.numeric(quantile(planting_doy, 0.9, na.rm = TRUE)),
    exact_date_share = mean(planting_precision == "day", na.rm = TRUE),
    season_length_observations = sum(!is.na(season_length_days)),
    season_length_median = as.numeric(median(season_length_days, na.rm = TRUE))
  ),
  by = .(admin1_name, admin2_legacy, survey_season, survey_season_name, crop_code, crop_name)
]
summary[!is.finite(season_length_median), season_length_median := NA_real_]

write_parquet(observations, file.path(output_dir, "KEN_ceepa_crop_observations.parquet"))
write_parquet(summary, file.path(output_dir, "KEN_ceepa_crop_calendar_summary.parquet"))

qa <- observations[, .(
  observations = .N,
  parsed = sum(planting_parse_status == "parsed"),
  eligible_annual = sum(planting_validation_eligible, na.rm = TRUE),
  unparsed = sum(planting_parse_status == "unparsed"),
  continuous = sum(planting_parse_status == "continuous")
)]
fwrite(qa, file.path(output_dir, "KEN_ceepa_parse_qa.csv"))

failed <- FALSE
finish_run(
  run,
  "success",
  list(
    households = nrow(survey),
    crop_observations = nrow(observations),
    eligible_annual = qa$eligible_annual,
    output_dir = output_dir
  )
)
cat("Wrote CEEPA Kenya validation products under", output_dir, "\n")
