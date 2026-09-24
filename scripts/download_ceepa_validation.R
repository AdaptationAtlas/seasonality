#!/usr/bin/env Rscript

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("download-ceepa-validation", config$country$iso3)
log_event(run, event = "run_started", message = "Downloading CEEPA validation archive")

failed <- TRUE
on.exit({
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
ceepa <- config$validation$ceepa
target_dir <- file.path(data_root, "static_raw", "validation", "ceepa")
archive <- file.path(target_dir, "CEEPA.zip")
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(archive) || !identical(unname(tools::md5sum(archive)), ceepa$md5)) {
  download.file(ceepa$download_url, archive, mode = "wb", quiet = FALSE)
}

actual_md5 <- unname(tools::md5sum(archive))
if (!identical(actual_md5, ceepa$md5)) {
  stop("CEEPA archive checksum mismatch: ", actual_md5)
}

utils::unzip(
  archive,
  files = c("CEEPASurvey.dta", "CEEPASurvey.txt", "Questionnaire.pdf", "SurveyManual.pdf"),
  exdir = target_dir,
  overwrite = TRUE
)

meta <- list(
  source = "CEEPA African agricultural survey",
  article_id = ceepa$article_id,
  doi = ceepa$doi,
  download_url = ceepa$download_url,
  md5 = actual_md5,
  downloaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
jsonlite::write_json(
  meta,
  file.path(target_dir, "source_metadata.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

failed <- FALSE
finish_run(run, "success", list(target_dir = target_dir, md5 = actual_md5))
cat("CEEPA validation source ready:", target_dir, "\n")
