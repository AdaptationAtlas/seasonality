#!/usr/bin/env Rscript

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("preflight", config$country$iso3)
log_event(run, event = "run_started", message = "Preflight checks started")

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
checks <- list(
  data_root = dir.exists(data_root),
  kenya_phenology = file.exists(file.path(
    data_root, "climate_derived", "glass_phenology", "countries",
    "KEN_seasonal-phenology_plus-rain.parquet"
  )),
  pixel_index = file.exists(file.path(
    data_root, "climate_derived", "glass_phenology", "pixel_index.parquet"
  )),
  admin1_boundaries = file.exists(file.path(
    data_root, "static_raw", "boundaries",
    "atlas_gaul24_a1_africa_simple-highres.parquet"
  )),
  kenya_chirps = file.exists(file.path(
    data_root, "climate_raw", "chirps", "chirps_v3_cog_countries",
    "KEN.parquet"
  )),
  gdal_translate = nzchar(Sys.which("gdal_translate"))
)

required_packages <- c(
  "arrow", "data.table", "DBI", "duckdb", "jsonlite", "lubridate", "terra"
)
package_checks <- setNames(
  vapply(required_packages, requireNamespace, logical(1), quietly = TRUE),
  paste0("package_", required_packages)
)
checks <- c(checks, as.list(package_checks))

admin2_file <- file.path(
  data_root, "static_raw", "boundaries",
  "atlas_gaul24_a2_africa_simple-highres.parquet"
)
admin2_valid <- FALSE
admin2_error <- NULL
if (file.exists(admin2_file) && requireNamespace("arrow", quietly = TRUE)) {
  admin2_valid <- tryCatch({
    arrow::read_parquet(admin2_file, as_data_frame = FALSE)
    TRUE
  }, error = function(e) {
    admin2_error <<- conditionMessage(e)
    FALSE
  })
}

for (nm in names(checks)) {
  ok <- isTRUE(checks[[nm]])
  log_event(
    run,
    level = if (ok) "INFO" else "ERROR",
    event = "check",
    message = nm,
    data = list(ok = ok)
  )
  cat(sprintf("%-28s %s\n", nm, if (ok) "OK" else "FAIL"))
}

log_event(
  run,
  level = if (admin2_valid) "INFO" else "WARN",
  event = "optional_check",
  message = "admin2_boundaries",
  data = list(ok = admin2_valid, error = admin2_error)
)
cat(sprintf("%-28s %s\n", "admin2_boundaries", if (admin2_valid) "OK" else "WARN"))

failed <- names(checks)[!vapply(checks, isTRUE, logical(1))]
if (length(failed)) {
  finish_run(run, "failed", list(failed_checks = failed))
  stop("Preflight failed: ", paste(failed, collapse = ", "), call. = FALSE)
}

finish_run(run, "success", list(admin2_ready = admin2_valid))
cat("Preflight passed. Runtime log:", run$log_file, "\n")
