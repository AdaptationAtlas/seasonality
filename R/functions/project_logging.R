# Structured runtime logging for seasonality pipeline stages.

project_root <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

read_project_config <- function(path = file.path(project_root(), "config", "project.json")) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read project configuration.")
  }
  if (!file.exists(path)) stop("Project config not found: ", path)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

new_run_context <- function(stage, country = "KEN", log_dir = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for structured logging.")
  }
  stopifnot(is.character(stage), length(stage) == 1L, nzchar(stage))

  if (is.null(log_dir)) {
    log_dir <- Sys.getenv(
      "SEASONALITY_LOG_DIR",
      unset = file.path(project_root(), "logs")
    )
  }
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  started_at <- Sys.time()
  run_id <- sprintf(
    "%s_%s_%s_%s",
    format(started_at, "%Y%m%dT%H%M%S"),
    country,
    gsub("[^A-Za-z0-9_-]", "-", stage),
    Sys.getpid()
  )

  structure(
    list(
      run_id = run_id,
      stage = stage,
      country = country,
      started_at = started_at,
      log_file = file.path(log_dir, paste0(run_id, ".jsonl"))
    ),
    class = "seasonality_run"
  )
}

log_event <- function(run, level = "INFO", event, message, data = list()) {
  stopifnot(inherits(run, "seasonality_run"))
  record <- c(
    list(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
      level = toupper(level),
      run_id = run$run_id,
      stage = run$stage,
      country = run$country,
      event = event,
      message = message
    ),
    data
  )
  line <- jsonlite::toJSON(record, auto_unbox = TRUE, null = "null", na = "null")
  cat(line, "\n", file = run$log_file, append = TRUE, sep = "")
  invisible(record)
}

finish_run <- function(run, status = "success", data = list()) {
  elapsed <- as.numeric(difftime(Sys.time(), run$started_at, units = "secs"))
  log_event(
    run,
    level = if (identical(status, "success")) "INFO" else "ERROR",
    event = "run_finished",
    message = paste("Run finished with status", status),
    data = c(list(status = status, elapsed_seconds = elapsed), data)
  )
}
