#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(DBI)
  library(duckdb)
})

source("R/functions/project_logging.R")

config <- read_project_config()
run <- new_run_context("kenya-rainfall-metrics", config$country$iso3)
log_event(run, event = "run_started", message = "Building Kenya rainfall seasonality metrics")

failed <- TRUE
connection <- NULL
on.exit({
  if (!is.null(connection)) try(DBI::dbDisconnect(connection, shutdown = TRUE), silent = TRUE)
  if (failed) finish_run(run, "failed")
}, add = TRUE)

data_root <- Sys.getenv(config$data_root_env, unset = config$data_root_default)
rain_file <- file.path(
  data_root, "climate_raw", "chirps", "chirps_v3_cog_countries", "KEN.parquet"
)
pheno_file <- file.path(
  data_root, "climate_derived", "glass_phenology", "countries",
  "KEN_seasonal-phenology.parquet"
)
output_dir <- file.path(data_root, "climate_derived", "glass_phenology")
monthly_file <- file.path(output_dir, "KEN_monthly_rainfall.parquet")
annual_file <- file.path(output_dir, "KEN_annual_rainfall_metrics.parquet")
baseline_file <- file.path(output_dir, "KEN_baseline_monthly_rainfall.parquet")

if (!file.exists(rain_file)) stop("Kenya CHIRPS file not found: ", rain_file)
if (!file.exists(pheno_file)) stop("Kenya phenology file not found: ", pheno_file)

sql_path <- function(path) gsub("'", "''", normalizePath(path, mustWork = FALSE), fixed = TRUE)
rain_sql <- sql_path(rain_file)
monthly_sql <- sql_path(monthly_file)
baseline_sql <- sql_path(baseline_file)
start_date <- sprintf("%d-01-01", config$baseline$start_year)
end_date <- sprintf("%d-12-31", config$baseline$end_year)

connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbExecute(connection, "SET threads TO 8")

monthly_query <- sprintf(
  paste(
    "SELECT pixel,",
    "CAST(substr(date, 1, 4) AS INTEGER) AS year,",
    "CAST(substr(date, 6, 2) AS INTEGER) AS month,",
    "SUM(value) AS monthly_rain",
    "FROM read_parquet('%s')",
    "WHERE date BETWEEN '%s' AND '%s'",
    "GROUP BY pixel, year, month",
    sep = "\n"
  ),
  rain_sql, start_date, end_date
)

DBI::dbExecute(
  connection,
  sprintf(
    "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE_OR_IGNORE TRUE)",
    monthly_query, monthly_sql
  )
)
log_event(run, event = "monthly_complete", message = "Monthly rainfall totals written")

baseline_query <- sprintf(
  paste(
    "SELECT pixel, month,",
    "AVG(monthly_rain) AS mean_monthly_rain,",
    "MEDIAN(monthly_rain) AS median_monthly_rain,",
    "quantile_cont(monthly_rain, 0.1) AS p10_monthly_rain,",
    "quantile_cont(monthly_rain, 0.9) AS p90_monthly_rain,",
    "COUNT(*) AS years_observed",
    "FROM read_parquet('%s')",
    "GROUP BY pixel, month",
    sep = "\n"
  ),
  monthly_sql
)
DBI::dbExecute(
  connection,
  sprintf(
    "COPY (%s) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD, OVERWRITE_OR_IGNORE TRUE)",
    baseline_query, baseline_sql
  )
)
log_event(run, event = "baseline_complete", message = "Baseline monthly climatology written")

annual_query <- sprintf(
  paste(
    "WITH monthly AS (",
    "  SELECT *, SUM(monthly_rain) OVER (PARTITION BY pixel, year) AS annual_rain",
    "  FROM read_parquet('%s')",
    "),",
    "components AS (",
    "  SELECT *,",
    "    2 * pi() * (month - 0.5) / 12 AS theta",
    "  FROM monthly",
    ")",
    "SELECT pixel, year,",
    "  COUNT(*) AS months_observed,",
    "  MAX(annual_rain) AS annual_rain,",
    "  SUM(abs(monthly_rain - annual_rain / 12)) / MAX(annual_rain) AS rainfall_si,",
    "  sqrt(pow(SUM(monthly_rain * cos(theta)), 2) +",
    "       pow(SUM(monthly_rain * sin(theta)), 2)) / MAX(annual_rain) AS rain_h1,",
    "  sqrt(pow(SUM(monthly_rain * cos(2 * theta)), 2) +",
    "       pow(SUM(monthly_rain * sin(2 * theta)), 2)) / MAX(annual_rain) AS rain_h2,",
    "  arg_max(month, monthly_rain) AS wettest_month",
    "FROM components",
    "GROUP BY pixel, year",
    sep = "\n"
  ),
  monthly_sql
)

annual <- as.data.table(DBI::dbGetQuery(connection, annual_query))
annual[, `:=`(
  baseline_median_rain = median(annual_rain, na.rm = TRUE),
  baseline_mad_rain = stats::mad(annual_rain, center = median(annual_rain, na.rm = TRUE), constant = 1, na.rm = TRUE)
), by = pixel]
annual[, wet_anomaly := fifelse(
  is.finite(baseline_mad_rain) & baseline_mad_rain > 0,
  (annual_rain - baseline_median_rain) / baseline_mad_rain,
  NA_real_
)]

pheno <- as.data.table(read_parquet(pheno_file, col_select = c("pixel", "admin1_name")))
pixel_admin <- unique(pheno)
annual <- pixel_admin[annual, on = .(pixel)]
setcolorder(annual, c("pixel", "admin1_name", "year"))
write_parquet(annual, annual_file)

monthly_rows <- DBI::dbGetQuery(
  connection,
  sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s')", monthly_sql)
)$n[[1]]
baseline_rows <- DBI::dbGetQuery(
  connection,
  sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s')", baseline_sql)
)$n[[1]]

DBI::dbDisconnect(connection, shutdown = TRUE)
connection <- NULL
failed <- FALSE
finish_run(
  run,
  "success",
  list(
    monthly_output = monthly_file,
    annual_output = annual_file,
    baseline_output = baseline_file,
    monthly_rows = monthly_rows,
    annual_rows = nrow(annual),
    baseline_rows = baseline_rows
  )
)

cat("Wrote", monthly_file, "\n")
cat("Wrote", annual_file, "\n")
cat("Wrote", baseline_file, "\n")
