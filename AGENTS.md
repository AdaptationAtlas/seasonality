# Seasonality project instructions

## Mission

Build reproducible Kenya seasonal-climate indicators from GLASS NDVI and CHIRPS rainfall. Primary outputs: onset date, end date, season length, historical variability, trends, extrema, coverage, and confidence at pixel, admin1, and admin2 levels.

## Current scope

- Country: Kenya (`KEN`). Do not expand Africa-wide workflow until Kenya method validates.
- Available phenology baseline: 2000–2024.
- Oceanic-driver analysis is future work. Preserve yearly time series and metadata needed for later linkage, but do not model drivers now.
- Source data live under `ANALOGUE_DATA_ROOT` (default `/Volumes/clim_dat`) and must not enter Git.

## Working rules

1. Read `docs/ROADMAP.md`, `docs/DECISIONS.md`, and latest `docs/WORKLOG.md` entry before material changes.
2. Preserve user work. Check `git status` and relevant diffs before edits.
3. Keep configuration outside analysis logic. Put project defaults in `config/project.json`.
4. Put reusable logic in `R/functions/`; orchestration scripts should stay thin.
5. Every executable stage must emit JSONL events through `R/functions/project_logging.R` and record inputs, parameters, outputs, duration, warnings, and status.
6. Never force one or two seasons where signal is not identifiable. Emit explicit detectability class, method source, coverage, and confidence.
7. Treat day-of-year as circular. Tests must cover year-boundary cases.
8. Prefer Parquet for tables and COG GeoTIFF for rasters. Include stable keys: `pixel`, `year`, `season_id`, `admin1_code`, `admin2_code` where applicable.
9. Run `Rscript scripts/preflight.R` and `Rscript tests/run_tests.R` before handing off code changes.
10. Append concise dated entry to `docs/WORKLOG.md` after each material work session. Record scientific decisions in `docs/DECISIONS.md`.

## Scientific guardrails

- Separate phenology detectability from model goodness-of-fit.
- Use NDVI-derived dates only where seasonal vegetation response is identifiable.
- In weak-NDVI but rainfall-seasonal areas, permit rainfall-derived dates with distinct `method_source` and validation status.
- In evergreen/aseasonal areas, return `not_identifiable`; do not manufacture onset dates.
- Derive stable season windows from baseline climatology before classifying individual wet/merged years.
- Report uncertainty and observation coverage beside every summary and trend.
- Validate against independent Kenya crop calendars or planting-date observations before decision-use claims.

## Common commands

```sh
Rscript scripts/preflight.R
Rscript tests/run_tests.R
```

## Git workflow

- Remote: `https://github.com/AdaptationAtlas/seasonality.git`.
- Use focused branches and commits. Never mix recovered legacy edits with new scaffolding without explicit review.
- Do not commit source data, generated rasters, runtime logs, `.Rhistory`, `.Rproj.user`, or `.DS_Store`.
