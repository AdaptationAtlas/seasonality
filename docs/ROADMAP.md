# Kenya seasonality roadmap

Status values: `TODO`, `IN PROGRESS`, `BLOCKED`, `DONE`.

## Outcome

Decision-ready Kenya admin1/admin2 time series describing onset, end, season length, baseline variability, trends, minima, maxima, coverage, method source, and uncertainty.

## Workstreams

| ID | Status | Deliverable | Acceptance check |
|---|---|---|---|
| P0.1 | DONE | Repository and data audit | Code, Git, Kenya inputs, and known defects documented |
| P0.2 | DONE | AI project scaffolding | `AGENTS.md`, config, worklog, decisions, preflight, tests, JSONL logging |
| P0.3 | TODO | Dependency lock | `renv.lock` restores clean environment without `pacman` auto-install |
| P1.1 | IN PROGRESS | Kenya detectability classifier | NDVI/rainfall builders, combined calibration table, and pathway function added; thresholds require calibration/validation |
| P1.2 | IN PROGRESS | Stable Kenya season windows | Circular rainfall windows and stable event reassignment built; bimodal acceptance threshold still requires calibration |
| P1.3 | TODO | Phenology estimator | Source-aware dates plus confidence; no forced dates for aseasonal pixels |
| P1.4 | IN PROGRESS | Independent validation | CEEPA 2003 comparison built for 61 well-supported county-seasons; next review humid/western outliers and score by ecological regime |
| P2.1 | TODO | Historical indicators | Annual onset/end/length plus baseline quantiles, variability, min/max, and robust trends |
| P2.2 | DONE | Admin2 boundaries | Replaced corrupt XML response with validated Parquet: 4,596 Africa rows, 290 Kenya rows |
| P2.3 | TODO | Admin1/admin2 aggregation | Area-aware summaries, coverage, confidence, and method composition |
| P3.1 | TODO | Reproducible pipeline | Staged `{targets}` pipeline with resumable Kenya-only execution |
| P3.2 | TODO | Decision products | Parquet time series, COG rasters, QA report, metadata/data dictionary |
| P4.1 | TODO | Interactive data viewer | Map with year/season/metric controls, source/confidence filters, admin1/admin2 drill-down, pixel time series, and downloads |
| P4.2 | TODO | Viewer deployment | Publish documented, versioned viewer using derived products only; include data/version timestamp and methodology links |

## Immediate sequence

1. Extract small representative Kenya sample: humid west/highlands, bimodal southeast/central, arid north/east.
2. Quantify NDVI amplitude, circular timing concentration, rainfall seasonality, missing events, and season-count stability.
3. Implement detectability classes and stable season windows as tested pure functions.
4. Compare NDVI and rainfall pathways against known Kenya calendars.
5. Produce pilot admin1 outputs; add admin2 after boundary repair.
6. Build interactive viewer after indicators and aggregation contracts stabilize.

## Definition of done

- Fresh environment can restore dependencies and run pipeline from documented commands.
- Every run records config, input fingerprints, warnings, outputs, duration, and status.
- Every indicator includes observation count, coverage, method source, and confidence.
- Circular date statistics and trends pass year-boundary tests.
- Kenya validation report states where onset is reliable, proxy-derived, or not identifiable.
- Viewer exposes spatial layers, historical time series, uncertainty, method source, and downloadable admin summaries.
