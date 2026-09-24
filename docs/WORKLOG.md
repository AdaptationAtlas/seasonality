# Work log

## 2026-09-24 — Recovery audit and project foundation

- Confirmed Git remote `origin` points to `AdaptationAtlas/seasonality`; authenticated GitHub account `peetmate`; local `main` matches `origin/main` before new work.
- Preserved existing uncommitted edits in dataset download, admin1 summary, extracted helper files, archive move, and generated output.
- Parsed all active R scripts successfully.
- Audited mounted climate-data tree and Kenya files.
- Found invalid 408-byte admin2 boundary file caused by old S3 key typo (`level=adm1` for an admin2 filename).
- Quantified Kenya humid-region detection gap: complete Greenup–Senescence pairs decline from 75.8% in arid pixels to 53.4% in humid pixels despite high median fit scores.
- Added AI working rules, config, roadmap, decisions, data inventory, structured JSONL runtime logger, preflight, and initial circular/logging tests.
- Changed setup script to respect pre-set `ANALOGUE_DATA_ROOT`, avoid overwriting setup metadata, and suppress full tree output unless requested.
- Added annual raw-NDVI signal metric builder and explicit source-path classifier. Thresholds remain calibration inputs, not assumed truths.
- Built 463,375 annual raw-NDVI metric rows for 18,535 Kenya pixels covering 2000–2024. Found low amplitude in both humid evergreen and arid sparse-vegetation settings, confirming need for multi-signal classification.
- Replaced corrupt admin2 placeholder with validated 17 MB Parquet from current S3 key; retained original XML error file for provenance.
- Added out-of-memory-safe DuckDB aggregation for Kenya CHIRPS: monthly totals, annual rainfall seasonality/harmonics, robust wet anomalies, and baseline monthly climatology.
- Added calibration-table assembly joining NDVI, rainfall, event timing, aridity, elevation, and land cover by shared raster pixel ID.
- Added total MapSPAM crop activity and crop-presence flag to distinguish agricultural pixels from evergreen forest, bare ground, and other weak-signal regimes.
- Added circular rainfall peak detection and candidate fixed windows. Window diagnostics retain uncertainty; no bimodal threshold is forced before calibration.
- Built fixed-window candidates for 18,253 pixels. Dominant crop-pixel rainfall peaks are April/November. Dry-valley strength is much lower in humid crop pixels (~0.28 median) than arid crop pixels (~0.81), confirming peak count alone cannot define seasons.

Next: build/test Kenya detectability classifier and stable baseline season windows on representative regions.
