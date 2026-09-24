# Data inventory

Audit date: 2026-09-24. Root: `/Volumes/clim_dat` (override with `ANALOGUE_DATA_ROOT`).

## Kenya-ready inputs

| Dataset | Path | Observed state |
|---|---|---|
| Kenya phenology + rainfall | `climate_derived/glass_phenology/countries/KEN_seasonal-phenology_plus-rain.parquet` | 59 MB; 672,900 rows; 18,535 pixels; 48 admin1; 2000–2024 |
| Kenya phenology | `climate_derived/glass_phenology/countries/KEN_seasonal-phenology.parquet` | 44 MB |
| Pixel index | `climate_derived/glass_phenology/pixel_index.parquet` | 12 MB |
| Pixel aridity/elevation | `climate_derived/glass_phenology/pixel_index_aridity_elevation.parquet` | 36 MB |
| Africa admin1 boundaries | `static_raw/boundaries/atlas_gaul24_a1_africa_simple-highres.parquet` | Valid; 719 rows |
| Africa admin2 boundaries | `static_raw/boundaries/atlas_gaul24_a2_africa_simple-highres.parquet` | Valid; 17 MB; 4,596 Africa rows; 290 Kenya rows |
| CEEPA agricultural survey | `static_raw/validation/ceepa/CEEPA.zip` | Valid Figshare archive; MD5 verified; 9,597 Africa households, including 816 Kenya households in 44 historical districts |

## Main upstream stores

- GLASS NDVI source and converted GeoTIFFs: `climate_raw/glass_ndvi*`
- GLASS fitted phenology: `climate_derived/glass_phenology/glass_nvdi_phenofit.parquet` (1.7 GB)
- CHIRPS v3 daily COGs: `climate_raw/chirps/chirps_v3_cog`
- Hobbins reference ET: `climate_raw/hobbins_ref_et/1983` through `2024`
- Aridity: `static_raw/aridity`
- SRTM elevation: `static_raw/strm`
- Existing run outputs: `output/2026-03-27` and `output/2026-03-28`

## Kenya recovery derivatives

- `KEN_annual_signal_metrics.parquet`: annual NDVI amplitude, coverage, and event counts.
- `KEN_baseline_season_signal_metrics.parquet`: baseline event coverage and circular timing concentration by candidate season.
- `KEN_monthly_rainfall.parquet`: CHIRPS monthly totals for 2000–2024.
- `KEN_annual_rainfall_metrics.parquet`: annual rainfall total, Walsh–Lawler seasonality, first/second harmonics, and robust wet anomaly.
- `KEN_baseline_monthly_rainfall.parquet`: monthly rainfall climatology.
- `KEN_detectability_inputs.parquet`: joined pixel-season-year calibration features including aridity, elevation, land cover, and total mapped crop activity.
- `KEN_candidate_season_windows.parquet`: two candidate rainfall peaks and fixed circular month windows per pixel, with valley-strength diagnostic.
- `KEN_stable_phenology_events.parquet`: fitted events relabelled into fixed baseline rainfall windows; includes explicit missing seasons and duplicate-candidate counts.
- `independent_validation/KEN_ceepa_crop_observations.parquet`: decoded Kenya household crop planting/harvest observations with precision and eligibility flags.
- `independent_validation/KEN_ceepa_crop_calendar_summary.parquet`: crop-calendar summaries by mapped current county, historical district, and named survey season.
- `independent_validation/KEN_ceepa_remote_validation.*`: county-season comparison of 2003 CEEPA planting/harvest dates with stable GLASS greenup/senescence events.

## Kenya audit signal

Complete Greenup–Senescence pair rate by aridity class:

| Class | Rows | Pixels | Complete pair rate |
|---|---:|---:|---:|
| Arid | 264,194 | 7,117 | 75.8% |
| Semi-arid | 270,553 | 7,203 | 65.7% |
| Sub-humid | 56,682 | 1,713 | 60.7% |
| Humid | 81,471 | 2,502 | 53.4% |

Median R² remains high (0.933 in humid class). Conclusion: fitted-curve goodness alone cannot establish seasonal detectability. Need amplitude/timing-concentration and rainfall-seasonality diagnostics.

Raw annual NDVI audit adds second warning: correctly calendar-aligned median annual amplitude is 0.160 in humid pixels and 0.136 in arid pixels, versus 0.227–0.230 in semi-arid/sub-humid pixels. Low amplitude therefore occurs both in evergreen humid systems and sparsely vegetated arid systems. Classification must combine amplitude, timing concentration, event coverage, rainfall seasonality, and land cover; amplitude alone is unsafe.

Initial independent comparison covers 61 county-seasons with at least 10 CEEPA observations and 10 quality-screened mapped-crop pixels. Median GLASS greenup occurs 10 days after reported planting; median absolute timing difference is 15.5 days. Remote quality-event coverage is 0.59 median. Larger western/humid discrepancies require targeted review rather than global threshold tuning. Planting-to-greenup lag and planting/harvest versus greenup/senescence definitions mean these are diagnostic differences, not interchangeable dates.

## Known integrity concerns

- Corrupt admin2 response preserved as `atlas_gaul24_a2_africa_simple-highres.invalid-20260924.xml`; it records old `NoSuchKey` error.
- Boundary download code points to newer S3 names while local filenames retain older names. New admin2 object was downloaded from `boundaries/atlas-region_admin2_simplified.parquet` and validated before replacement.
- `R/1-2_process_nvdi.R` sources remote R code at runtime, harming reproducibility.
- Large pixel processing writes one `.RData` file per pixel and notes memory leak.
- `nvdi` misspelling is embedded in existing path/key names. Preserve compatibility until planned migration.
