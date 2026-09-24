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

## Main upstream stores

- GLASS NDVI source and converted GeoTIFFs: `climate_raw/glass_ndvi*`
- GLASS fitted phenology: `climate_derived/glass_phenology/glass_nvdi_phenofit.parquet` (1.7 GB)
- CHIRPS v3 daily COGs: `climate_raw/chirps/chirps_v3_cog`
- Hobbins reference ET: `climate_raw/hobbins_ref_et/1983` through `2024`
- Aridity: `static_raw/aridity`
- SRTM elevation: `static_raw/strm`
- Existing run outputs: `output/2026-03-27` and `output/2026-03-28`

## Kenya audit signal

Complete Greenup–Senescence pair rate by aridity class:

| Class | Rows | Pixels | Complete pair rate |
|---|---:|---:|---:|
| Arid | 264,194 | 7,117 | 75.8% |
| Semi-arid | 270,553 | 7,203 | 65.7% |
| Sub-humid | 56,682 | 1,713 | 60.7% |
| Humid | 81,471 | 2,502 | 53.4% |

Median R² remains high (0.933 in humid class). Conclusion: fitted-curve goodness alone cannot establish seasonal detectability. Need amplitude/timing-concentration and rainfall-seasonality diagnostics.

Raw annual NDVI audit adds second warning: median annual amplitude is 0.160 in humid pixels and 0.127 in arid pixels, versus 0.220–0.229 in semi-arid/sub-humid pixels. Low amplitude therefore occurs both in evergreen humid systems and sparsely vegetated arid systems. Classification must combine amplitude, timing concentration, event coverage, rainfall seasonality, and land cover; amplitude alone is unsafe.

## Known integrity concerns

- Corrupt admin2 response preserved as `atlas_gaul24_a2_africa_simple-highres.invalid-20260924.xml`; it records old `NoSuchKey` error.
- Boundary download code points to newer S3 names while local filenames retain older names. New admin2 object was downloaded from `boundaries/atlas-region_admin2_simplified.parquet` and validated before replacement.
- `R/1-2_process_nvdi.R` sources remote R code at runtime, harming reproducibility.
- Large pixel processing writes one `.RData` file per pixel and notes memory leak.
- `nvdi` misspelling is embedded in existing path/key names. Preserve compatibility until planned migration.
