# Kenya seasonality

Reproducible estimation of growing-season onset, end, and length from GLASS NDVI and CHIRPS rainfall. Current development focuses on Kenya and historical baseline 2000–2024.

Project handles four states explicitly:

1. NDVI-seasonal: vegetation phenology supplies dates.
2. Weak NDVI but rainfall-seasonal: rainfall may supply proxy dates after validation.
3. Wet-year merged season: baseline windows retain stable season identity.
4. Evergreen/aseasonal: onset is reported as not identifiable.

## Setup

Source data remain outside repository. Default root is `/Volumes/clim_dat`; override when needed:

```sh
export ANALOGUE_DATA_ROOT=/path/to/climate-data
Rscript scripts/preflight.R
Rscript tests/run_tests.R
```

Runtime logs are written as JSON Lines under `logs/` and ignored by Git.

## Project records

- [Roadmap](docs/ROADMAP.md)
- [Scientific and technical decisions](docs/DECISIONS.md)
- [Data inventory](docs/DATA_INVENTORY.md)
- [Work log](docs/WORKLOG.md)
- [AI coding instructions](AGENTS.md)

## Current pipeline

| Stage | Script | Purpose |
|---|---|---|
| Setup | `R/0-1_setup_folders.R` | Resolve external data root and expected folders |
| Download | `R/0-2_download_datasets.R` | Acquire climate/static inputs |
| Rainfall | `R/1-1_process_chirps.R` | Prepare CHIRPS time series |
| Phenology | `R/1-2_process_nvdi.R` | Fit GLASS NDVI curves and extract candidate events |
| Summaries | `R/2-1_admin1_pheno_summaries.r` | Harmonize seasons, apply QC, and create outputs |

Current scripts remain legacy orchestration. Roadmap refactors them into tested, resumable Kenya stages before scientific production use.
