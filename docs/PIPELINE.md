# Kenya recovery pipeline

All commands run from repository root. Large inputs and outputs live below `ANALOGUE_DATA_ROOT` (default `/Volumes/clim_dat`), not in Git.

## Data flow

```text
GLASS GeoTIFFs ──> annual NDVI metrics ──────────────┐
                                                     ├─> stable detectability inputs ─> pathway classes
CHIRPS daily ─────> rainfall metrics ─> windows ─────┤                 │
                                      │              │                 ├─> calibration report/panels
legacy phenology ─────────────────────┴─> stable events               │
                                                                       │
CEEPA survey ─────> planting/harvest observations ─────────────────────┴─> independent validation
```

## Run order

```bash
Rscript scripts/preflight.R
Rscript scripts/build_kenya_signal_metrics.R
Rscript scripts/build_kenya_rainfall_metrics.R
Rscript scripts/build_kenya_season_windows.R
Rscript scripts/build_kenya_stable_phenology.R
Rscript scripts/assemble_kenya_detectability_inputs.R
Rscript scripts/classify_kenya_detectability.R
Rscript scripts/build_kenya_calibration_report.R
Rscript scripts/render_kenya_review_panels.R
```

Independent validation:

```bash
Rscript scripts/download_ceepa_validation.R
Rscript scripts/build_kenya_ceepa_validation.R
Rscript scripts/score_kenya_ceepa_validation.R
Rscript scripts/analyze_kenya_validation_regimes.R
```

Run checks after code or configuration changes:

```bash
Rscript tests/run_tests.R
Rscript scripts/preflight.R
```

## Current safeguards

- GLASS `AYYYYDDD` filenames use strict calendar parsing.
- Rainfall windows build directly from source metadata; no circular dependency on classifier outputs.
- Raw fitted season numbers are retained only as provenance. Stable `season_id` drives analysis.
- Missing/weak humid signals can remain `not_identifiable`; dates are never forced.
- Detectability thresholds remain `provisional` in `config/project.json` until human panel review and further validation are complete.
- Runtime logs write to ignored `logs/*.jsonl`; scientific decisions and material results enter tracked documentation.
