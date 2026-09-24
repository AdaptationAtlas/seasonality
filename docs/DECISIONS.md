# Decision log

## D001 — Kenya-first scope

- Date: 2026-09-24
- Status: accepted
- Decision: Refactor and validate Kenya workflow before returning to Africa-wide processing.
- Reason: Kenya contains arid, bimodal, highland, humid, and evergreen regimes needed to resolve known method failure while keeping iteration tractable.

## D002 — No forced phenology

- Date: 2026-09-24
- Status: accepted
- Decision: `not_identifiable` is valid output. Do not force onset or season count where NDVI and rainfall lack sufficient seasonal structure.
- Reason: Evergreen/humid areas can remain green and wet; invented dates create false precision.

## D003 — Source-aware fallback

- Date: 2026-09-24
- Status: proposed; requires validation
- Decision: Use three pathways: NDVI-derived when identifiable; rainfall-derived proxy when NDVI is weak but rainfall seasonal; no date when both are weak/aseasonal. Mark wet-year merged seasons separately.
- Reason: Current Kenya data show lower complete NDVI phenology coverage in humid pixels (53.4%) than arid pixels (75.8%). Fit scores alone remain high and therefore do not measure detectability.

## D004 — Stable baseline season windows

- Date: 2026-09-24
- Status: proposed; requires validation
- Decision: Learn season windows from robust baseline years and assign each year within those windows. Avoid re-clustering national seasons independently during each summary run.
- Reason: Wet years can bridge bimodal seasons; free clustering permits label swaps and merged-season artefacts.

## D005 — Data outside Git; logs split by purpose

- Date: 2026-09-24
- Status: accepted
- Decision: Keep raster/Parquet source data and generated outputs outside Git. Track concise scientific/project history in `docs/WORKLOG.md`; write machine runtime events as ignored JSONL files.
- Reason: Data are large; both human provenance and machine-debuggable logs are needed.

## D006 — Detectability requires multiple signals

- Date: 2026-09-24
- Status: accepted
- Decision: Never classify detectability from NDVI amplitude or fit score alone. Combine annual amplitude, timing concentration, event-year coverage, rainfall seasonality, wet anomaly, and land-cover context.
- Reason: Kenya median NDVI amplitude is low in both humid (0.160) and arid (0.127) classes for different ecological reasons; fit scores remain deceptively high.

## D007 — Preserve full domain; flag agricultural relevance

- Date: 2026-09-24
- Status: accepted
- Decision: Keep all Kenya pixels in scientific derivatives, but include land cover, mapped crop activity, and crop-presence flags. Decision summaries can filter or weight agricultural pixels explicitly.
- Reason: Hard masking would hide ecological transitions and uncertain crop-map omissions; unflagged all-land summaries would let forest, bare ground, and shrubland dominate agricultural interpretation.

## D008 — Rainfall windows are candidates until calibrated

- Date: 2026-09-24
- Status: accepted
- Decision: Derive fixed circular windows from baseline rainfall peaks and troughs, but accept bimodality only after peak/valley strength, NDVI timing, crop relevance, and validation agree.
- Reason: Two mathematical rainfall peaks occur for 18,253 phenology pixels, but median dry-valley strength differs sharply: about 0.81 in arid crop pixels versus 0.28 in humid crop pixels. Counting peaks alone would overstate meaningful bimodality.

## D009 — Viewer follows stable data contracts

- Date: 2026-09-24
- Status: accepted
- Decision: Build interactive viewer only after annual indicators, confidence fields, and admin aggregation schemas stabilize. Viewer must use derived products, not scan raw climate archives.
- Reason: Separating computation from presentation keeps viewer responsive, reproducible, and cheap to deploy while preventing UI choices from shaping scientific methods.
