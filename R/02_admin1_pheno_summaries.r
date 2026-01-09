pacman::p_load(circular, changepoint, data.table, arrow, lubridate, ggplot2, terra)
source("R/00_setup_folders.R")

# Create Functions ####
## --- find automatic split point between the two clusters ---
assign_season <- function(doy_vec, gap_factor = 2) {
  x <- sort(doy_vec)
  gaps <- diff(x)

  if (length(gaps) ==  0) {
    # only one observation → certainly one season
    return(rep(1L, length(doy_vec)))
  }

  max_gap <- max(gaps)
  med_gap <- median(gaps)

  # decide if distribution is bimodal
  is_bimodal <- max_gap > gap_factor * med_gap

  if (!is_bimodal) {
    # unimodal → assign season 1 only
    return(rep(1L, length(doy_vec)))
  }

  # bimodal → find threshold
  idx <- which.max(gaps)
  split <- floor((x[idx] + x[idx + 1]) / 2)

  as.integer(doy_vec > split) + 1L
}

# x = vector of Date or POSIXct
split_seasons_cpt <- function(x, max_seasons = 3L) {
  doy_raw <- yday(x)

  # trivial cases
  n <- length(doy_raw)
  if (n <=  2L || max_seasons ==  1L) {
    return(rep(1L, n))
  }

  # order by DOY (seasons must be contiguous in DOY space)
  o        <- order(doy_raw)
  doy_sort <- doy_raw[o]

  # at most max_seasons - 1 change points
  res <- cpt.mean(
    doy_sort,
    method  = "BinSeg",      # binary segmentation, fast
    Q       = max_seasons - 1,
    penalty = "BIC"          # lets BIC choose how many cps are justified
  )

  cps <- cpts(res)           # indices of change points in the *sorted* vector

  # assign segment IDs to sorted indices
  seg_sorted <- findInterval(seq_len(n), vec = cps) + 1L
  # map back to original order
  seg <- integer(n)
  seg[o] <- seg_sorted

  seg
}

doy_to_circular <- function(doy) {
  circular(2 * pi * (doy / 365), units = "radians", modulo = "2pi")
}

quantile_circular_safe <- function(
    doy,
    probs = c(0.05, 0.5, 0.95),
    n_days = 365) {
  r <- range(doy, na.rm = TRUE)
  # if the cluster is very wide, assume it crosses year boundary and rotate by 180 d # nolint: line_length_linter.
  shift <- if ((r[2] - r[1]) > (n_days / 2)) n_days / 2 else 0

  doy_rot <- ((doy - shift - 1L) %% n_days) + 1L
  qs_rot  <- quantile(doy_rot, probs, na.rm = TRUE)

  # rotate back
  ((qs_rot + shift - 1L) %% n_days) + 1L
}

season_length_doy <- function(sos, eos, n_days = 365) {
  (eos - sos + n_days) %% n_days
}

# Choose & load data ####


# Load pixel index #####
coords_index <- read_parquet(
  file.path(dirs$nvdi_phenology, "pixel_index.parquet")
)

pixel_map <- rast(
x = coords_index,
type = "xyz", # indicates x, y, value columns
crs = "EPSG:4326" # adjust if needed
)

# Assign values using the pixel column
values(pixel_map) <- coords_index$pixel

# Load / compute DEM + Aridity Index by pixel #####
# Set output file for DEM and Aridity Index by pixel
dem_ai_file<-file.path(dirs$nvdi_phenology, "pixel_index_aridity_elevation.parquet")

if(file.exists(dem_ai_file)){
  px_env <- setDT(read_parquet(dem_ai_file))
} else {
  # Load DEM
  dem <- rast(file.path(dirs$srtm, "DEM_SRTM_Africa.tif"))

  # Load Aridity Index
  ai_r <- rast(list.files(
    dirs$aridity, 
    pattern = "\\.tif$",
    full.names = TRUE
  ))

  # Align to pixel_map
  dem  <- resample(dem,  pixel_map, method = "bilinear")
  ai_r <- resample(ai_r, pixel_map, method = "bilinear")

  # Zonal mean by pixel ID
  dem_px <- terra::zonal(dem,  pixel_map, fun = "mean", na.rm = TRUE)
  ai_px  <- terra::zonal(ai_r, pixel_map, fun = "mean", na.rm = TRUE)

  setDT(dem_px); setDT(ai_px)
  setnames(dem_px, c("zone", "elevation"))
  setnames(ai_px,  c("zone", "aridity"))

  px_env <- merge(dem_px, ai_px, by = "zone", all = TRUE)
  setnames(px_env, "zone", "pixel")

  px_env[elevation>0 & aridity]

    # Divide aridity by 10000
  px_env[, aridity := aridity / 10000]

  # Aridity bins (UNEP-style defaults; adjust if you want)
  px_env[, aridity_bin := cut(
    aridity,
    breaks = c(-Inf, 0.05, 0.20, 0.50, 0.65, Inf),
    labels = c("hyper-arid", "arid", "semi-arid", "sub-humid", "humid")
  )]

  # Elevation bins (from earlier suggestion)
  px_env[, elev_bin := cut(
    elevation,
    breaks = c(-Inf, 800, 1500, 2200, Inf),
    labels = c("<800", "800-1500", "1500-2200", ">2200")
  )]
  write_parquet(px_env, dem_ai_file)
}

# Load Phenology data #####
# Folder for country phenology data
pheno_dat <- file.path(dirs$nvdi_phenology, "countries")

# Choose ISO3 country code
files<-list.files(pheno_dat, "plus-rain", full.names = TRUE)
iso3_choice <- "KEN"
file_choice <- grep(iso3_choice, files, value = TRUE)

# Load data
dat <- setDT(read_parquet(file_choice))

# Merge dem + aridity + coords #####
dat <- px_env[dat, on = .(pixel)]

# Create combined zone ID
dat[, zone_id := paste(admin1_name, aridity_bin, elev_bin, sep = "|")]

# Prepare Phenology Data #####
# Split year and season from flag
parts <- tstrsplit(dat$flag, "_")
set(dat, j = "year",   value = parts[[1]])
set(dat, j = "season", value = parts[[2]])

# Count season 2 occurrences per pixel
dat[, season2_n :=  sum(season == 2), by = pixel]

dat[, Greenup   := as.IDate(Greenup)]
dat[, Senescence := as.IDate(Senescence)]
dat[, DER.sos   := as.IDate(DER.sos)]
dat[, DER.eos   := as.IDate(DER.eos)]
dat[, TRS2.sos  := as.IDate(TRS2.sos)]
dat[, TRS2.eos  := as.IDate(TRS2.eos)]

# Step 1: choose best available season source per row ######

# Priority: Greenup > DER > TRS2 (requires both sos/eos to be present)
dat[, season_source := fifelse(!is.na(Greenup) & !is.na(Senescence), "greenup",
                        fifelse(!is.na(DER.sos) & !is.na(DER.eos), "der",
                        fifelse(!is.na(TRS2.sos) & !is.na(TRS2.eos), "trs2", NA_character_)))]

# unify timing + rainfall metrics based on season_source
dat[, sos := fifelse(season_source == "greenup", Greenup,
              fifelse(season_source == "der", DER.sos,
              fifelse(season_source == "trs2", TRS2.sos, as.IDate(NA))))]

dat[, eos := fifelse(season_source == "greenup", Senescence,
              fifelse(season_source == "der", DER.eos,
              fifelse(season_source == "trs2", TRS2.eos, as.IDate(NA))))]

dat[, rain_total := fifelse(season_source == "greenup", rain_greenup_scenescence,
                     fifelse(season_source == "der", rain_der,
                     fifelse(season_source == "trs2", rain_trs2, NA_real_)))]

dat[, rain_p30 := fifelse(season_source == "greenup", rain_greenup_p30,
                   fifelse(season_source == "der", rain_der_p30,
                   fifelse(season_source == "trs2", rain_trs2_p30, NA_real_)))]

dat[, cdd_p45 := fifelse(season_source == "greenup", cdd_greenup_p45,
                  fifelse(season_source == "der", cdd_der_p45,
                  fifelse(season_source == "trs2", cdd_trs2_p45, NA_real_)))]

# Set season start / end dates
dat[, doy :=  yday(sos)]
dat[, doy_eos :=  yday(eos)]
dat[, slen :=  as.numeric(eos - sos)]

# Harmonize season numbering within each admin1 region
dat[!is.na(sos),
    season_harmonized :=  split_seasons_cpt(x = sos, max_seasons = 2L),
    by = admin1_name]

# Step 2: Gates / QC flags (zone-based, aridity-aware) ######

#  Parameters

# Quantile gates (zone-specific thresholds):
#   p30_q_fail / p30_q_weak  : lower-tail cutoffs for early rain (p30).
#                              p30 <= q_fail => strong failure signal
#                              p30 <= q_weak => weak signal
#                              Higher values = stricter (more fails).
#   rtot_q_fail / rtot_q_weak: lower-tail cutoffs for total seasonal rain.
#                              Higher values = stricter (more fails).
#   cdd_q_fail / cdd_q_weak  : upper-tail cutoffs for dry-spell length (CDD).
#                              cdd >= q_fail => strong failure signal
#                              cdd >= q_weak => weak signal
#                              Lower values = stricter (more fails).

p30_q_fail  <- 0.10
p30_q_weak  <- 0.20
cdd_q_fail  <- 0.90
cdd_q_weak  <- 0.80
rtot_q_fail <- 0.10
rtot_q_weak <- 0.20

# Hard limits (absolute thresholds, applied only in arid/semi/hyper-arid bins):
  #   p30_hard_*  : minimum early rain. Raise -> stricter; lower -> looser.
  #   rtot_hard_* : minimum total rain. Raise -> stricter; lower -> looser.
  #   cdd_hard_*  : maximum dry-spell length. Lower -> stricter; higher -> looser.
  #   In humid/sub-humid bins, hard-fail is disabled to avoid over-penalizing.

p30_hard_fail_arid <- 5
p30_hard_weak_arid <- 15

cdd_hard_fail_arid <- 35
cdd_hard_weak_arid <- 25

rtot_hard_fail_arid <- 50
rtot_hard_weak_arid <- 120

# NDVI fit quality (protect Greenup-based seasons):
  #   ndvi_nse_min / ndvi_r2_min: if below these, Greenup seasons can be failed
  #   by rainfall gates. Higher thresholds = stricter NDVI quality requirement.
ndvi_nse_min <- 0.60
ndvi_r2_min  <- 0.60

# min_n:
  #   Minimum sample size per zone (zone_id + season_harmonized) to compute
  #   stable quantiles. If n < min_n, thresholds are set to NA (you may want
  #   a fallback to broader zones).

min_n <- 30  # minimum rows per zone for stable quantiles

#  1) Zone-specific thresholds (zone_id + season_harmonized) #######
zone_thr <- dat[
  !is.na(zone_id) & !is.na(season_harmonized),
  .(
    n_zone = .N,
    p30_fail_thr  = as.numeric(quantile(rain_p30,   p30_q_fail,  na.rm = TRUE)),
    p30_weak_thr  = as.numeric(quantile(rain_p30,   p30_q_weak,  na.rm = TRUE)),
    cdd_fail_thr  = as.numeric(quantile(cdd_p45,    cdd_q_fail,  na.rm = TRUE)),
    cdd_weak_thr  = as.numeric(quantile(cdd_p45,    cdd_q_weak,  na.rm = TRUE)),
    rtot_fail_thr = as.numeric(quantile(rain_total, rtot_q_fail, na.rm = TRUE)),
    rtot_weak_thr = as.numeric(quantile(rain_total, rtot_q_weak, na.rm = TRUE))
  ),
  by = .(zone_id, season_harmonized)
]

# If small sample sizes, set thresholds to NA (optional fallback)
zone_thr[n_zone < min_n,
  c("p30_fail_thr","p30_weak_thr","cdd_fail_thr","cdd_weak_thr","rtot_fail_thr","rtot_weak_thr") := NA_real_
]

# Join thresholds back
dat[zone_thr, on = .(zone_id, season_harmonized), `:=`(
  p30_fail_thr  = i.p30_fail_thr,
  p30_weak_thr  = i.p30_weak_thr,
  cdd_fail_thr  = i.cdd_fail_thr,
  cdd_weak_thr  = i.cdd_weak_thr,
  rtot_fail_thr = i.rtot_fail_thr,
  rtot_weak_thr = i.rtot_weak_thr
)]

# 2) NDVI reliability indicator #######
dat[, ndvi_ok := {
  has_fit <- !is.na(NSE) & !is.na(R2)
  fifelse(has_fit, (NSE >= ndvi_nse_min & R2 >= ndvi_r2_min), TRUE)
}]

# 3) Quantile-based flags #######
dat[, `:=`(
  q_p30_fail  = !is.na(rain_p30)   & !is.na(p30_fail_thr)  & (rain_p30   <= p30_fail_thr),
  q_p30_weak  = !is.na(rain_p30)   & !is.na(p30_weak_thr)  & (rain_p30   <= p30_weak_thr),

  q_cdd_fail  = !is.na(cdd_p45)    & !is.na(cdd_fail_thr)  & (cdd_p45    >= cdd_fail_thr),
  q_cdd_weak  = !is.na(cdd_p45)    & !is.na(cdd_weak_thr)  & (cdd_p45    >= cdd_weak_thr),

  q_rtot_fail = !is.na(rain_total) & !is.na(rtot_fail_thr) & (rain_total <= rtot_fail_thr),
  q_rtot_weak = !is.na(rain_total) & !is.na(rtot_weak_thr) & (rain_total <= rtot_weak_thr)
)]

# 4) Hard-limit flags (aridity-aware) #######
dat[, `:=`(
  hard_p30_fail  = FALSE,
  hard_p30_weak  = FALSE,
  hard_cdd_fail  = FALSE,
  hard_cdd_weak  = FALSE,
  hard_rtot_fail = FALSE,
  hard_rtot_weak = FALSE
)]

dat[aridity_bin %chin% c("hyper-arid","arid","semi-arid"), `:=`(
  hard_p30_fail  = !is.na(rain_p30)   & (rain_p30   < p30_hard_fail_arid),
  hard_p30_weak  = !is.na(rain_p30)   & (rain_p30   < p30_hard_weak_arid),

  hard_cdd_fail  = !is.na(cdd_p45)    & (cdd_p45    > cdd_hard_fail_arid),
  hard_cdd_weak  = !is.na(cdd_p45)    & (cdd_p45    > cdd_hard_weak_arid),

  hard_rtot_fail = !is.na(rain_total) & (rain_total < rtot_hard_fail_arid),
  hard_rtot_weak = !is.na(rain_total) & (rain_total < rtot_hard_weak_arid)
)]

# In humid/sub-humid: never hard-fail on dryness
dat[aridity_bin %chin% c("sub-humid","humid"), `:=`(
  hard_p30_fail  = FALSE,
  hard_cdd_fail  = FALSE,
  hard_rtot_fail = FALSE
)]

# 5) Combine gates (CDD can fail on its own) #######
dat[, gate_moist_fail := (q_cdd_fail | hard_cdd_fail) |
                         ((q_p30_fail | hard_p30_fail) & (q_rtot_fail | hard_rtot_fail))]

dat[, gate_moist_weak := (q_p30_weak | hard_p30_weak) |
                         (q_cdd_weak | hard_cdd_weak) |
                         (q_rtot_weak | hard_rtot_weak)]

# 6) Source-aware final QC #######
dat[, `:=`(season_fail = FALSE, season_weak = FALSE)]

dat[season_source == "greenup", `:=`(
  season_fail = gate_moist_fail & !ndvi_ok,
  season_weak = gate_moist_weak | (gate_moist_fail & ndvi_ok)
)]

dat[season_source %chin% c("der","trs2"), `:=`(
  season_fail = gate_moist_fail,
  season_weak = gate_moist_weak & !season_fail
)]

# 7) Labels (simple) #######
dat[, season_qc := fifelse(season_fail, "fail",
                    fifelse(season_weak, "weak", "ok"))]


# admin1 doy start of season averages
admin1_tab <- dat[!is.na(season_harmonized),
  {
    ci <- quantile_circular_safe(doy, probs = c(0.1, 0.5, 0.9))
    names(ci) <- paste0("sos_", names(ci))
    ci_eos <- quantile_circular_safe(doy_eos, probs = c(0.1, 0.5, 0.9))
    names(ci_eos) <- paste0("eos_", names(ci_eos))
    c(
      as.list(ci),
      as.list(ci_eos),
      list(
        rain_mean = mean(rain_trs2, na.rm = TRUE),
        slen_mean = mean(slen, na.rm = TRUE)
      )
    )
  },
  by = .(admin1_name, season_harmonized)
][, slen_calc :=  season_length_doy(`sos_50%`, `eos_50%`)]

# Labeling thresholds (to avoid over-interpreting tiny differences)
# If the relative difference between seasons is smaller than these thresholds,
# labels will be set to "similar" instead of forcing a max/min label.
label_diff_pct_rain   <-  0 # e.g. 0.05 = 5%, 0.10 = 10%
label_diff_pct_length <- 0


admin1_tab[, `:=`(
  # wetter / drier label
  rain_label = {
    if (.N < 2L || all(is.na(rain_mean))) {
      rep(NA_character_, .N)                        # only one season or no data
    } else {
      r <- rain_mean
      rmax <- max(r, na.rm = TRUE)
      rmin <- min(r, na.rm = TRUE)

      # relative difference (avoid dividing by 0)
      rel_diff <- if (is.finite(rmax) && rmax > 0) (rmax - rmin) / rmax else NA_real_

      if (!is.na(rel_diff) && rel_diff < label_diff_pct_rain) {
        rep("similar", .N)
      } else {
        wetter  <- which.max(r)
        drier   <- which.min(r)
        lab <- character(.N)
        lab[wetter] <- "wetter"
        lab[drier]  <- "drier"
        lab
      }
    }
  },

  # longer / shorter label (based on median-derived season length)
  length_label = {
    sl <- slen_mean
    if (.N < 2L || all(is.na(sl))) {
      rep(NA_character_, .N)
    } else {
      slmax <- max(sl, na.rm = TRUE)
      slmin <- min(sl, na.rm = TRUE)

      rel_diff <- if (is.finite(slmax) && slmax > 0) (slmax - slmin) / slmax else NA_real_

      if (!is.na(rel_diff) && rel_diff < label_diff_pct_length) {
        rep("similar", .N)
      } else {
        longer  <- which.max(sl)
        shorter <- which.min(sl)
        lab <- character(.N)
        lab[longer]  <- "longer"
        lab[shorter] <- "shorter"
        lab
      }
    }
  },

  # first vs second by timing (earlier/later in year, using median SOS)
  timing_label = {
    sos <- `sos_50%`
    if (.N < 2L || all(is.na(sos))) {
      rep(NA_character_, .N)
    } else {
      earlier <- which.min(sos)
      later   <- which.max(sos)
      lab <- character(.N)
      lab[earlier] <- "first"
      lab[later]   <- "second"
      lab
    }
  }
), by = admin1_name]

admin1_tab[,.N,by=.(rain_label,length_label,timing_label)]

admin1 <- "Kwale"

dat |>
  subset(admin1_name == admin1 & !is.na(season_harmonized) & !is.na(TRS2.sos)) |>
  transform(doy = yday(TRS2.sos),
            season_harmonized = factor(season_harmonized)) |>
  ggplot(aes(x = doy, fill = season_harmonized)) +
  geom_histogram(position = "identity", bins = 30, alpha = 0.45) +
  scale_x_continuous(limits = c(1, 366)) +
  labs(
    title = admin1,
    x = "Day of Year (TRS2.sos)",
    fill = "season_harmonized"
  )+
  theme_bw()

dat[!is.na(season_harmonized),
    {
      qs <- quantile_circular_safe(doy, probs = c(0.1, 0.5, 0.90))
      .(lower = qs[1], median = qs[2], upper = qs[3])
    },
    by = .(admin1_name, season_harmonized)]


dat[, n_seasons :=  .N / (length(unique(year)) * length(unique(pixel))), by = admin1_name # nolint: line_length_linter.
][, admin_N :=  length(unique(pixel)), by = admin1_name]

unique(dat[, .(admin1_name, n_seasons, admin_N)][order(n_seasons)])
