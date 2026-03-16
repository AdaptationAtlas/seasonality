pacman::p_load(circular, changepoint, data.table, arrow, lubridate, ggplot2, terra)
source("R/00_setup_folders.R")

# Create Functions ----
  # Find Automatic Split Point Between The Two Clusters ----
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
    probs = c(0.1, 0.5, 0.9),
    n_days = 365
  ) {
  doy <- doy[!is.na(doy)]
  if (!length(doy)) {
    out <- rep(NA_real_, length(probs))
    names(out) <- paste0(probs * 100, "%")
    return(out)
  }

  # circular centre from mean angle
  ang <- 2 * pi * (doy / n_days)
  mu <- atan2(mean(sin(ang)), mean(cos(ang)))
  if (mu < 0) mu <- mu + 2 * pi
  centre <- mu * n_days / (2 * pi)

  # unwrap around centre so values stay close together
  x <- doy
  x[x < centre - n_days / 2] <- x[x < centre - n_days / 2] + n_days
  x[x > centre + n_days / 2] <- x[x > centre + n_days / 2] - n_days

  qs <- quantile(x, probs = probs, na.rm = TRUE)

  # wrap back to day-of-year
  qs <- ((qs - 1) %% n_days) + 1
  qs
  }
  season_length_doy <- function(sos, eos, n_days = 365) {
    (eos - sos + n_days) %% n_days
  }

# Choose And Load Data ----
# Load Pixel Index ----
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

  # Load Or Compute DEM And Aridity Index By Pixel ----
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

# Load Phenology Data ----
  # Folder for country phenology data
  pheno_dat <- file.path(dirs$nvdi_phenology, "countries")

  # Choose ISO3 country code
  files<-list.files(pheno_dat, "plus-rain", full.names = TRUE)
  iso3_choice <- "KEN"
  file_choice <- grep(iso3_choice, files, value = TRUE)

  # Load data
  dat <- setDT(read_parquet(file_choice))

  # Merge DEM, Aridity, And Coordinates ----
  dat <- px_env[dat, on = .(pixel)]

  # Create combined zone ID
  dat[, zone_id := paste(admin1_name, aridity_bin, elev_bin, sep = "|")]

  # Prepare Phenology Data ----
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

  # Step 1: Choose Best Available Season Source Per Row ----

    # Priority: Greenup > DER > TRS2 (requires both sos/eos to be present)
    dat[, season_source := fifelse(!is.na(Greenup) & !is.na(Senescence), "greenup",  
    fifelse(!is.na(DER.sos) & !is.na(DER.eos), "der",
    fifelse(!is.na(TRS2.sos) & !is.na(TRS2.eos), "trs2", NA_character_)))
    ]

    # unify timing + rainfall metrics based on season_source
    dat[, sos := fifelse(season_source == "greenup", Greenup,
    fifelse(season_source == "der", DER.sos,
    fifelse(season_source == "trs2", TRS2.sos, as.IDate(NA))))
    ]

    dat[, eos := fifelse(season_source == "greenup", Senescence,
    fifelse(season_source == "der", DER.eos,
    fifelse(season_source == "trs2", TRS2.eos, as.IDate(NA))))]

    dat[, rain_total := fifelse(season_source == "greenup", rain_greenup_scenescence,
    fifelse(season_source == "der", rain_der,
    fifelse(season_source == "trs2", rain_trs2, NA_real_)))
    ]

    dat[, rain_p30 := fifelse(season_source == "greenup", rain_greenup_p30,
    fifelse(season_source == "der", rain_der_p30,
    fifelse(season_source == "trs2", rain_trs2_p30, NA_real_)))
    ]

    dat[, cdd_p45 := fifelse(season_source == "greenup", cdd_greenup_p45,
    fifelse(season_source == "der", cdd_der_p45,
    fifelse(season_source == "trs2", cdd_trs2_p45, NA_real_)))
    ]

    # Set season start / end dates
    dat[, doy :=  yday(sos)]
    dat[, doy_eos :=  yday(eos)]
    dat[, slen :=  as.numeric(eos - sos)]

    # Harmonize season numbering within each admin1 region
    dat[!is.na(sos),
      season_harmonized :=  split_seasons_cpt(x = sos, max_seasons = 2L),
      by = admin1_name
      ]

    # Add the length of the preceding observed season for each pixel.
    # This is based on chronological order of season start date within pixel.
    dat[order(pixel, sos, eos, year, season), preceding_slen := shift(slen), by = pixel]

# Step 2: Gates And QC Flags (Zone-Based, Aridity-Aware) 

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

# Step 2.0: Long season flags ----
  # Calculate the proportion of season 1 vs season 2 occurrences for each pixel across years
  dat[,season1_prop:=sum(season_harmonized==1,na.rm=TRUE)/length(unique(year)),by=.(admin1_name,pixel)]
  dat[,season2_prop:=sum(season_harmonized==2,na.rm=TRUE)/length(unique(year)),by=.(admin1_name,pixel)]
  dat[,season3_prop:=sum(season=3,na.rm=TRUE)/length(unique(year)),by=.(admin1_name,pixel)]

  # Add flags for long-preceding season which could push the onset of the following season to an odd time.
  dat[season2_prop>=0.2, s_long_preceding := fifelse(preceding_slen > 30*8, TRUE, FALSE)]
  dat[season2_prop<0.2, s_long_preceding := fifelse(preceding_slen > 30*11, TRUE, FALSE)]
  dat[season2_prop>=0.2, s_long := fifelse(slen > 30*8, TRUE, FALSE)]
  dat[season2_prop<0.2, s_long := fifelse(slen > 30*11, TRUE, FALSE)]

# Step 2.1: Zone-Specific Thresholds (zone_id + season_harmonized) ----
  zone_thr <- dat[!is.na(zone_id) & !is.na(season_harmonized),
    .(
      n_zone = .N,
      p30_fail_thr  = as.numeric(quantile(rain_p30,   p30_q_fail,  na.rm = TRUE)),
      p30_weak_thr  = as.numeric(quantile(rain_p30,   p30_q_weak,  na.rm = TRUE)),
      cdd_fail_thr  = as.numeric(quantile(cdd_p45,    cdd_q_fail,  na.rm = TRUE)),
      cdd_weak_thr  = as.numeric(quantile(cdd_p45,    cdd_q_weak,  na.rm = TRUE)),
      rtot_fail_thr = as.numeric(quantile(rain_total, rtot_q_fail, na.rm = TRUE)),
      rtot_weak_thr = as.numeric(quantile(rain_total, rtot_q_weak, na.rm = TRUE))
    ), by = .(zone_id, season_harmonized)
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

# Step 2.2: NDVI Reliability Indicator ----
  dat[, ndvi_ok := {
    has_fit <- !is.na(NSE) & !is.na(R2)
    fifelse(has_fit, (NSE >= ndvi_nse_min & R2 >= ndvi_r2_min), TRUE)
  }]

# Step 2.3: Quantile-Based Flags ----

  dat[, `:=`(
    q_p30_fail  = !is.na(rain_p30)   & !is.na(p30_fail_thr)  & (rain_p30   <= p30_fail_thr),
    q_p30_weak  = !is.na(rain_p30)   & !is.na(p30_weak_thr)  & (rain_p30   <= p30_weak_thr),

    q_cdd_fail  = !is.na(cdd_p45)    & !is.na(cdd_fail_thr)  & (cdd_p45    >= cdd_fail_thr),
    q_cdd_weak  = !is.na(cdd_p45)    & !is.na(cdd_weak_thr)  & (cdd_p45    >= cdd_weak_thr),

    q_rtot_fail = !is.na(rain_total) & !is.na(rtot_fail_thr) & (rain_total <= rtot_fail_thr),
    q_rtot_weak = !is.na(rain_total) & !is.na(rtot_weak_thr) & (rain_total <= rtot_weak_thr)
  )]

# Step 2.4: Hard-Limit Flags (Aridity-Aware) ----
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

# Step 2.5: Combine Gates (CDD Can Fail On Its Own) ----
  dat[, gate_moist_fail := (q_cdd_fail | hard_cdd_fail) |((q_p30_fail | hard_p30_fail) & (q_rtot_fail | hard_rtot_fail))]

  dat[, gate_moist_weak := (q_p30_weak | hard_p30_weak) |(q_cdd_weak | hard_cdd_weak) |(q_rtot_weak | hard_rtot_weak)]

# Step 2.6: Source-Aware Final QC ----
  dat[, `:=`(season_fail = FALSE, season_weak = FALSE)]

  dat[season_source == "greenup", `:=`(
    season_fail = gate_moist_fail & !ndvi_ok,
    season_weak = gate_moist_weak | (gate_moist_fail & ndvi_ok)
  )]

  dat[season_source %chin% c("der","trs2"), `:=`(
    season_fail = gate_moist_fail,
    season_weak = gate_moist_weak & !season_fail
  )]

# Step 3: Labels ----
  n_years<-length(unique(focal_dat$year))

  dat[, season_qc := fifelse(season_fail, "fail",
                    fifelse(season_weak, "weak", "ok"))]

  
  # admin1 doy start of season averages
  dat[,pixel_n:=length(unique(pixel))*length(unique(year)),by=.(admin1_name)]
  dat[,no_detect:=n_years-.N,by=.(admin1_name,pixel,season_harmonized)]
  dat[,failed:=sum(season_qc=="fail",na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]
  dat[,weak:=sum(season_qc=="weak",na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]
  dat[,ndvi_weak:=sum(!ndvi_ok,na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]

# Use all non NA pixels ----
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
          rain_mean = floor(mean(rain_total, na.rm = TRUE)),
          slen_mean = round(mean(slen, na.rm = TRUE),1),
          s_weak=round(sum(season_qc=="weak",na.rm=TRUE)/.N,3),
          s_fail=round(sum(season_qc=="fail",na.rm=TRUE)/.N,3),
          s_long_prec = round(sum(s_long_preceding,na.rm=TRUE)/.N,3),
          s_long = round(sum(s_long,na.rm=TRUE)/.N,3),
          s_pixel_prop = .N/pixel_n
        )
      )
    },
    by = .(admin1_name, season_harmonized,pixel_n)
  ][, slen_calc :=  season_length_doy(`sos_50%`, `eos_50%`)]

# Conservative Version ----
   dat_c<-copy(dat)

  dat_c<-dat_c[!is.na(season_harmonized) & s_long_preceding!=TRUE & season_qc!="fail"]
  # Set eos & rain to NA if the season is flagged as long, since the end date may be unreliable in this case (could be pushed by an excessively long preceding season).
  dat_c[s_long==TRUE,c("doy_eos", "rain_total"):=NA]

  admin1_tab_c <- dat_c[,
    {
      ci <- quantile_circular_safe(doy, probs = c(0.1, 0.5, 0.9))
      names(ci) <- paste0("sos_", names(ci))
      ci_eos <- quantile_circular_safe(doy_eos, probs = c(0.1, 0.5, 0.9))
      names(ci_eos) <- paste0("eos_", names(ci_eos))
      c(
        as.list(ci),
        as.list(ci_eos),
        list(
          rain_mean = floor(mean(rain_total, na.rm = TRUE)),
          slen_mean = round(mean(slen, na.rm = TRUE),1),
          s_weak=round(sum(season_qc=="weak",na.rm=TRUE)/.N,3),
          s_long = round(sum(s_long,na.rm=TRUE)/.N,3),
          s_pixel_prop = .N/pixel_n
        )
      )
    },
    by = .(admin1_name, season_harmonized,pixel_n)
  ][, slen_calc :=  season_length_doy(`sos_50%`, `eos_50%`)]

  # Join into a list
  admin1_tabs <- list(
    all = admin1_tab,
    conservative = admin1_tab_c
  )

n<-30
rbind(
  admin1_tab[n,!c("s_fail","s_long_prec")],
  admin1_tab_c[admin1_name==admin1_tab$admin1_name[n] & 
  season_harmonized==admin1_tab$season_harmonized[n]]
)

# Plot differences between all and conservative summaries ----
  target_vars <- c("sos", "eos", "slen")

  admin1_plot_all <- admin1_tab[
    season_harmonized == 1,
    .(
      admin1_name,
      dataset = "all",
      sos_10 = `sos_10%`,
      sos_50 = `sos_50%`,
      sos_90 = `sos_90%`,
      eos_10 = `eos_10%`,
      eos_50 = `eos_50%`,
      eos_90 = `eos_90%`,
      slen_10 =  NA,
      slen_50 = `slen_calc`,
      slen_90 = NA
    )
  ]

  admin1_plot_cons <- admin1_tab_c[
    season_harmonized == 1,
    .(
      admin1_name,
      dataset = "conservative",
      sos_10 = `sos_10%`,
      sos_50 = `sos_50%`,
      sos_90 = `sos_90%`,
      eos_10 = `eos_10%`,
      eos_50 = `eos_50%`,
      eos_90 = `eos_90%`,
      slen_10 = NA,
      slen_50 = `slen_calc`,
      slen_90 = NA
    )
  ]

  compare_data_wide <- rbindlist(
    list(admin1_plot_all, admin1_plot_cons),
    use.names = TRUE,
    fill = TRUE
  )

  unwrap_doy <- function(x, anchor, n_days = 365) {
    y <- x
    y[y < anchor - n_days / 2] <- y[y < anchor - n_days / 2] + n_days
    y[y > anchor + n_days / 2] <- y[y > anchor + n_days / 2] - n_days
    y
  }

  compare_data <- rbindlist(
    list(
      compare_data_wide[
        , .(
          admin1_name,
          dataset,
          metric = "SOS",
          q10 = sos_10,
          q50 = sos_50,
          q90 = sos_90
        )
      ],
      compare_data_wide[
        , .(
          admin1_name,
          dataset,
          metric = "EOS",
          q10 = eos_10,
          q50 = eos_50,
          q90 = eos_90
        )
      ],
      compare_data_wide[
        , .(
          admin1_name,
          dataset,
          metric = "Season length",
          q10 = NA,
          q50 = slen_50,
          q90 = NA
        )
      ]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  compare_data[, dataset := factor(dataset, levels = c("all", "conservative"))]
  compare_data[, metric := factor(metric, levels = c("SOS", "EOS", "Season length"))]

  # unwrap circular variables within each admin1-metric combination
  compare_data[
    metric %chin% c("SOS", "EOS"),
    anchor := median(q50, na.rm = TRUE),
    by = .(admin1_name, metric)
  ]

  compare_data[
    metric %chin% c("SOS", "EOS"),
    `:=`(
      q10_plot = unwrap_doy(q10, anchor),
      q50_plot = unwrap_doy(q50, anchor),
      q90_plot = unwrap_doy(q90, anchor)
    )
  ]

  # non-circular variable: keep original values
  compare_data[
    metric == "Season length",
    `:=`(
      q10_plot = NA,
      q50_plot = q50,
      q90_plot = NA
    )
  ]

  ggplot(
    compare_data,
    aes(
      x = metric,
      y = q50_plot,
      colour = dataset,
      shape = dataset
    )
  ) +
    geom_errorbar(
      aes(ymin = q10_plot, ymax = q90_plot),
      width = 0.15,
      linewidth = 0.7,
      position = position_dodge(width = 0.4)
    ) +
    geom_point(
      size = 2.4,
      position = position_dodge(width = 0.4)
    ) +
    facet_wrap(~admin1_name, scales = "free_y") +
    labs(
      title = "Season 1: all vs conservative summaries",
      x = NULL,
      y = "Median with 10-90% interval",
      colour = "Dataset",
      shape = "Dataset"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

# Label wetter/drier, longer/shorter, first/second by admin1 ----
  # Labeling thresholds (to avoid over-interpreting tiny differences)
  # If the relative difference between seasons is smaller than these thresholds,
  # labels will be set to "similar" instead of forcing a max/min label.
  label_diff_pct_rain   <-  0 # e.g. 0.05 = 5%, 0.10 = 10%
  label_diff_pct_length <- 0

  admin1_tabs<-lapply(admin1_tabs,function(x){
    x[, `:=`(
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
    return(x)
  })
  names(admin1_tabs) <- c("all", "conservative")

# Show differences between seasons in terms of rain and season length ----
      season_diff_base <- dat |>
      subset(
        admin1_name %in% admin1_units &
          !is.na(season_harmonized) &
          !is.na(year) &
          !is.na(pixel) &
          !is.na(rain_greenup_scenescence) &
          !is.na(Greenup) &
          !is.na(Senescence)
      ) |>
      transform(
        year = as.integer(year),
        season_length = as.numeric(Senescence - Greenup)
      ) |>
      setDT()

    # 1) Admin1-year means: one mean value per season per admin1-year
    season_diff_admin <- season_diff_base[
      , .(
        rain = mean(rain_greenup_scenescence, na.rm = TRUE),
        season_length = mean(season_length, na.rm = TRUE)
      ),
      by = .(admin1_name, year, season_harmonized)
    ][
      , dcast(
        .SD,
        admin1_name + year ~ season_harmonized,
        value.var = c("rain", "season_length")
      )
    ][
      !is.na(rain_1) & !is.na(rain_2) &
        !is.na(season_length_1) & !is.na(season_length_2)
    ][
      , `:=`(
        rain_diff = rain_2 - rain_1,
        norm_rain_diff = fifelse(
          (rain_1 + rain_2) > 0,
          (rain_2 - rain_1) / (rain_1 + rain_2),
          NA_real_
        ),
        season_length_diff = season_length_2 - season_length_1,
        method = "Admin-year mean"
      )
    ]

    # 2) Pixel-year means: one mean value per season per pixel-year
    season_diff_pixel <- season_diff_base[
      , .(
        rain = mean(rain_greenup_scenescence, na.rm = TRUE),
        season_length = mean(season_length, na.rm = TRUE)
      ),
      by = .(admin1_name, pixel, year, season_harmonized)
    ][
      , dcast(
        .SD,
        admin1_name + pixel + year ~ season_harmonized,
        value.var = c("rain", "season_length")
      )
    ][
      !is.na(rain_1) & !is.na(rain_2) &
        !is.na(season_length_1) & !is.na(season_length_2)
    ][
      , .(
        rain_diff = mean(rain_2 - rain_1, na.rm = TRUE),
        norm_rain_diff = mean(
          fifelse(
            (rain_1 + rain_2) > 0,
            (rain_2 - rain_1) / (rain_1 + rain_2),
            NA_real_
          ),
          na.rm = TRUE
        ),
        season_length_diff = mean(season_length_2 - season_length_1, na.rm = TRUE)
      ),
      by = .(admin1_name, year)
    ][
      , method := "Pixel-year mean"
    ]

    season_diff_all <- rbindlist(
      list(season_diff_admin, season_diff_pixel),
      use.names = TRUE,
      fill = TRUE
    )

    # Rainfall plot data: absolute and normalized versions
      rain_plot_dat_abs <- season_diff_all[
      !is.na(rain_diff),
      .(
        n = .N,
        mean_diff = mean(rain_diff),
        sd_diff = sd(rain_diff)
      ),
      by = .(method, admin1_name)
    ][
      , `:=`(
        se_diff = sd_diff / sqrt(n),
        ci_low = mean_diff - 1.96 * sd_diff / sqrt(n),
        ci_high = mean_diff + 1.96 * sd_diff / sqrt(n)
      )
    ]

    rain_levels_abs <- rain_plot_dat_abs[
      method == "Admin-year mean"
    ][
      order(mean_diff),
      unique(admin1_name)
    ]

    rain_plot_dat_abs[
      , admin1_name := factor(admin1_name, levels = rain_levels_abs)
    ]

    p_rain_abs <- ggplot(rain_plot_dat_abs, aes(x = admin1_name, y = mean_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
      coord_flip() +
      facet_wrap(~method, ncol = 2) +
      labs(
        title = "Season 2 - Season 1 rainfall difference by Admin1",
        x = "admin1_name",
        y = "Mean rainfall difference (mm)"
      ) +
      theme_bw()

   
    rain_plot_dat_norm <- season_diff_all[
      !is.na(norm_rain_diff),
      .(
        n = .N,
        mean_diff = mean(norm_rain_diff),
        sd_diff = sd(norm_rain_diff)
      ),
      by = .(method, admin1_name)
    ][
      , `:=`(
        se_diff = sd_diff / sqrt(n),
        ci_low = mean_diff - 1.96 * sd_diff / sqrt(n),
        ci_high = mean_diff + 1.96 * sd_diff / sqrt(n)
      )
    ]

    rain_levels_norm <- rain_plot_dat_norm[
      method == "Admin-year mean"
    ][
      order(mean_diff),
      unique(admin1_name)
    ]

    rain_plot_dat_norm[
      , admin1_name := factor(admin1_name, levels = rain_levels_norm)
    ]

    p_rain_norm <- ggplot(rain_plot_dat_norm, aes(x = admin1_name, y = mean_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_col(fill = "seagreen3", alpha = 0.8) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
      coord_flip() +
      facet_wrap(~method, ncol = 2) +
      labs(
        title = "Season 2 - Season 1 rainfall difference normalized by annual total",
        x = "admin1_name",
        y = "Mean normalized rainfall difference"
      ) +
      theme_bw()

    gridExtra::grid.arrange(
      p_rain_abs,
      p_rain_norm,
      ncol = 1
    )

    slen_plot_dat <- season_diff_all[
      , .(
        n = sum(!is.na(season_length_diff)),
        mean_diff = mean(season_length_diff, na.rm = TRUE),
        sd_diff = sd(season_length_diff, na.rm = TRUE)
      ),
      by = .(method, admin1_name)
    ][
      , `:=`(
        se_diff = sd_diff / sqrt(n),
        ci_low = mean_diff - 1.96 * se_diff,
        ci_high = mean_diff + 1.96 * se_diff
      )
    ]

    slen_levels <- slen_plot_dat[
      order(mean_diff),
      unique(admin1_name)
    ]

    slen_plot_dat[
      , admin1_name := factor(admin1_name, levels = slen_levels)
    ]

    ggplot(slen_plot_dat, aes(x = admin1_name, y = mean_diff)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_col(fill = "darkorange", alpha = 0.8) +
      geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
      coord_flip() +
      facet_wrap(~method) +
      labs(
        title = "Season 2 - Season 1 season length difference by Admin1",
        x = "admin1_name",
        y = "Mean season length difference (days)"
      ) +
      theme_bw()



# Circular helper functions ----
    circ_dist <- function(x, y, n_days = 365) {
      d <- abs(x - y)
      pmin(d, n_days - d)
    }

    assign_to_focus <- function(doy, sos_50_1, sos_50_2, n_days = 365) {
      d1 <- circ_dist(doy, sos_50_1, n_days)
      d2 <- circ_dist(doy, sos_50_2, n_days)
      fifelse(d1 <= d2, 1L, 2L)
    }

  forward_circ_dist <- function(start, end, n_days = 365) {
    start <- rep_len(start, max(length(start), length(end)))
    end   <- rep_len(end,   max(length(start), length(end)))
    (end - start + n_days) %% n_days
  }

  in_circular_window <- function(x, start, end, n_days = 365) {
    n <- max(length(x), length(start), length(end))

    x     <- rep_len(x, n)
    start <- rep_len(start, n)
    end   <- rep_len(end, n)

    out <- rep(NA, n)

    ok <- !is.na(x) & !is.na(start) & !is.na(end)
    nowrap <- ok & (start <= end)
    wrap   <- ok & (start > end)

    out[nowrap] <- x[nowrap] >= start[nowrap] & x[nowrap] <= end[nowrap]
    out[wrap]   <- x[wrap] >= start[wrap] | x[wrap] <= end[wrap]

    out
  }

  interval_within_circular_window <- function(doy, doy_eos, sos, eos, n_days = 365) {
    n <- max(length(doy), length(doy_eos), length(sos), length(eos))

    doy     <- rep_len(doy, n)
    doy_eos <- rep_len(doy_eos, n)
    sos     <- rep_len(sos, n)
    eos     <- rep_len(eos, n)

    out <- rep(NA, n)

    ok <- !is.na(doy) & !is.na(doy_eos) & !is.na(sos) & !is.na(eos)
    if (!any(ok)) return(out)

    start_ok <- in_circular_window(doy[ok], sos[ok], eos[ok], n_days = n_days)
    cand_len <- forward_circ_dist(doy[ok], doy_eos[ok], n_days = n_days)
    max_len  <- forward_circ_dist(doy[ok], eos[ok], n_days = n_days)

    out[ok] <- start_ok & (cand_len <= max_len)
    out
  }

# Create pixel-level maps ----
focal_dat<-copy(dat_c) 
focal_admin<-copy(admin1_tab_c)

# Remove % from summary column names for safer joins and downstream handling
focal_admin <- copy(focal_admin)

setnames(
  focal_admin,
  old = c("sos_10%", "sos_50%", "sos_90%", "eos_10%", "eos_50%", "eos_90%"),
  new = c("sos_10", "sos_50", "sos_90", "eos_10", "eos_50", "eos_90"),
  skip_absent = TRUE
)

# Keep only the columns needed for the join, then reshape to wide so each row gets
# the admin1 summaries for both season 1 and season 2.
focal_admin_join <- unique(
  focal_admin[
    , .(admin1_name, season_harmonized, sos_10, sos_50, sos_90, eos_10, eos_50, eos_90)
  ]
)

focal_admin_wide <- dcast(
  focal_admin_join,
  admin1_name ~ season_harmonized,
  value.var = c("sos_10", "sos_50", "sos_90", "eos_10", "eos_50", "eos_90")
)

# Merge the wide admin1 summaries into the pixel-level data so each row carries
# season-1 and season-2 reference values regardless of its own season label.
focal_dat[
  focal_admin_wide,
  on = .(admin1_name),
  `:=`(
    s1_sos_10 = i.sos_10_1,
    s1_sos_50 = i.sos_50_1,
    s1_sos_90 = i.sos_90_1,
    s1_eos_10 = i.eos_10_1,
    s1_eos_50 = i.eos_50_1,
    s1_eos_90 = i.eos_90_1,
    s2_sos_10 = i.sos_10_2,
    s2_sos_50 = i.sos_50_2,
    s2_sos_90 = i.sos_90_2,
    s2_eos_10 = i.eos_10_2,
    s2_eos_50 = i.eos_50_2,
    s2_eos_90 = i.eos_90_2
  )
]

# Flag rows where the observed SOS or EOS falls entirely within the other season's
# 10-90% timing window and not within its own season's 10-90% timing window.
focal_dat[
  , `:=`(
    own_sos_match = fifelse(
      season_harmonized == 1L,
      in_circular_window(doy, s1_sos_10, s1_sos_90),
      fifelse(
        season_harmonized == 2L,
        in_circular_window(doy, s2_sos_10, s2_sos_90),
        NA
      )
    ),
    other_sos_match = fifelse(
      season_harmonized == 1L,
      in_circular_window(doy, s2_sos_10, s2_sos_90),
      fifelse(
        season_harmonized == 2L,
        in_circular_window(doy, s1_sos_10, s1_sos_90),
        NA
      )
    ),
    own_eos_match = fifelse(
      season_harmonized == 1L,
      in_circular_window(doy_eos, s1_eos_10, s1_eos_90),
      fifelse(
        season_harmonized == 2L,
        in_circular_window(doy_eos, s2_eos_10, s2_eos_90),
        NA
      )
    ),
    other_eos_match = fifelse(
      season_harmonized == 1L,
      in_circular_window(doy_eos, s2_eos_10, s2_eos_90),
      fifelse(
        season_harmonized == 2L,
        in_circular_window(doy_eos, s1_eos_10, s1_eos_90),
        NA
      )
    )
  )
]

focal_dat[
  , `:=`(
    sos_mismatch = !is.na(own_sos_match) & !is.na(other_sos_match) & !own_sos_match & other_sos_match,
    eos_mismatch = !is.na(own_eos_match) & !is.na(other_eos_match) & !own_eos_match & other_eos_match
  )
]

focal_dat[, c("own_sos_match", "other_sos_match", "own_eos_match", "other_eos_match") := NULL]

# Remove rows with a mismatch
focal_dat <- focal_dat[,mismatch:=sum(sos_mismatch | eos_mismatch),by=.(admin1_name,pixel,season_harmonized)][!(sos_mismatch | eos_mismatch)]

# Histograms of greenup-based SOS by admin1 and season ----

focal_dat |>
    subset(
        !is.na(season_harmonized) &
        !is.na(Greenup)
    ) |>
    transform(
      doy = yday(Greenup),
      admin1_name = factor(admin1_name),
      season_harmonized = factor(season_harmonized)
    ) |>
    ggplot(aes(x = doy, fill = season_harmonized)) +
    geom_histogram(aes(y = after_stat(density)),position = "identity", bins = 30, alpha = 0.45) +
    facet_wrap(~admin1_name) +
    scale_x_continuous(limits = c(1, 366)) +
    labs(
      title = "Greenup Start Of Season By Admin1",
      x = "Day of Year (greenup-based)",
      fill = "season_harmonized"
    ) +
    theme_bw()

# Average focal_dat over years including CIs at 10,50,90% quantiles, by admin1 and season ----
n_years<-length(unique(focal_dat$year))

  focal_dat_avg <- focal_dat[
    !is.na(season_harmonized) &
      !is.na(doy) &
      !is.na(doy_eos),
    {
      sos_ci <- quantile_circular_safe(doy, probs = c(0.1, 0.5, 0.9))
      names(sos_ci) <- paste0("sos_", names(sos_ci))
      eos_ci <- quantile_circular_safe(doy_eos, probs = c(0.1, 0.5, 0.9))
      names(eos_ci) <- paste0("eos_", names(eos_ci))
      rain_ci<-quantile(rain_total, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
      names(rain_ci) <- paste0("rain_", names(rain_ci))
      slen_ci<-quantile(slen, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
      names(slen_ci) <- paste0("slen_", names(slen_ci))

      .(
        sos_10 = sos_ci[1],
        sos_50 = sos_ci[2],
        sos_90 = sos_ci[3],
        eos_10 = eos_ci[1],
        eos_50 = eos_ci[2],
        eos_90 = eos_ci[3],
        rain_10 = rain_ci[1],
        rain_50 = rain_ci[2],
        rain_90 = rain_ci[3],
        slen_10 = slen_ci[1],
        slen_50 = slen_ci[2],
        slen_90 = slen_ci[3],
        n=.N,
        n_prop = .N/n_years,
        failed_seasons = mean(failed,na.rm=TRUE),
        weak_seasons = sum(season_qc=="weak",na.rm=TRUE),
        mismatch = mean(mismatch,na.rm=TRUE),
        no_detect = mean(no_detect,na.rm=TRUE),
        ndvi_ok =  sum(ndvi_ok,na.rm=TRUE)
      )
    },
    by = .(admin1_name,pixel, season_harmonized)
  ]

focal_dat_avg[,failed_prop:=failed_seasons/(n_years-mismatch-no_detect)
][,weak_prop:=weak_seasons/(n_years-mismatch-no_detect)
][,mismatch_prop:=mismatch/(n_years-failed_seasons-no_detect)
][,no_detect_prop:=no_detect/n_years
][,ndvi_ok_prop:=ndvi_ok/(n_years-mismatch-no_detect)]

# Merge lat-lon from pixel index
  focal_dat_avg <- merge(
    focal_dat_avg,
    coords_index,
    by = "pixel",
    all.x = TRUE
  )

  focal_dat_avg[is.na(eos_90) & !is.na(eos_10)]
  bursummary(focal_dat)

# Make variables long (sos,eos,rain,slen) for easier plotting
  focal_dat_avg_long <- melt(
    focal_dat_avg[,!c("failed_seasons", "weak_seasons", "mismatch", "no_detect","ndvi_ok")],
    id.vars = c("admin1_name", "pixel", "season_harmonized","x", "y")
  )

# Make a raster stack for each season where the layers are variables.
  focal_dat_avg_long[
    , variable := factor(variable)
  ]

  raster_list <- list()

  for (season in unique(focal_dat_avg_long$season_harmonized)) {
    season_data <- focal_dat_avg_long[season_harmonized == season]

    # Create a terra SpatRaster with one layer per variable
    var_levels <- unique(as.character(season_data$variable))

    r_list <- lapply(var_levels, function(v) {
      # cat(season,v,"\n")
      v_dat <- season_data[variable == v, .(x, y, value)]
      r_v <- terra::rast(v_dat, type = "xyz", crs = "EPSG:4326")
      names(r_v) <- v
      r_v
    })

    r <- terra::rast(r_list)
    names(r) <- paste0("season",season,"_",var_levels)

    raster_list[[as.character(season)]] <- r
  }

names(raster_list) <- paste("season", names(raster_list), sep = "_")

# Plot palettes and helpers ----
  # Keep these helpers self-contained so they still work when run out of order in the terminal.
  get_circ_pal <- function(palette = c("hcl_soft", "hcl_month", "phenology"), n = 365) {
    palette <- match.arg(palette)

    if (palette == "hcl_soft") {
      # Balanced continuous cyclic HCL palette
      return(grDevices::hcl(
        h = seq(0, 360, length.out = n + 1)[1:n],
        c = 75,
        l = 62
      ))
    }

    if (palette == "hcl_month") {
      # Discrete month-like cyclic palette, repeated across the year
      month_cols <- grDevices::hcl(
        h = seq(15, 345, length.out = 12),
        c = 70,
        l = 64
      )
      return(rep(month_cols, length.out = n))
    }

    if (palette == "phenology") {
      # Phenology-style seasonal wheel: muted but more interpretable month separation
      # Jan-Feb cool blues -> Mar-May greens -> Jun-Aug yellows/oranges -> Sep-Dec reds/purples
      anchors <- c(
        "#4C78A8", "#5AA9E6", "#4ECDC4", "#7BC96F",
        "#C7E77F", "#F4D35E", "#EE964B", "#F95738",
        "#D1495B", "#A05195", "#6C5B7B", "#4C78A8"
      )
      ramp <- grDevices::colorRampPalette(anchors, space = "Lab")
      return(ramp(n))
    }
  }

  get_seq_pal <- function(palette = c("ylgnbu", "viridis", "magma", "terrain"), n = 100) {
    palette <- match.arg(palette)

    if (palette == "ylgnbu") {
      return(hcl.colors(n, palette = "YlGnBu", rev = TRUE))
    }

    if (palette == "viridis") {
      return(viridisLite::viridis(n))
    }

    if (palette == "magma") {
      return(viridisLite::magma(n))
    }

    if (palette == "terrain") {
      return(terrain.colors(n))
    }
  }

    month_at  <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
    month_lab <- month.abb

    is_circular_layer <- function(layer_name) {
      grepl("(^|_)(sos|eos)_", layer_name)
    }

  plot_circular_layer <- function(r, main = NULL, circ_pal = get_circ_pal("phenology")) {
    terra::plot(
      r,
      col = circ_pal,
      range = c(1, 366),
      main = main,
      plg = list(at = month_at, labels = month_lab)
    )
  }

  plot_normal_layer <- function(r, main = NULL, seq_pal = get_seq_pal("ylgnbu")) {
    terra::plot(
      r,
      col = seq_pal,
      main = main
    )
  }

   plot_season_stack <- function(
    r_stack,
    season_name = NULL,
    circ_palette = c("phenology", "hcl_soft", "hcl_month"),
    seq_palette = c("ylgnbu", "viridis", "magma", "terrain")
  ) {
    circ_palette <- match.arg(circ_palette)
    seq_palette <- match.arg(seq_palette)

    layer_names <- names(r_stack)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    n_layers <- terra::nlyr(r_stack)
    ncol_plot <- if (n_layers <= 4) 2 else 3
    nrow_plot <- ceiling(n_layers / ncol_plot)
    par(mfrow = c(nrow_plot, ncol_plot), mar = c(2.5, 2.5, 2.5, 4))

    circ_pal <- get_circ_pal(circ_palette)
    seq_pal <- get_seq_pal(seq_palette)

    # For selected non-circular variables, keep legend limits fixed across all
    # layers in the provided stack that share the same underlying variable.
    base_layer_names <- sub("^season[0-9]+_", "", layer_names)
    base_var_names <- sub("_(10|50|90)$", "", base_layer_names)
    shared_scale_max <- setNames(vector("list", length(layer_names)), layer_names)

    for (i in seq_along(layer_names)) {
      base_name <- base_var_names[i]

      if (base_name %in% c("rain", "slen")) {
        same_var_idx <- which(base_var_names == base_name)

        same_var_vals <- unlist(lapply(layer_names[same_var_idx], function(nm) {
          terra::values(r_stack[[nm]], mat = FALSE)
        }))

        shared_scale_max[[layer_names[i]]] <- max(same_var_vals, na.rm = TRUE)
      } else {
        shared_scale_max[[layer_names[i]]] <- NA_real_
      }
    }

    for (layer_name in layer_names) {
      layer_r <- r_stack[[layer_name]]
      layer_title <- if (is.null(season_name)) layer_name else paste(season_name, layer_name, sep = " - ")
      base_name <- sub("_(10|50|90)$", "", sub("^season[0-9]+_", "", layer_name))

      if (is_circular_layer(layer_name)) {
        plot_circular_layer(layer_r, main = layer_title, circ_pal = circ_pal)

      } else if (grepl("prop", layer_name)) {
        terra::plot(
          layer_r,
          col = seq_pal,
          range = c(0, 1),
          main = layer_title
        )

      } else if (base_name %in% c("rain", "slen") && is.finite(shared_scale_max[[layer_name]])) {
        terra::plot(
          layer_r,
          col = seq_pal,
          range = c(0, shared_scale_max[[layer_name]]),
          main = layer_title
        )

      } else {
        plot_normal_layer(layer_r, main = layer_title, seq_pal = seq_pal)
      }
    }
  }

# Plot raster outputs ----

 # Palette options for plot_season_stack()
  # --------------------------------------
  # circ_palette (for circular variables such as SOS and EOS):
  #   "phenology"  : recommended default; seasonal colour wheel similar to
  #                   palettes used in phenology papers (blue→green→yellow→red→purple)
  #   "hcl_soft"   : smooth continuous cyclic HCL palette, good for gradients
  #   "hcl_month"  : discrete month-like colours (12 repeating hues)
  #
  # seq_palette (for non‑circular variables such as rain, slen, proportions):
  #   "ylgnbu"     : sequential blue‑green palette (default; good for rainfall)
  #   "viridis"    : perceptually uniform scientific palette
  #   "magma"      : dark‑to‑light palette with strong contrast
  #   "terrain"    : classic terrain-style palette
  #
  # Example usage:
  #   plot_season_stack(raster_list$season_1, circ_palette = "phenology")
  #   plot_season_stack(raster_list$season_1, circ_palette = "hcl_soft", seq_palette = "viridis")

plot_season_stack(raster_list$season_1, season_name = NULL, circ_palette = "phenology", seq_palette = "ylgnbu")
plot_season_stack(raster_list$season_2, season_name = NULL, circ_palette = "phenology", seq_palette = "ylgnbu")

x <- c(raster_list$season_1[[grep("sos", names(raster_list$season_1))]],
       raster_list$season_2[[grep("sos", names(raster_list$season_2))]])
plot_season_stack(x, season_name = NULL, circ_palette = "phenology", seq_palette = "ylgnbu")

x <- c(raster_list$season_1[[grep("eos", names(raster_list$season_1))]],
       raster_list$season_2[[grep("eos", names(raster_list$season_2))]])
plot_season_stack(x, season_name = NULL, circ_palette = "phenology", seq_palette = "ylgnbu")

