# TO DO:

# DZA -> clearly there is an issue with merging close together season at the year end boundary,
# we can see them being split between Dec/Jan but merged elsewhere in the year
# this is an issue with section 1.1
# ERI - Debubawi Keih Bahri not playing ball


pacman::p_load(circular,changepoint, data.table,  arrow, lubridate,  ggplot2,  terra)

source("R/0-1_setup_folders.R")

# Create Functions ----
{

  source("R/functions/circular_seasons.r")

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
    o  <- order(doy_raw)
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

  quantile_circular_safe <- function(doy,probs = c(0.1, 0.5, 0.9),n_days = 365) {
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

  circ_dist <- function(x, y, n_days = 365) {
    d <- abs(x - y)
    pmin(d, n_days - d)
  }

  # circular helpers
  circ_mean_doy <- function(x, period = 365) {
    a <- 2 * pi * x / period
    ang <- atan2(mean(sin(a), na.rm = TRUE), mean(cos(a), na.rm = TRUE))
    if (ang < 0) ang <- ang + 2 * pi
    ang * period / (2 * pi)
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
}
# Plot palettes and helpers ----
{
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

  mask_raster_list_by_nprop <- function(raster_list, n_prop_min = 0.1) {

    out <- lapply(names(raster_list), function(snm) {
      r <- raster_list[[snm]]
      season_num <- sub("^season_", "", snm)

      n_prop_name <- paste0("season", season_num, "_n_prop")
      core_pattern <- paste0("^season", season_num, "_(sos|eos|rain|slen)_")

      if (!n_prop_name %in% names(r)) {
        stop(sprintf("Layer '%s' not found in raster_list[['%s']]", n_prop_name, snm))
      }

      # mask condition: keep only where n_prop >= threshold
      mask_ok <- r[[n_prop_name]] >= n_prop_min

      # core layers to mask
      core_idx <- grep(core_pattern, names(r), ignore.case = TRUE)
      core_layers <- r[[core_idx]]

      core_layers_masked <- terra::mask(
        core_layers,
        mask_ok,
        maskvalues = 0,
        updatevalue = NA
      )

      # keep non-core layers unchanged
      other_idx <- setdiff(seq_along(names(r)), core_idx)
      other_layers <- if (length(other_idx) > 0) r[[other_idx]] else NULL

      if (is.null(other_layers)) {
        core_layers_masked
      } else {
        c(core_layers_masked, other_layers)
      }
    })

    names(out) <- names(raster_list)
    out
  }

  build_season_stack <- function(
    raster_list,
    vars = c("sos"),   # any of: sos, eos, rain, slen
    include_n = TRUE
  ) {

    vars <- tolower(vars)
    valid_vars <- c("sos", "eos", "rain", "slen")

    if (!all(vars %in% valid_vars)) {
      stop("vars must be any of: 'sos', 'eos', 'rain', 'slen'")
    }

    season_ids <- as.integer(sub("^season_", "", names(raster_list)))
    season_names <- names(raster_list)[order(season_ids)]

    x <- do.call(
      c,
      lapply(season_names, function(snm) {
        r <- raster_list[[snm]]
        season_num <- sub("^season_", "", snm)

        # match requested variables
        var_pattern <- paste0("_(", paste(vars, collapse = "|"), ")_")
        var_idx <- grep(var_pattern, names(r), ignore.case = TRUE)

        if (length(var_idx) == 0L) {
          stop(sprintf("No matching layers (%s) in %s", paste(vars, collapse = ","), snm))
        }

        layers <- r[[var_idx]]

        # optionally include n layer
        if (include_n) {
          n_name <- paste0("season", season_num, "_n")
          if (!n_name %in% names(r)) {
            stop(sprintf("Layer '%s' not found in %s", n_name, snm))
          }
          layers <- c(layers, r[[n_name]])
        }

        layers
      })
    )

    # move *_n layers to end
    if (include_n) {
      nm <- names(x)
      is_n <- grepl("_n$", nm)
      x <- x[[c(which(!is_n), which(is_n))]]
    }

    x
  }

  plot_season_vars <- function(
    raster_list,
    vars = c("sos"),
    include_n = TRUE,
    season_name = NULL,
    circ_palette = "phenology",
    seq_palette = "ylgnbu"
  ) {

    x <- build_season_stack(
      raster_list = raster_list,
      vars = vars,
      include_n = include_n
    )

    plot_season_stack(
      x,
      season_name = season_name,
      circ_palette = circ_palette,
      seq_palette = seq_palette
    )

    invisible(x)
  }

}
# Choose And Load Data ----
# Load Pixel Index ----
{
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
}
# Set custom overrides to fix season sequencing ----
{
  season_overrides <- data.table::data.table(
    iso3 = c("AGO"),
    admin1_name = c(NA_character_),   # NA = apply to all admin1 in country (or specify names)
    rule = c("sos_window_swap"),
    season_target = c(1L),
    window_start = c(5*30),  # Jul
    window_end   = c(9*30)   # Sept
  )

  apply_season_overrides <- function(dat, overrides, iso3_selected) {

    if (is.null(overrides) || nrow(overrides) == 0) return(dat)

    ov <- overrides[iso3 == iso3_selected]
    if (nrow(ov) == 0) return(dat)

    for (i in seq_len(nrow(ov))) {

      rule <- ov$rule[i]
      target <- ov$season_target[i]
      start <- ov$window_start[i]
      end   <- ov$window_end[i]
      admin_filter <- ov$admin1_name[i]

      # restrict to admin if provided
      dat_sub <- if (is.na(admin_filter)) {
        dat
      } else {
        dat[admin1_name == admin_filter]
      }

      if (rule == "sos_window_swap") {

        map <- dat_sub[
          !is.na(season_harmonized),
          .(sos_med = quantile_circular_safe(doy, probs = 0.5)),
          by = .(admin1_name, season_harmonized)
        ]

        map[
          season_harmonized == target &
          in_circular_window(sos_med, start, end),
          flip := TRUE
        ]

        flip_admins <- unique(map[flip == TRUE, admin1_name])

        if (length(flip_admins) > 0) {
          dat[
            admin1_name %in% flip_admins,
            season_harmonized := 3L - season_harmonized
          ]
        }
    }
   }

  dat
  }
}
# Load Phenology Data ----
  # Folder for country phenology data
  pheno_dat <- file.path(dirs$nvdi_phenology, "countries")
  # Choose ISO3 country code
  files<-list.files(pheno_dat, "plus-rain", full.names = TRUE)

  iso3_choices<-unlist(tstrsplit(basename(files),"_",keep=1))

# Set parameters ----
  # use DER
  use_DER<-F

  #  Parameters
  {
    # QC parameter definitions ---------------------------------------------------
    #
    # Quantile-based thresholds are computed within zone_id × season_harmonized.
    # They are relative thresholds, not absolute physical values.
    #
    # Variable definitions:
    #   p30         = rainfall accumulated in the first 30 days of the season
    #                 Unit: mm
    #   rtot        = total rainfall over the season
    #                 Unit: mm
    #   cdd         = consecutive dry days over the focal seasonal window
    #                 Unit: days
    #   NSE         = Nash-Sutcliffe Efficiency of NDVI seasonal fit
    #                 Unit: unitless, typically (-Inf, 1]
    #   R2          = coefficient of determination of NDVI seasonal fit
    #                 Unit: unitless, [0, 1]
    #
    # Quantile gates (zone-specific thresholds):
    #   p30_q_fail / p30_q_weak
    #       Meaning: lower-tail quantiles of early-season rainfall (p30).
    #       Interpretation:
    #         p30 <= zone-specific p30_q_fail threshold  -> strong failure signal
    #         p30 <= zone-specific p30_q_weak threshold  -> weak stress signal
    #       Unit of underlying variable: mm
    #       Parameter unit: probability / quantile in [0, 1]
    #       Higher values = stricter (more seasons flagged)
    #
    #   rtot_q_fail / rtot_q_weak
    #       Meaning: lower-tail quantiles of total seasonal rainfall (rtot).
    #       Interpretation:
    #         rtot <= zone-specific rtot_q_fail threshold -> strong failure signal
    #         rtot <= zone-specific rtot_q_weak threshold -> weak stress signal
    #       Unit of underlying variable: mm
    #       Parameter unit: probability / quantile in [0, 1]
    #       Higher values = stricter
    #
    #   cdd_q_fail / cdd_q_weak
    #       Meaning: upper-tail quantiles of consecutive dry days (CDD).
    #       Interpretation:
    #         cdd >= zone-specific cdd_q_fail threshold -> strong failure signal
    #         cdd >= zone-specific cdd_q_weak threshold -> weak stress signal
    #       Unit of underlying variable: days
    #       Parameter unit: probability / quantile in [0, 1]
    #       Lower values = stricter, because they pull the threshold downward
    #
    # Absolute hard limits (applied only in hyper-arid / arid / semi-arid bins):
    #   p30_hard_fail_arid / p30_hard_weak_arid
    #       Meaning: absolute minimum early rainfall in first 30 days
    #       Unit: mm
    #       Interpretation:
    #         p30 < threshold -> fail / weak
    #       Higher values = stricter
    #
    #   rtot_hard_fail_arid / rtot_hard_weak_arid
    #       Meaning: absolute minimum total seasonal rainfall
    #       Unit: mm
    #       Interpretation:
    #         rtot < threshold -> fail / weak
    #       Higher values = stricter
    #
    #   cdd_hard_fail_arid / cdd_hard_weak_arid
    #       Meaning: absolute maximum dry-spell length
    #       Unit: days
    #       Interpretation:
    #         cdd > threshold -> fail / weak
    #       Lower values = stricter
    #
    # NDVI fit quality thresholds:
    #   ndvi_nse_min
    #       Meaning: minimum acceptable Nash-Sutcliffe Efficiency for NDVI fit
    #       Unit: unitless
    #       Higher values = stricter
    #
    #   ndvi_r2_min
    #       Meaning: minimum acceptable R² for NDVI fit
    #       Unit: unitless
    #       Higher values = stricter
    #
    # Sample-size threshold:
    #   min_n
    #       Meaning: minimum number of observations required within a
    #       zone_id × season_harmonized group to compute stable quantile thresholds
    #       Unit: count of observations
    #       Higher values = more conservative / more stable, but less coverage

    p30_q_fail  <- 0.10
    p30_q_weak  <- 0.20
    cdd_q_fail  <- 0.90
    cdd_q_weak  <- 0.80
    rtot_q_fail <- 0.10
    rtot_q_weak <- 0.20

    p30_hard_fail_arid <- 5
    p30_hard_weak_arid <- 15

    cdd_hard_fail_arid <- 35
    cdd_hard_weak_arid <- 25

    rtot_hard_fail_arid <- 50
    rtot_hard_weak_arid <- 120


    ndvi_nse_min <- 0.60
    ndvi_r2_min  <- 0.60

    min_n <- 30

    params<-list()
  params$params<-list(
    use_DER = use_DER,
    p30_q_fail  = p30_q_fail,
    p30_q_weak  = p30_q_weak,
    cdd_q_fail  = cdd_q_fail,
    cdd_q_weak  = cdd_q_weak,
    rtot_q_fail = rtot_q_fail,
    rtot_q_weak = rtot_q_weak,
    p30_hard_fail_arid = p30_hard_fail_arid,
    p30_hard_weak_arid = p30_hard_weak_arid,
    cdd_hard_fail_arid = cdd_hard_fail_arid,
    cdd_hard_weak_arid = cdd_hard_weak_arid,
    rtot_hard_fail_arid = rtot_hard_fail_arid,
    rtot_hard_weak_arid = rtot_hard_weak_arid,
    ndvi_nse_min = ndvi_nse_min,
    ndvi_r2_min  = ndvi_r2_min,
    min_n = min_n
  )

# Set run version & output dirs----
version<-Sys.Date()

output_dir <- file.path(dirs$output,version)
output_dir_rast <- file.path(dirs$output,version,"rast")

if (!dir.exists(output_dir_rast)) dir.create(output_dir_rast, recursive = TRUE, showWarnings = FALSE)
}

# Loop through countries ----

admin_3_season<-list(AGO=c("Uíge","Lunda Norte","Lunda Sul"),
                     COD = c("Lomami","Kasaï-Central","Kinshasa","Kwilu","Kasaï-Central",
                             "Maniema","Sud-Kivu"),
                     CAF = "Equateur",
                     DJI=c("Obock","Djiboutii","Ali Sabieh","Arta","Tadjoura","Dikhil"),
                     COG=c("Pool"),
                     EGY= c("Bīr Ṭawīl"),
                     ERI =c("Debubawi Keih Bahri"),
                     RWA = c("Iburasirazuba","Amajyaruguru","Amajyepfo","Iburengerazuba","Umujyi Wa Kigali"),
                     SOM = c("Awdal"),
                     TZA = c("Tanga","Mara"))

skip_season_harmonization_admin1<-list(AGO= c("Luanda","Namibe"),
                                       BDI = c("Bujumbura","Rumonge","Bubanza"),
                                       BEN = c("Borgou", "Donga","Alibori","Atacora"),
                                       COG= c("Kouilou","Lekoumou","Niari"),
                                       CIV = c("Woroba","Savanes","Denguélé"),
                                       CMR = c("Extrême-Nord","Nord","Adamaoua"),
                                       DZA = c("Adrar","Alger","Biskra","Batna","M'Sila","Djelfa","El Bayadh",
                                               "Illizi","Khenchela","Laghouat","Naama","Tamanrasset",
                                               "Tebessa","Tiaret","Ghardaia","Bechar","El Oued","Guelma","Medea",
                                               "Mila","Setif","Souk-Ahras","Oum El Bouaghi","Sidi Bel Abbes",
                                               "Tindouf","Tissemsilt","Tlemcen","Ouargla","Constantine","Bordj Bou Arrer"),
                                       EGY = c("Aswan","North Sinai","Matrouh","Red Sea","Bejaia","Al-, Ismailia","Assiut",
                                               "Behera","Beni Suef","Giza","Luxor","Menia","New Valley","Qena","Suez",
                                               "Suhag"),
                                       ERI = c("Anseba","Debub","Gash Barka","Maekel"),
                                       ETH = c("Addis Ababa","Dire Dawa","Amhara","Gambela","Tigray"),
                                       GAB = c("Estuaire","Moyen-Ogooue","Ngounie","Nyanga","Ogooue-Maritime"),
                                       GHA = c("Northern"),
                                       GIN = c("Faranah","Kankan","Labe","Mamou","Kindia"),
                                       GMB = c("North Bank","West Coast"),
                                       KEN = c("Kakamega","Bungoma","Baringo","Busia","Elgeyo-Marakwet","Kericho","Nakuru","Nandi","Narok","Nyandarua","Uasin Gishu","West Pokot","Siaya","Kisumu"),
                                       LBR = c("Grand Gedeh","Grand Kru","Rivercess","Sinoe","River Gee","Maryland"),
                                       MLI = c("Kayes","Mopti","Ségou"),
                                       MRT = c("Assaba","Gorgol","Hodh El Gharbi","Nouakchott-Sud"),
                                       NER = c("Dosso","Maradi","Tillabéri"),
                                       NGA = c("Adamawa","Benue","Borno","Ebonyi","Ekiti","Enugu",
                                               "Federal Capital Territory","Gombe","Jigawa","Kaduna",
                                               "Kebbi","Kwara","Nassarawa","Niger","Oyo","Plateau",
                                               "Sokoto","Taraba","Yobe"),
                                       SDN = c('Aj Jazirah',"Gedaref","Hala'Ib Triangle","Kassala",
                                               "Northern","Sennar","West Kordofan","White Nile"),
                                       SOM = c("Banadir"),
                                       SSD = c("Eastern Equatoria","Jonglei","Lakes","Northern Bahr El Ghazal",
                                               "Unity","Upper Nile","Warrap"),
                                       STP = c("São Tomé"),
                                       TCD = c("Barh-El-Gazel","Batha","Chari-Baguirmi","Hadjer Lamis",
                                               "Hadjer-Lamis","Lac","Mayo-Kebbi Est","Salamat","Sila",
                                               "Tandjilé","Wadi Fira"),
                                       TGO = c("Kara","Maritime","Centrale","Plateaux"),
                                       TUN = c("Bizerte","Gabes","Kebili","Le Kef","Manouba",
                                               "Nabeul","Mednine","Tataouine","Tunis","Zaghouan"),
                                       TZA = c("Arusha","Dodoma","Katavi","Lindi","Manyara","Mbeya","Morogoro",
                                               "Mtwara","Njombe","Pwani","Rukwa","Ruvuma","Shinyanga","Simiyu",
                                               "Singida","Songwe"),
                                       UGA = c("Central","Eastern","Northern","Western")
                                       )

force1<-c("BFA","BWA","CAF","COM","CPV","ESH","GNB","LBY","LSO","MAR","MDG","MOZ","MWI","NAM","SEN","SWZ","SYC","ZAF","ZMB","ZWE")

force1_admin1<-list(BDI = c("Mairie De Bujumbura","Cibitoke","Bururi"),
                    COD=c("Kongo Central","Haut-Katanga","Tanganyika","Bas-Uélé"))

# Check integrity of chirps data before starting
if(F){
  chirps_check<-rbindlist(lapply(1:length(iso3_choices),function(k){
  iso3_selected<-iso3_choices[k]
  cat("Processing",k,iso3_selected, "...\n")

  file_choice <- grep(iso3_selected, files, value = TRUE)
  dat_raw <- setDT(read_parquet(file_choice))

  chirps<-dat_raw[!is.na(Greenup) & !is.na(Senescence),sum(is.na(rain_greenup_scenescence))/.N]
  data.table(country=iso3_selected,chirps=chirps)
}))
}


# Skip reassigment of seasons in section 3" "National consolidation from scratch on pixel-season summaries"
skip_nat_cons<-list(KEN = T)

for(k in 1:length(iso3_choices)){

  iso3_selected<-iso3_choices[k]
  cat("Processing",k,iso3_selected, "...\n")

  # Load data ----
  {
  file_choice <- grep(iso3_selected, files, value = TRUE)
  dat_raw <- setDT(read_parquet(file_choice))
  }
  # Merge DEM, Aridity, And Coordinates ----
  {
    dat_raw <- px_env[dat_raw, on = .(pixel)]

  # Create combined zone ID
    dat_raw[, zone_id := paste(admin1_name, aridity_bin, elev_bin, sep = "|")]

  # Crate raster temnplate for later plotting
  template_xyz<-unique(merge(dat_raw[,"pixel"],coords_index,all.x=T))[,pixel:=NULL]
  template_r <- terra::rast(template_xyz, type = "xyz", crs = "EPSG:4326")
  terra::values(template_r) <- NA_real_
}
  # Prepare Phenology Data ----
  {
    # Split year and season from flag
    parts <- tstrsplit(dat_raw$flag, "_")
    set(dat_raw, j = "year",   value = parts[[1]])
    set(dat_raw, j = "season", value = parts[[2]])

    dat_raw[,years_pixels:=uniqueN(year)*uniqueN(pixel),by=admin1_name]
    raw_seasons_tab<-suppressWarnings(dcast(dat_raw[!is.na(Greenup)],admin1_name+years_pixels~season))

    # Count season 2 occurrences per pixel
    # Ensure date columns are in  IDate format
    dat_raw[, Greenup := as.IDate(Greenup)]
    dat_raw[, Senescence := as.IDate(Senescence)]
    dat_raw[, DER.sos   := as.IDate(DER.sos)]
    dat_raw[, DER.eos   := as.IDate(DER.eos)]
    dat_raw[, TRS2.sos  := as.IDate(TRS2.sos)]
    dat_raw[, TRS2.eos  := as.IDate(TRS2.eos)]
  }
  # Step 1: Choose Best Available Season Source Per Row ----
  {
      # Priority: Greenup > DER > TRS2 (requires both sos/eos to be present)
      if(use_DER){
        dat_raw[, season_source := fifelse(!is.na(Greenup) & !is.na(Senescence), "greenup",
      fifelse(!is.na(DER.sos) & !is.na(DER.eos), "der",
      fifelse(!is.na(TRS2.sos) & !is.na(TRS2.eos), "trs2", NA_character_)))
      ]

      # unify timing + rainfall metrics based on season_source
        dat_raw[, sos := fifelse(season_source == "greenup", Greenup,
      fifelse(season_source == "der", DER.sos,
      fifelse(season_source == "trs2", TRS2.sos, as.IDate(NA))))
      ]

        dat_raw[, eos := fifelse(season_source == "greenup", Senescence,
      fifelse(season_source == "der", DER.eos,
      fifelse(season_source == "trs2", TRS2.eos, as.IDate(NA))))]

        dat_raw[, rain_total := fifelse(season_source == "greenup", rain_greenup_scenescence,
      fifelse(season_source == "der", rain_der,
      fifelse(season_source == "trs2", rain_trs2, NA_real_)))
      ]

      dat_raw[, rain_total := fifelse(season_source == "greenup", rain_greenup_scenescence,
                                      fifelse(season_source == "der", rain_der,
                                              fifelse(season_source == "trs2", rain_trs2, NA_real_)))
      ]

        dat_raw[, rain_pp28 := fifelse(season_source == "greenup", rain_gs_pp28,
      fifelse(season_source == "der", rain_der_pp28,
      fifelse(season_source == "trs2", rain_trs2_pp28, NA_real_)))
      ]

        dat_raw[, cdd_p45 := fifelse(season_source == "greenup", cdd_greenup_p45,
      fifelse(season_source == "der", cdd_der_p45,
      fifelse(season_source == "trs2", cdd_trs2_p45, NA_real_)))
      ]
      }else{
        dat_raw[, season_source := fifelse(!is.na(Greenup) & !is.na(Senescence), "greenup",NA)]
        dat_raw[, sos := fifelse(season_source == "greenup", Greenup,NA)]
        dat_raw[, eos := fifelse(season_source == "greenup", Senescence,NA)]
        dat_raw[, rain_total := fifelse(season_source == "greenup", rain_gs_pp28,NA)]
        dat_raw[, rain_pp28 := fifelse(season_source == "greenup", rain_greenup_scenescence,NA)]
        dat_raw[, rain_p30 := fifelse(season_source == "greenup", rain_greenup_p30,NA)]
        dat_raw[, cdd_p45 := fifelse(season_source == "greenup", cdd_greenup_p45,NA)]
      }

    dat_raw[,rain_total:=rain_total+rain_pp28]

      if(dat_raw[,all(is.na(rain_total))]){
        cat(iso3_selected,": No CHIRPS rainfall data is available for this location\n")
        no_chirps<-T
      }else{
        no_chirps<-F
      }

      # Set season start / end dates
    dat_raw[, doy :=  yday(sos)]
    dat_raw[, doy_eos :=  yday(eos)]
    dat_raw[, slen :=  as.numeric(eos - sos)]
      # Add the length of the preceding observed season for each pixel.
      # This is based on chronological order of season start date within pixel.
    dat_raw[order(pixel, sos, eos, year, season), preceding_slen := shift(slen), by = pixel]
    }
  # Step 1.2: Relable events close together under the same season ----
  {
    # one row per pixel-year-season for timing comparison
    py_season <- dat_raw[
      !is.na(pixel) & !is.na(year) & !is.na(season) & !is.na(sos),
      .(
        sos = sos[1],
        eos = eos[1],
        doy = doy[1],
        doy_eos = doy_eos[1],
        slen = slen[1]
      ),
      by = .(admin1_name, pixel, year, season)
    ]

    # reshape to wide for comparing s1 and s2 within the same pixel-year
    py_wide <- dcast(
      py_season,
      admin1_name + pixel + year ~ season,
      value.var = c("sos", "eos", "doy", "doy_eos", "slen")
    )

    # parameters for merging close seasons
    min_gap_sos_days <- 60L
    min_gap_eos_to_s2_days <- 15L

    # gap from season-1 SOS to season-2 SOS
    py_wide[, gap_s1s2 := fifelse(
      !is.na(doy_1) & !is.na(doy_2),
      forward_circ_dist(doy_1, doy_2),
      NA_real_
    )]

    # optional extra rule: season 2 should not begin almost immediately after season 1 ends
    py_wide[, gap_eos1_s2 := fifelse(
      !is.na(doy_eos_1) & !is.na(doy_2),
      forward_circ_dist(doy_eos_1, doy_2),
      NA_real_
    )]

    # valid second season only if sufficiently separated
    py_wide[, valid_s2 := fifelse(
      is.na(doy_2),
      FALSE,
      gap_s1s2 >= min_gap_sos_days &
        (is.na(gap_eos1_s2) | gap_eos1_s2 >= min_gap_eos_to_s2_days)
    )]

    # merge validity back
    dat_raw[py_wide, on = .(admin1_name, pixel, year), valid_s2 := i.valid_s2]

    # cleaned season label:
    # - season 1 always kept if detected
    # - season 2 kept only if valid_s2 == TRUE
    dat_raw[, season_clean := fifelse(
      season == 1L, 1L,
      fifelse(season == 2L & valid_s2 == TRUE, 2L, NA_integer_)
    )]

  # Rotated 180  version
  py_season <- dat_raw[
    !is.na(pixel) & !is.na(year) & !is.na(season) & !is.na(sos),
    .(
      sos = sos[1],
      eos = eos[1],
      slen = slen[1]
    ),
    by = .(admin1_name, pixel, year, season)
  ]

  py_season[,sos_180:=sos+180
            ][,eos_180:=eos+180
              ][,doy:=yday(sos_180)
                ][,doy_eos:=yday(eos_180)]

  # reshape to wide for comparing s1 and s2 within the same pixel-year
  py_wide <- dcast(
    py_season,
    admin1_name + pixel + year ~ season,
    value.var = c("sos_180", "eos_180", "doy", "doy_eos", "slen")
  )

  # parameters for merging close seasons
  min_gap_sos_days <- 60L
  min_gap_eos_to_s2_days <- 15L

  # gap from season-1 SOS to season-2 SOS
  py_wide[, gap_s1s2 := fifelse(
    !is.na(doy_1) & !is.na(doy_2),
    forward_circ_dist(doy_1, doy_2),
    NA_real_
  )]

  # optional extra rule: season 2 should not begin almost immediately after season 1 ends
  py_wide[, gap_eos1_s2 := fifelse(
    !is.na(doy_eos_1) & !is.na(doy_2),
    forward_circ_dist(doy_eos_1, doy_2),
    NA_real_
  )]

  # valid second season only if sufficiently separated
  py_wide[, valid_s2 := fifelse(
    is.na(doy_2),
    FALSE,
    gap_s1s2 >= min_gap_sos_days &
      (is.na(gap_eos1_s2) | gap_eos1_s2 >= min_gap_eos_to_s2_days)
  )]

  # merge validity back
  dat_raw[py_wide[,.(admin1_name,pixel,year,valid_s2)], on = .(admin1_name, pixel, year), valid_s2_180 := i.valid_s2]

  # cleaned season label:
  # - season 1 always kept if detected
  # - season 2 kept only if valid_s2 == TRUE
  dat_raw[, season_clean := fifelse(
    season == 1L, 1L,
    fifelse(season == 2L & valid_s2_180 == TRUE, 2L, NA_integer_)
  )]

  # if you want to keep rows with invalid season 2 in the table but treat them as absent,
  # leave them as NA. If instead you want everything forced into season 1 in mono-season systems,
  # do that later, after summaries.

  # optional diagnostics
  dat_raw[, invalid_s2_close := season == 2L & valid_s2 == FALSE]
    }

  # Step 1.3: Harmonize season labels within each admin1 using original split_seasons_cpt ----

  dat<-dat[rain_total>10]
  dat<-copy(dat_raw)
  dat[!is.na(season_clean), n_seasons_clean := uniqueN(season_clean), by = admin1_name]

  if(iso3_selected %in% force1){
    dat[,season_harmonized:=1]
  }else{
      admin_skip<-iso3_selected %in% names(skip_season_harmonization_admin1)

      if(admin_skip){
        rm_admin<-unlist(skip_season_harmonization_admin1[iso3_selected])
        dat1<-dat[admin1_name %in% rm_admin]
        dat1[,season_harmonized:=season_clean]
        dat<-dat[!admin1_name %in% rm_admin]
      }

      #dat[n_seasons_clean>1 & !is.na(sos) & !admin1_name %in% admin_3_season[[iso3_selected]], season_harmonized :=  split_seasons_cpt(x = sos, max_seasons = 2L), by = admin1_name]
      #dat[n_seasons_clean>1 & !is.na(sos) & admin1_name %in% admin_3_season[[iso3_selected]], season_harmonized :=  split_seasons_cpt(x = sos, max_seasons = 3L), by = admin1_name]

      dat[n_seasons_clean>1 & !is.na(sos) & !admin1_name %in% admin_3_season[[iso3_selected]],
          season_harmonized :=assign_admin_seasons(sos, eos, max_seasons = 2L)$season, by = admin1_name]

      dat[n_seasons_clean>1 & !is.na(sos) & admin1_name %in% admin_3_season[[iso3_selected]],
          season_harmonized :=assign_admin_seasons(sos, eos, max_seasons = 3L)$season, by = admin1_name]


      dat[n_seasons_clean==1 & !is.na(sos) & !admin1_name %in% admin_3_season, season_harmonized :=  1]

      if(admin_skip){
        dat<-rbind(dat,dat1)
      }
  }

  if(iso3_selected %in% force1){
    dat[,season_harmonized:=1]
  }

  if(iso3_selected %in% names(force1_admin1)){
    dat[admin1_name %in% unlist(force1_admin1[iso3_selected]),season_harmonized:=1]
  }

  # Exceptions
 # skip_admin<-skip_season_close[[iso3_selected]]
  #dat[admin1_name %in% skip_admin,season_harmonized:=season]

  # Step 1.4: Plot histograms to validate season cleaning ----
  {
    print(raw_seasons_tab)

    bin_width <- 10
    breaks <- seq(0.5, 360.5, by = bin_width)
    mids   <- breaks[-length(breaks)] + bin_width / 2

    polar_dat <- copy(
      dat[!is.na(season_harmonized) & !is.na(doy),
          .(admin1_name, season_harmonized, doy)]
    )

    # wrap DOY cleanly to 1..365
    polar_dat[, doy_wrap := ((doy - 1) %% 365) + 1]

    # assign bins
    polar_dat[, bin_id := cut(
      doy_wrap,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE,
      labels = FALSE
    )]

    polar_dat[, bin_mid := mids[bin_id]]

    # count by admin1 x season x bin
    polar_counts <- polar_dat[
      ,
      .(count = .N),
      by = .(admin1_name, season_harmonized, bin_mid)
    ]

    # complete missing bins
    polar_template <- CJ(
      admin1_name = unique(polar_dat$admin1_name),
      season_harmonized = unique(polar_dat$season_harmonized),
      bin_mid = mids,
      unique = TRUE
    )

    polar_counts <- polar_template[
      polar_counts,
      on = .(admin1_name, season_harmonized, bin_mid)
    ][
      is.na(count), count := 0
    ]

    # density within each admin1 x season
    polar_counts[, density := count / sum(count), by = .(admin1_name, season_harmonized)]

    (g_polar <- ggplot(
      polar_counts,
      aes(x = bin_mid, y = count, fill = factor(season_harmonized))
    ) +
      geom_col(
        position = "identity",
        alpha = 0.45,
        width = bin_width
      ) +
      coord_polar(start = -pi/2, clip = "off") +
      facet_wrap(~admin1_name, scales = "free") +
      scale_x_continuous(
        limits = c(0.5, 360.5),
        breaks = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
        labels = month.abb,
        expand = c(0, 0)
      ) +
      labs(
        title = paste0(iso3_selected, " - Polar seasonal detections"),
        x = NULL,
        y = "Density",
        fill = "Season"
      ) +
      theme_bw())


  }
  # Step 2: Gates And QC Flags (Zone-Based, Aridity-Aware)
  # Step 2.0: Long season flags ----
  {
    # Calculate the proportion of season 1 vs season 2 occurrences for each pixel across years
    dat[,season_prop:=.N/years_pixels,by=.(admin1_name,season_harmonized)]
    dat[order(pixel, sos, eos, year, season), preceding_season_prop := shift(slen), by = pixel]

    # Add flags for long-preceding season which could push the onset of the following season to an odd time.
    dat[preceding_season_prop>=0.2, s_long_preceding := fifelse(preceding_slen > 30*8, TRUE, FALSE)]
    dat[preceding_season_prop<0.2, s_long_preceding := fifelse(preceding_slen > 30*11, TRUE, FALSE)]
    dat[preceding_season_prop>=0.2, s_long := fifelse(slen > 30*8, TRUE, FALSE)]
    dat[preceding_season_prop<0.2, s_long := fifelse(slen > 30*11, TRUE, FALSE)]
  }
  # Step 2.1: Zone-Specific Thresholds (zone_id + season_harmonized) ----
  {
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

    params$zone_thresholds<-zone_thr

    # Join thresholds back
    dat[zone_thr, on = .(zone_id, season_harmonized), `:=`(
      p30_fail_thr  = i.p30_fail_thr,
      p30_weak_thr  = i.p30_weak_thr,
      cdd_fail_thr  = i.cdd_fail_thr,
      cdd_weak_thr  = i.cdd_weak_thr,
      rtot_fail_thr = i.rtot_fail_thr,
      rtot_weak_thr = i.rtot_weak_thr
    )]
  }
  # Step 2.2: NDVI Reliability Indicator ----
  {
    dat[, ndvi_ok := {
      has_fit <- !is.na(NSE) & !is.na(R2)
      fifelse(has_fit, (NSE >= ndvi_nse_min & R2 >= ndvi_r2_min), TRUE)
    }]
  }
  # Step 2.3: Quantile-Based Flags ----
  {
    dat[, `:=`(
    q_p30_fail  = !is.na(rain_p30)   & !is.na(p30_fail_thr)  & (rain_p30   <= p30_fail_thr),
    q_p30_weak  = !is.na(rain_p30)   & !is.na(p30_weak_thr)  & (rain_p30   <= p30_weak_thr),

    q_cdd_fail  = !is.na(cdd_p45)    & !is.na(cdd_fail_thr)  & (cdd_p45    >= cdd_fail_thr),
    q_cdd_weak  = !is.na(cdd_p45)    & !is.na(cdd_weak_thr)  & (cdd_p45    >= cdd_weak_thr),

    q_rtot_fail = !is.na(rain_total) & !is.na(rtot_fail_thr) & (rain_total <= rtot_fail_thr),
    q_rtot_weak = !is.na(rain_total) & !is.na(rtot_weak_thr) & (rain_total <= rtot_weak_thr)
    )]
  }
  # Step 2.4: Hard-Limit Flags (Aridity-Aware) ----
  {
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
  }
  # Step 2.5: Combine Gates (CDD Can Fail On Its Own) ----
  {
    dat[, gate_moist_fail := (q_cdd_fail | hard_cdd_fail) |((q_p30_fail | hard_p30_fail) & (q_rtot_fail | hard_rtot_fail))]

    dat[, gate_moist_weak := (q_p30_weak | hard_p30_weak) |(q_cdd_weak | hard_cdd_weak) |(q_rtot_weak | hard_rtot_weak)]
  }
  # Step 2.6: Source-Aware Final QC ----
  {
    dat[, `:=`(season_fail = FALSE, season_weak = FALSE)]

    dat[season_source == "greenup", `:=`(
      season_fail = gate_moist_fail & !ndvi_ok,
      season_weak = gate_moist_weak | (gate_moist_fail & ndvi_ok)
    )]

    dat[season_source %chin% c("der","trs2"), `:=`(
      season_fail = gate_moist_fail,
      season_weak = gate_moist_weak & !season_fail
    )]
  }
  # Step 3: Labels ----
  {
    n_years<-length(unique(dat$year))

    dat[, season_qc := fifelse(season_fail, "fail",
                      fifelse(season_weak, "weak", "ok"))]


    # admin1 doy start of season averages ----
    dat[,pixel_n:=length(unique(pixel))*length(unique(year)),by=.(admin1_name)]
    dat[,no_detect:=n_years-.N,by=.(admin1_name,pixel,season_harmonized)
    ][no_detect<0, no_detect:=0]
    dat[,failed:=sum(season_qc=="fail",na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]
    dat[,weak:=sum(season_qc=="weak",na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]
    dat[,ndvi_weak:=sum(!ndvi_ok,na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]
    dat[,s_long_p:=sum(s_long_preceding,na.rm=TRUE),by=.(admin1_name,pixel,season_harmonized)]

    # Refine order of seasons ----

    # Across all admin1 units on average which season comes first if >1 season.
    season_order <- dat[!is.na(season_harmonized)
      , .(median_sos = quantile_circular_safe(doy, probs = 0.5)),
      by = season_harmonized
    ][order(median_sos), season_harmonized]

    # Create mapping: old season -> new ordered index
    map <- setNames(seq_along(season_order), season_order)

    # Apply mapping
    dat[, season_harmonized := as.integer(map[as.character(season_harmonized)])]

    # National consolidation from scratch on pixel-season summaries ----
    if(!is.null(skip_nat_cons[[iso3_selected]])){
      if(!skip_nat_cons[[iso3_selected]]){
    # summarise one object per pixel x current season label
    pix_season <- dat[
      !is.na(pixel) &
        !is.na(season_harmonized) &
        !is.na(doy) &
        !is.na(doy_eos),
      .(
        pix_sos = quantile_circular_safe(doy, probs = 0.5),
        pix_eos = quantile_circular_safe(doy_eos, probs = 0.5),
        n_obs   = .N
      ),
      by = .(pixel, season_old = season_harmonized)
    ]

    # decide k = number of national seasons actually present
    k <- pix_season[, uniqueN(season_old)]
    k <- min(max(k, 1L), 3L)

    # build circular feature space
    pix_season[
      ,
      `:=`(
        sos_cos = cos(2 * pi * pix_sos / 365),
        sos_sin = sin(2 * pi * pix_sos / 365),
        eos_cos = cos(2 * pi * pix_eos / 365),
        eos_sin = sin(2 * pi * pix_eos / 365)
      )
    ]

    feat <- as.matrix(
      pix_season[, .(sos_cos, sos_sin, eos_cos, eos_sin)]
    )

    # cluster nationally from scratch
    set.seed(1)
    km <- kmeans(feat, centers = k, nstart = 50)

    pix_season[, cluster_raw := km$cluster]

    # reorder clusters by national median SOS so labels become 1,2,3 in timing order
    cluster_order <- pix_season[
      ,
      .(median_sos = quantile_circular_safe(pix_sos, probs = 0.5)),
      by = cluster_raw
    ][
      order(median_sos)
    ][
      ,
      season_new := seq_len(.N)
    ]

    pix_season <- cluster_order[
      pix_season,
      on = .(cluster_raw)
    ]

    # write national reassigned season labels back to all rows
    dat <- dat[
      pix_season[, .(pixel, season_old, season_new)],
      on = .(pixel, season_harmonized = season_old)
    ]

    dat[!is.na(season_new), season_harmonized := season_new]

    # cleanup
    dat[
      ,
      c(
        "pix_sos", "pix_eos", "n_obs",
        "sos_cos", "sos_sin", "eos_cos", "eos_sin",
        "cluster_raw", "median_sos", "season_new"
      ) := NULL
    ]

    }
    }
    # Apply any overrides ----
    # dat <- apply_season_overrides(dat, season_overrides, iso3_selected)

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

    params$dat$admin_pixel_year<-dat
    params$dat$admin<-admin1_tab

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

    params$dat_c$admin_pixel_year<-dat_c
    params$dat_c$admin<-admin1_tab_c


    # Join into a list
    admin1_tabs <- list(
      all = admin1_tab,
      conservative = admin1_tab_c
    )
  }
  # Plot differences between all and conservative summaries ----
  {
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

    g<-ggplot(
      compare_data,
      aes(
        x = metric,
        y = q50_plot,
        colour = dataset,
        shape = dataset
      )) +
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
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    params$plots$all_vs_cons<-g

    plot(g)

  }
  # Label wetter/drier, longer/shorter, first/second/third by admin1 ----
  {
    # Labeling thresholds (to avoid over-interpreting tiny differences)
    label_diff_pct_rain   <- 0
    label_diff_pct_length <- 0

    admin1_tabs <- lapply(admin1_tabs, function(x) {
      x[, `:=`(

        # wetter / drier label: only the wettest and driest are labelled
        rain_label = {
          r <- rain_mean

          if (.N < 2L || all(is.na(r))) {
            rep(NA_character_, .N)

          } else {
            r_non_na <- r[!is.na(r)]

            if (length(r_non_na) < 2L) {
              rep(NA_character_, .N)

            } else {
              rmax <- max(r_non_na)
              rmin <- min(r_non_na)

              rel_diff <- if (is.finite(rmax) && rmax > 0) (rmax - rmin) / rmax else NA_real_

              if (!is.na(rel_diff) && rel_diff < label_diff_pct_rain) {
                rep("similar", .N)

              } else {
                lab <- rep(NA_character_, .N)

                wet_idx <- which(r == rmax)
                dry_idx <- which(r == rmin)

                # only assign if unique; otherwise ties become "similar"
                if (length(wet_idx) == 1L) lab[wet_idx] <- "wettest"
                if (length(dry_idx) == 1L) lab[dry_idx] <- "driest"

                # if both extremes tied away, optionally mark all as similar
                if (all(is.na(lab))) lab <- rep("similar", .N)

                lab
              }
            }
          }
        },

        # longer / shorter label: only the longest and shortest are labelled
        length_label = {
          sl <- slen_mean

          if (.N < 2L || all(is.na(sl))) {
            rep(NA_character_, .N)

          } else {
            sl_non_na <- sl[!is.na(sl)]

            if (length(sl_non_na) < 2L) {
              rep(NA_character_, .N)

            } else {
              slmax <- max(sl_non_na)
              slmin <- min(sl_non_na)

              rel_diff <- if (is.finite(slmax) && slmax > 0) (slmax - slmin) / slmax else NA_real_

              if (!is.na(rel_diff) && rel_diff < label_diff_pct_length) {
                rep("similar", .N)

              } else {
                lab <- rep(NA_character_, .N)

                long_idx  <- which(sl == slmax)
                short_idx <- which(sl == slmin)

                if (length(long_idx)  == 1L) lab[long_idx]  <- "longest"
                if (length(short_idx) == 1L) lab[short_idx] <- "shortest"

                if (all(is.na(lab))) lab <- rep("similar", .N)

                lab
              }
            }
          }
        },

        # timing label: assign first / second / third by increasing SOS
        timing_label = {
          timing_names <- c("first", "second", "third", "fourth", "fifth")
          lab <- rep(NA_character_, .N)
          ok <- !is.na(season_harmonized)
          lab[ok] <- timing_names[season_harmonized[ok]]
          lab
        }

      ), by = admin1_name]

      x
    })

    names(admin1_tabs) <- c("all", "conservative")
  }
  # Show differences between seasons in terms of rain and season length ----
  {
    if(dat[!is.na(season_harmonized),uniqueN(season_harmonized)==2]){
      season_diff_base <- dat |>
      subset(
          !is.na(season_harmonized) &
          !is.na(year) &
          !is.na(pixel) &
  #        !is.na(rain_greenup_scenescence) &
          !is.na(Greenup) &
          !is.na(Senescence)
      ) |>
      transform(
        year = as.integer(year),
        season_length = as.numeric(Senescence - Greenup)
      ) |>
      setDT()
    }
  }
  # Admin1-year means: one mean value per season per admin1-year
  {
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
    #  !is.na(rain_1) & !is.na(rain_2) &
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

    # Pixel-year means: one mean value per season per pixel-year
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

    if(!no_chirps){
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

   g<-gridExtra::grid.arrange(
      p_rain_abs,
      p_rain_norm,
      ncol = 1
    )

   params$plots$diff_rain<-g

   plot(g)
    }

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
        ci_low = mean_diff - 1.96 *  sd_diff / sqrt(n),
        ci_high = mean_diff + 1.96 * sd_diff / sqrt(n)
      )
    ]

    slen_levels <- slen_plot_dat[
      order(mean_diff),
      unique(admin1_name)
    ]

    slen_plot_dat[
      , admin1_name := factor(admin1_name, levels = slen_levels)
    ]

    g<-ggplot(slen_plot_dat, aes(x = admin1_name, y = mean_diff)) +
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

    params$plots$diff_slen<-g
    plot(g)
  }
  # Create pixel-level maps ----
  {
    focal_dat  <- copy(dat_c)
    focal_admin <- copy(admin1_tab_c)

    # Remove % from summary column names for safer joins and downstream handling
    setnames(
      focal_admin,
      old = c("sos_10%", "sos_50%", "sos_90%", "eos_10%", "eos_50%", "eos_90%"),
      new = c("sos_10",  "sos_50",  "sos_90",  "eos_10",  "eos_50",  "eos_90"),
      skip_absent = TRUE
    )

    # Keep only the admin1-level season timing columns needed for matching
    focal_admin_join <- unique(
      focal_admin[
        ,
        .(admin1_name, season_harmonized, sos_10, sos_50, sos_90, eos_10, eos_50, eos_90)
      ]
    )

    # ------------------------------------------------------------------
    # Mismatch detection for 1 to 3 seasons
    # ------------------------------------------------------------------

    # Attach each row to its OWN admin1-season reference window
    own_ref <- copy(focal_admin_join)
    setnames(
      own_ref,
      old = c("sos_10", "sos_50", "sos_90", "eos_10", "eos_50", "eos_90"),
      new = c("own_sos_10", "own_sos_50", "own_sos_90",
              "own_eos_10", "own_eos_50", "own_eos_90")
    )

    focal_dat <- own_ref[focal_dat, on = .(admin1_name, season_harmonized)]

    # Row id
    focal_dat[, row_id := .I]

    # Own-window match
    focal_dat[
      ,
      `:=`(
        own_sos_match = fifelse(
          !is.na(doy) & !is.na(own_sos_10) & !is.na(own_sos_90),
          in_circular_window(doy, own_sos_10, own_sos_90),
          NA
        ),
        own_eos_match = fifelse(
          !is.na(doy_eos) & !is.na(own_eos_10) & !is.na(own_eos_90),
          in_circular_window(doy_eos, own_eos_10, own_eos_90),
          NA
        )
      )
    ]

    # Other-season reference table with explicit names
    other_ref <- copy(focal_admin_join)
    setnames(
      other_ref,
      old = c("season_harmonized", "sos_10", "sos_50", "sos_90", "eos_10", "eos_50", "eos_90"),
      new = c("other_season",      "other_sos_10", "other_sos_50", "other_sos_90",
              "other_eos_10",      "other_eos_50", "other_eos_90")
    )

    # Compare each row against all OTHER seasons in the same admin1
    other_matches <- focal_dat[
      other_ref,
      on = .(admin1_name),
      allow.cartesian = TRUE,
      nomatch = 0L
    ][
      season_harmonized != other_season
    ][
      ,
      .(
        other_sos_match = any(
          !is.na(other_sos_10) & !is.na(other_sos_90) &
            in_circular_window(doy, other_sos_10, other_sos_90),
          na.rm = TRUE
        ),
        other_eos_match = any(
          !is.na(other_eos_10) & !is.na(other_eos_90) &
            in_circular_window(doy_eos, other_eos_10, other_eos_90),
          na.rm = TRUE
        )
      ),
      by = row_id
    ]

    # Join back
    focal_dat[
      other_matches,
      on = .(row_id),
      `:=`(
        other_sos_match = i.other_sos_match,
        other_eos_match = i.other_eos_match
      )
    ]

    # No other season -> FALSE
    focal_dat[
      ,
      `:=`(
        other_sos_match = fifelse(is.na(other_sos_match), FALSE, other_sos_match),
        other_eos_match = fifelse(is.na(other_eos_match), FALSE, other_eos_match)
      )
    ]

    # Final mismatch flags
    focal_dat[
      ,
      `:=`(
        sos_mismatch = !is.na(own_sos_match) & !own_sos_match & other_sos_match,
        eos_mismatch = !is.na(own_eos_match) & !own_eos_match & other_eos_match
      )
    ]

    # Remove helper columns used only for mismatch detection
    focal_dat[
      ,
      c(
        "row_id",
        "own_sos_10", "own_sos_50", "own_sos_90",
        "own_eos_10", "own_eos_50", "own_eos_90",
        "own_sos_match", "other_sos_match",
        "own_eos_match", "other_eos_match"
      ) := NULL
    ]

    # Remove rows with a mismatch, but keep a count per admin1/pixel/season
    focal_dat <- focal_dat[
      ,
      mismatch := sum(sos_mismatch | eos_mismatch),
      by = .(admin1_name, pixel, season_harmonized)
    ][
      !(sos_mismatch | eos_mismatch)
    ]

    # ------------------------------------------------------------------
    # Histograms of greenup-based SOS by admin1 and season
    # ------------------------------------------------------------------
    {
      bin_width <- 10
      breaks <- seq(0.5, 360.5, by = bin_width)
      mids   <- breaks[-length(breaks)] + bin_width / 2

      polar_dat <- copy(
        focal_dat[!is.na(season_harmonized) & !is.na(doy),
            .(admin1_name, season_harmonized, doy)]
      )

      # wrap DOY cleanly to 1..365
      polar_dat[, doy_wrap := ((doy - 1) %% 365) + 1]

      # assign bins
      polar_dat[, bin_id := cut(
        doy_wrap,
        breaks = breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
      )]

      polar_dat[, bin_mid := mids[bin_id]]

      # count by admin1 x season x bin
      polar_counts <- polar_dat[
        ,
        .(count = .N),
        by = .(admin1_name, season_harmonized, bin_mid)
      ]

      # complete missing bins
      polar_template <- CJ(
        admin1_name = unique(polar_dat$admin1_name),
        season_harmonized = unique(polar_dat$season_harmonized),
        bin_mid = mids,
        unique = TRUE
      )

      polar_counts <- polar_template[
        polar_counts,
        on = .(admin1_name, season_harmonized, bin_mid)
      ][
        is.na(count), count := 0
      ]

      # density within each admin1 x season
      polar_counts[, density := count / sum(count), by = .(admin1_name, season_harmonized)]

      (g_polar <- ggplot(
        polar_counts,
        aes(x = bin_mid, y = count, fill = factor(season_harmonized))
      ) +
          geom_col(
            position = "identity",
            alpha = 0.45,
            width = bin_width
          ) +
          coord_polar(start = -pi/2, clip = "off") +
          facet_wrap(~admin1_name, scales = "free") +
          scale_x_continuous(
            limits = c(0.5, 360.5),
            breaks = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349),
            labels = month.abb,
            expand = c(0, 0)
          ) +
          labs(
            title = paste0(iso3_selected, " - Polar seasonal detections"),
            x = NULL,
            y = "Density",
            fill = "Season"
          ) +
          theme_bw())

      plot(g_polar)
    }

    # ------------------------------------------------------------------
    # Average focal_dat over years including quantiles, by admin1/pixel/season
    # ------------------------------------------------------------------

    n_years <- length(unique(focal_dat$year))

    focal_dat_avg <- focal_dat[
      !is.na(season_harmonized) &
        !is.na(doy) &
        !is.na(doy_eos),
      {
        sos_ci  <- quantile_circular_safe(doy, probs = c(0.1, 0.5, 0.9))
        names(sos_ci) <- paste0("sos_", names(sos_ci))

        eos_ci  <- quantile_circular_safe(doy_eos, probs = c(0.1, 0.5, 0.9))
        names(eos_ci) <- paste0("eos_", names(eos_ci))

        rain_ci <- quantile(rain_total, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
        names(rain_ci) <- paste0("rain_", names(rain_ci))

        slen_ci <- quantile(slen, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
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
          n = .N,
          n_prop = .N / n_years,
          failed_seasons = mean(failed, na.rm = TRUE),
          weak_seasons = sum(season_qc == "weak", na.rm = TRUE),
          mismatch = mean(mismatch, na.rm = TRUE),
          no_detect = mean(no_detect, na.rm = TRUE),
          s_long_p = mean(s_long_p, na.rm = TRUE),
          ndvi_weak = mean(ndvi_weak, na.rm = TRUE)
        )
      },
      by = .(admin1_name, pixel, season_harmonized)
    ]

    focal_dat_avg[
      ,
      failed_prop := failed_seasons / (n_years - mismatch - no_detect - s_long_p)
    ][
      ,
      weak_prop := weak_seasons / (n_years - mismatch - no_detect - s_long_p)
    ][
      ,
      mismatch_prop := mismatch / (n_years - failed_seasons - no_detect - s_long_p)
    ][
      ,
      no_detect_prop := no_detect / n_years
    ][
      ,
      s_long_p_prop := s_long_p / n_years
    ][
      ,
      ndvi_weak_prop := ndvi_weak / (n_years - mismatch - no_detect - s_long_p)
    ]

    # Merge lat-lon from pixel index
    focal_dat_avg <- merge(
      focal_dat_avg,
      coords_index,
      by = "pixel",
      all.x = TRUE
    )

    if (no_chirps) {
      focal_dat_avg[, c("rain_10", "rain_50", "rain_90", "failed_prop", "weak_prop") := NULL]
    }

    # ------------------------------------------------------------------
    # Make variables long (sos,eos,rain,slen) for easier plotting
    # ------------------------------------------------------------------

    focal_dat_avg_long <- melt(
      focal_dat_avg[
        ,
        !c("failed_seasons", "weak_seasons", "mismatch", "no_detect", "ndvi_weak", "s_long_p")
      ],
      id.vars = c("admin1_name", "pixel", "season_harmonized", "x", "y")
    )

    focal_dat_avg_long[
      ,
      variable := factor(variable)
    ]

    params$dat_final$pixel_avg <- focal_dat_avg_long

    # ------------------------------------------------------------------
    # Make a raster stack for each season where the layers are variables
    # ------------------------------------------------------------------

    raster_list <- list()

    for (season in sort(unique(na.omit(focal_dat_avg_long$season_harmonized)))) {

      season_data <- focal_dat_avg_long[season_harmonized == season]
      var_levels  <- unique(as.character(season_data$variable))

      r_list <- lapply(var_levels, function(v) {
        v_dat <- season_data[variable == v, .(x, y, value)]
        r_v <- template_r
        cell_ids <- terra::cellFromXY(r_v, as.matrix(v_dat[, .(x, y)]))
        vals <- terra::values(r_v, mat = FALSE)
        vals[cell_ids] <- v_dat$value
        terra::values(r_v) <- vals

        names(r_v) <- paste0("season", season, "_", v)
        r_v
      })

      raster_list[[as.character(season)]] <- terra::rast(r_list)
    }

    names(raster_list) <- paste("season", names(raster_list), sep = "_")

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

  plot_season_stack(raster_list$season_1, season_name = NULL, circ_palette = "phenology", seq_palette = "magma")
  plot_season_stack(raster_list$season_2, season_name = NULL, circ_palette = "phenology", seq_palette = "magma")

  plot_season_vars(raster_list, vars = "sos")

  n_prop_min <- 0.1

  raster_list_masked <- mask_raster_list_by_nprop(
    raster_list = raster_list,
    n_prop_min = n_prop_min
  )

  plot_season_vars(raster_list_masked, vars = "sos")

  # Save outputs ----
  {
    for (season_name in names(raster_list)) {
      r_stack <- raster_list[[season_name]]

      for (i in seq_len(terra::nlyr(r_stack))) {
        layer_name <- names(r_stack)[i]
        cat("Saving", iso3_selected, layer_name, "...\r")
        r_layer <- r_stack[[i]]
        output_path <- file.path(output_dir_rast, paste0(iso3_selected, "_", layer_name, ".tif"))
        terra::writeRaster(r_layer, output_path, overwrite = TRUE)
      }
    }
  }
    save(params, file = file.path(output_dir,"parameters_data_plots.RData"))

}
