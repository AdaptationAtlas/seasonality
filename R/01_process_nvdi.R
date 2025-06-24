# 1) NDVI Green-Up Detection Workflow for Africa ####
# Project: Africa Agriculture Adaptation Atlas (AAAA)
# PI: Pete Steward | p.steward@cgiar.org | ORCID: 0000-0003-3985-4911
# Organization: Alliance of Bioversity International & CIAT

#' @title NDVI-Based Green-Up Analysis Workflow
#' @description This script prepares and analyzes 8-day NDVI data (GLASS) to detect green-up events,
#' assign seasonal timing (including possible bimodal seasonality), and estimate onset, peak, and greendown.
#' It follows the phenological detection principles used in literature (Jönsson & Eklundh, Zhang et al., etc.).
#'
#' References:
#' - Jönsson, P., & Eklundh, L. (2004). TIMESAT—A program for analyzing time-series of satellite sensor data. Computers & Geosciences. [DOI: 10.1016/j.cageo.2004.05.006]
#' - Zhang, X., et al. (2003). Monitoring vegetation phenology using MODIS. Remote Sensing of Environment. [DOI: 10.1016/S0034-4257(03)00084-6]
#' - Chen, J., et al. (2023). GLASS NDVI and EVI V1.0 Dataset Description. [DOI: 10.6084/m9.figshare.22220050]

# Libraries ####
#' @import data.table terra future.apply progressr phenofit

## 1.1) Loading Required Packages ####

pacman::p_load(
  arrow,
  dbscan,
  data.table,    # fast data manipulation
  terra,         # raster operations
  fs,            # filesystem operations
  progressr,     # progress bar integration
  future,        # parallel processing backend
  future.apply,  # parallelized apply functions
  pbapply,
  glue,          # string interpolation
  phenofit,      # smoothing and phenology algorithms
  stringr,        # string manipulation
  miceadds
)

# phenofit
# DOI: 10.1002/joc.6081
source("R/01_calc_pixel_greenup.R")

# 2) Parameters ####
input_dir <- dirs$glass_ndvi_tif
outdir <- dirs$nvdi_phenology
my_tempdir<-file.path(outdir,"temp")
if(!dir.exists(my_tempdir)){
  dir.create(my_tempdir,recursive =T)
}

# 3) Process NDVI Stack ####

clean_ndvi <- function(file, out_dir, crop_extent = NULL, mask = NULL) {
  r <- terra::rast(file)+0

  # Clip out invalid values
  r[r < 0 | r > 1] <- NA

  # Crop and mask if needed
  if (!is.null(crop_extent)) r <- terra::crop(r, crop_extent)
  if (!is.null(mask))        r <- terra::mask(r, mask)

  # Write output
  terra::writeRaster(r, file, overwrite = TRUE,
              datatype = "FLT4S",
              gdal = c("COMPRESS=ZSTD",
                       "PREDICTOR=2",
                       "TILED=YES",
                       "BIGTIFF=YES",
                       "COPY_SRC_OVERVIEWS=YES",
                       "BLOCKXSIZE=512",
                       "BLOCKYSIZE=512",
                       "NUM_THREADS=ALL_CPUS",
                       "COG=YES"))
}

tif_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

invisible(
  pbapply::pblapply(1:length(tif_files),FUN = function(i) {
          clean_ndvi(file=tif_files[i], output_dir,crop_extent=af_v)
    })
)

r_stack<-rast(tif_files)
dates<-gsub("A","",unlist(tstrsplit(basename(tif_files),"[.]",keep=3)))
years<-as.numeric(substr(dates,1,4))
days<-as.numeric(substr(dates,5,7))
dates<-as.Date(paste0(years, "-01-01")) + days - 1

# Image per year
nptperyear<-sum(grepl("2021",years))

# 4) Apply NDVI Smoothing & Green-Up Extraction Function to Raster Stack ####

# Convert raster to matrix
vals <- terra::as.matrix(r_stack)  # ncell × nlayer
coords <- terra::xyFromCell(r_stack, 1:ncell(r_stack))

# Save coords
coords_index<-data.table(coords)
coords_index[,pixel:=.I]
write_parquet(coords_index,file.path(dirs$nvdi_phenology,"pixel_index.parquet"))

good_cells<-app(r_stack,sum)
good_cells<-which(!is.na(good_cells[]))

vals_good <- vals[good_cells, , drop = FALSE]   # (n_good  ×  n_time)
rm(vals)                                        # free RAM
gc()

vals_good <- asplit(vals_good, 1)   # 1 → split over rows
for(i in 1:length(vals_good)){
  cat(i,"/",length(vals_good),"      \r")
  vals_good[[i]]<-list(vals=vals_good[[i]],index=good_cells[i])
}

out_files<-file.path(my_tempdir,paste0(good_cells,".RData"))
to_do<-!file.exists(out_files)

if(length(to_do)>0){
vals_good<-vals_good[to_do]

# Dev note: This has an irritating memory leak, I had to restart many times, needs revising

future::plan(multisession, workers = parallel::detectCores()-2)

handlers(global = TRUE)                      # enable handlers once
handlers("progress")                         # progress bar in console

with_progress({
  p <- progressor(along = vals_good)        # initialise bar

 future_lapply(vals_good,function(dat) {
      pixel <- dat$index
      file <- file.path(my_tempdir, paste0(pixel, ".RData"))
      if (!file.exists(file)) {
        res <- calc_pixel_greenup(
          ndvi_ts        = as.numeric(dat$vals),
          dates          = dates,
          nptperyear     = nptperyear,
          w_low          = 0.1,
          methods_to_fit = c("Beck", "Elmore", "Zhang"),
          min_frac       = 0.20,
          min_years      = 3
        )
        if (!is.null(res)) {
          res$pixel <- pixel
        }
        save(res, file = file)
        rm(res,pixel,dat,file)
      }
      p()
      return(NULL)
    })
  })

plan(sequential)
}

analysis_params<-list(fun=calc_pixel_greenup,
                      nptperyear     = nptperyear,
                      w_low          = 0.1,
                      methods_to_fit = c("Beck", "Elmore", "Zhang"),
                      min_frac       = 0.20,
                      min_years      = 3)

jsonlite::write_json(analysis_params,file.path(dirs$nvdi_phenology,"phenofit_analysis_params.json"))

meta <- tibble::tribble(
  ~column,       ~type,       ~description,
  "flag",        "chr",       "Season ID: `YYYY_n` where n = 1 (1st) or 2 (2nd) season in that calendar year",
  "origin",      "Date",      "Origin date supplied to check_input() – here ‘YYYY-01-01’",
  "TRS2.sos",    "Date",      "Start-of-season when fitted curve crosses 20 % of seasonal amplitude",
  "TRS2.eos",    "Date",      "End-of-season for 20 % threshold",
  "TRS5.sos",    "Date",      "Start-of-season (50 % threshold)",
  "TRS5.eos",    "Date",      "End-of-season (50 % threshold)",
  "DER.sos",     "Date",      "Date of maximum positive derivative (fastest green-up)",
  "DER.pos",     "Date",      "Date of seasonal peak (maximum NDVI)",
  "DER.eos",     "Date",      "Date of maximum negative derivative (fastest senescence)",
  "UD",          "Date",      "‘Upturn’ – logistic breakpoint before rapid green-up",
  "SD",          "Date",      "‘Stable’ – plateau after green-up / before senescence",
  "DD",          "Date",      "‘Downturn’ – logistic breakpoint where rapid decline begins",
  "RD",          "Date",      "‘Recovery/Dormancy’ – logistic breakpoint after decline",
  "Greenup",     "Date",      "Four-stage logistic: date of green-up (Zhang model only)",
  "Maturity",    "Date",      "Four-stage logistic: date of maturity (plateau)",
  "Senescence",  "Date",      "Four-stage logistic: start of senescence",
  "Dormancy",    "Date",      "Four-stage logistic: dormancy date",
  "meth",        "chr",       "Curve-fitting method selected as ‘best’ (Beck, Elmore, Zhang…)",
  "R2",          "num",       "Coefficient of determination for selected model",
  "NSE",         "num",       "Nash–Sutcliffe Efficiency",
  "R",           "num",       "Pearson correlation between fitted and observed NDVI",
  "RMSE",        "num",       "Root-mean-square error of fit (lower = better)",
  "pvalue",      "num",       "p-value of the correlation (R)",
  "pixel",       "int",       "Raster cell index (row number in original `r_stack`)"
)


fwrite(meta,file.path(dirs$nvdi_phenology,"glass_pheno_raw_metadata.csv"))


# 6) Clustering Green-Up into Seasons ####
# This could use k-means or DBSCAN to group events into season 1, 2, etc.
# Also consider alternate clustering by peak rainfall, season duration, or greendown interval
# Future implementation placeholder

pheno_file<-file.path(dirs$nvdi_phenology,"glass_nvdi_phenofit.parquet")
coords_index<-read_parquet(file.path(dirs$nvdi_phenology,"pixel_index.parquet"))


if(!file.exists(pheno_file)){
out_files<-list.files(my_tempdir,".RData",full.names=T)


plan(multisession, workers = 10L)                # 10 cores on a 128-GB Mac
options(future.globals.maxSize = +Inf)           # lift the 500-MB ceiling

handlers(global  = TRUE)                         # enable once per session
handlers("progress")                             # nice txt bar in console
pboptions(type = "none")                         # turn off pbapply chatter

# helper that loads ONE .RData file (object ‘res’ inside)
load_one <- function(f){
  e <- new.env(parent = emptyenv())
  load(f, envir = e)                             # loads ‘res’
  e$res                                          # return it
}

# Split paths into ~20k-file chunks
chunk_size  <- 20000L
file_chunks <- split(out_files,
                     ceiling(seq_along(out_files) / chunk_size))

# Parallel read with progress bar
green_dat <- with_progress({
  p <- progressor(steps = length(file_chunks))   # one tick per chunk

  rbindlist(
    future_lapply(seq_along(file_chunks), function(idx){
      grp <- file_chunks[[idx]]

      # read every file in this chunk *serially* inside the worker
      ans <- rbindlist(lapply(grp, load_one),
                       use.names = TRUE, fill = TRUE)

      p(sprintf("chunk %d / %d", idx, length(file_chunks)))  # update bar
      ans                                                    # returned to master
    }, future.seed = FALSE),      # we don't need RNGs inside
    use.names = TRUE, fill = TRUE)
})

write_parquet(green_dat,pheno_file)
}else{
  green_dat<-read_parquet(pheno_file)
}

# create a *provisional* season label
add_season_id <- function(dt,
                          sos_col   = "TRS2.sos",
                          eps       = 70L,      # days
                          minPts    = 3L) {

  ## helper that returns a named integer vector (season ID)
  class_one_pixel <- function(sos_dates) {
    doy <- yday(sos_dates)
    # Drop NAs – they cannot be clustered
    ok  <- !is.na(doy)
    if (sum(ok) < minPts) return(rep(NA_integer_, length(doy)))

    cl <- dbscan(matrix(doy[ok]), eps = eps, minPts = minPts)$cluster

    if (all(cl == 0)) {              # DBSCAN found only noise
      return(rep(NA_integer_, length(doy)))
    }

    # Map DBSCAN cluster ID (1,2, …) → season number (1 = earliest)
    centers <- tapply(doy[ok], cl, median, na.rm = TRUE)
    ordering <- order(centers)              # earliest first
    season_map <- setNames(seq_along(ordering), names(centers)[ordering])

    out <- rep(NA_integer_, length(doy))
    out[ok] <- season_map[ as.character(cl) ]
    out
  }
  ## apply per pixel
  dt[ , `:=`(
    season = class_one_pixel(get(sos_col)),
    year   = year(get(sos_col))    # still handy to keep
  ), by = pixel]

  setattr(dt, "sorted", "pixel")
  dt[]
}

sos_col = "TRS2.sos"
eos_col = "TRS2.eos"

green_dat[, qscore := R2 - RMSE]   # ranges roughly −1 … +1
# Interpretation guide:
#   qscore > 0.8  → excellent fit
#   qscore 0.6-0.8→ usable but inspect if critical
#   qscore < 0.6  → drop for robust analyses

# Parameters for seasonal filtering and detections
q_score_min<-0.6
eps<-c(50,70,100)
minPts<-c(3,6,9,12)

summary_params<-expand.grid(q_score_min=q_score_min,eps=eps,minPts=minPts)

for(i in 1:nrow(summary_params)){
  cat("Summarizing sos using parameter set",i,"/",nrow(summary_params),"          \n")
  summary_name<-paste0("glass_sos_q",summary_params$q_score_min[i],
                       "_eps",summary_params$eps[i],
                       "_minpts",summary_params$minPts[i])


  green_dat_seasons<-add_season_id (dt=green_dat,
                        sos_col = "TRS2.sos",
                        eps     = summary_params$eps[i],   # ~70 days between peaks
                        minPts  = summary_params$minPts[i])


  # Summarize & Map

  green_dat_summary<-green_dat_seasons[!is.na(get(sos_col))  &
                                         qscore >= summary_params$q_score_min[i] ,
                                       .(sos_doy_med = as.integer(median(lubridate::yday(get(sos_col))), na.rm = TRUE),
                                         sos_doy_wmean = weighted.mean(yday(get(sos_col)), w = R2, na.rm = TRUE),
                                         sos_doy_mean = mean(yday(get(sos_col)), na.rm = TRUE),
                                         sos_doy_sd = mean(yday(get(sos_col)), na.rm = TRUE),
                                         eos_doy_med = as.integer(median(lubridate::yday(get(eos_col))), na.rm = TRUE),
                                         eos_doy_wmean = weighted.mean(yday(get(eos_col)), w = R2, na.rm = TRUE),
                                         eos_doy_mean = mean(yday(get(eos_col)), na.rm = TRUE),
                                         len_mean = mean(as.numeric(get(eos_col) - get(sos_col)), na.rm = TRUE),
                                         len_sd = sd(as.numeric(get(eos_col) - get(sos_col)), na.rm = TRUE)),
                                       by = .(pixel, season)][,sos_doy_cv:=sos_doy_sd/sos_doy_mean]

  # add lat/lon
  green_dat_summary<-merge(green_dat_summary,coords_index,all.x=T)

  # create rasters
  tmpl <- terra::rast(list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)[1])

  metrics <- setdiff(names(green_dat_summary),c("pixel","season","x","y"))

  make_layer <- function(dt_season, col_name, template)
  {
    pts <- terra::vect(
      dt_season[, .(x, y, value = get(col_name))],
      geom = c("x", "y"), crs = terra::crs(template)
    )
    lyr <- terra::rasterize(pts, template, field = "value", fun = "first")
    names(lyr) <- col_name
    lyr
  }

  season_ids <- sort(na.omit(unique(green_dat_summary$season)))

  # Limit to 3 seasons
  season_ids<-season_ids[season_ids<=3]

  stacked_seasons<-rast(lapply(season_ids,function(s){
    dt_s   <- green_dat_summary[season == s]

    ## make one layer per metric, then merge them into a multiband raster
    layers <- lapply(metrics, make_layer, dt_season = dt_s, template = tmpl)
    layers<-terra::rast(layers)
    names(layers)<-paste0("s",s,"_",names(layers))
    layers
  }))

  write_parquet(green_dat_summary,file.path(dirs$nvdi_phenology,paste0(summary_name,".parquet")))
  writeRaster(stacked_seasons,file.path(dirs$nvdi_phenology,paste0(summary_name,".tif")),overwrite=T)

}


# 8) Save Metadata & Future Enhancements ####
# - Store greenup variability (SD)
# - Store per-season cluster
# - Cross-validate against rainfall onset
# - Compare duration between green-up and green-down as proxy for season length

# End of Script ####
