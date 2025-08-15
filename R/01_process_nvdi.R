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

#' - First run:
#' 00_setup_folders.R
#' 00_download_glass_nvdi.R
#' 00_download_datasets.R
#' 01_process_chirps.R

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
  miceadds,
  geoarrow,
  sf,
  miceadds,
  s3fs
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
}


sos_cols = c("TRS2.sos","Greenup")
eos_cols = c("TRS2.eos","Senescence")


# Interpretation guide:
#   qscore > 0.8  → excellent fit
#   qscore 0.6-0.8→ usable but inspect if critical
#   qscore < 0.6  → drop for robust analyses

q_score_min<-0.6

summary_params<-data.table(q_score_min=q_score_min,sos_col=sos_cols,eos_cols=eos_cols)

# Create helper functions
day2rad <- function(d) 2 * pi * d / 365                # DOY  →  radians

circ_mean <- function(d, w = NULL) {
  ang <- day2rad(d)
  if (is.null(w)) w <- rep(1, length(ang))
  S <- sum(w * sin(ang), na.rm = TRUE)
  C <- sum(w * cos(ang), na.rm = TRUE)
  (atan2(S, C) %% (2 * pi)) * 365 / (2 * pi)           # back to DOY
}

circ_sd <- function(d, w = NULL) {
  ang <- day2rad(d)
  if (is.null(w)) w <- rep(1, length(ang))
  S <- sum(w * sin(ang), na.rm = TRUE)
  C <- sum(w * cos(ang), na.rm = TRUE)
  R <- sqrt(S^2 + C^2) / sum(w[!is.na(ang)])
  sqrt(-2 * log(R)) * 365 / (2 * pi)                   # DOY units
}

season_len <- function(sos, eos) {
  # Compute DOY with NA handling
  doy_sos <- suppressWarnings(lubridate::yday(sos))
  doy_eos <- suppressWarnings(lubridate::yday(eos))

  # Initialize length as NA
  len <- rep(NA_integer_, length(sos))

  # Identify valid indices
  valid_idx <- which(!is.na(doy_sos) & !is.na(doy_eos))

  # Compute season length for valid entries
  len[valid_idx] <- doy_eos[valid_idx] - doy_sos[valid_idx]
  len[valid_idx][len[valid_idx] < 0] <- len[valid_idx][len[valid_idx] < 0] + 365  # Wrap-around

  return(as.integer(len))
}

field_descriptions <- list(
  pixel             = "Unique numeric identifier for the pixel location.",
  season            = "Season number within the year (1 = first detected season, 2 = second season if bimodal).",
  n_seasons         = "Total number of valid seasons detected for the pixel across all years.",
  sos_doy_med       = "Median Start of Season (SoS) expressed as day of year (DOY) across all valid years.",
  sos_doy_mean      = "Mean Start of Season (SoS) in DOY.",
  sos_doy_wmean     = "Weighted mean of SoS in DOY, weighted by seasonal R².",
  sos_doy_sd        = "Standard deviation of SoS in DOY.",
  eos_doy_med       = "Median End of Season (EoS) expressed as day of year (DOY).",
  eos_doy_mean      = "Mean End of Season (EoS) in DOY.",
  eos_doy_wmean     = "Weighted mean of EoS in DOY, weighted by seasonal R²",
  eos_doy_sd        = "Standard deviation of EoS in DOY.",
  len_mean          = "Mean growing season length in days (EoS − SoS), adjusted for year wraparound.",
  len_sd            = "Standard deviation of growing season length.",
  useable_seasons   = "Number of seasons that passed filtering thresholds and were used in statistics.",
  sos_doy_cv        = "Coefficient of variation (CV) of SoS in DOY, calculated as sd/mean.",
  x                 = "Longitude of the pixel centroid.",
  y                 = "Latitude of the pixel centroid.",
  failed_seasons    = "Proportion of candidate seasons rejected due to low confidence or failed fitting (range 0 to 1)."
)

fwrite(field_descriptions,file.path(dirs$nvdi_phenology,"glass_sos_metadata.csv"))

for(i in 1:nrow(summary_params)){

  cat("Summarizing sos using parameter set",i,"/",nrow(summary_params),"          \n")

  sos_col<-summary_params$sos_col[i]
  eos_col<-summary_params$eos_col[i]

  summary_name<-paste0("glass_sos_q",summary_params$q_score_min[i],
                       "_sos-",sos_col,
                       "_eos-",eos_col)


  green_dat<-read_parquet(pheno_file)

  # Parameters for seasonal filtering and detections
  green_dat[, qscore := R2 - RMSE]   # ranges roughly −1 … +1


  # create a *provisional* season label
  green_dat[,season:=tstrsplit(flag,"_",keep=2)]
  green_dat<-green_dat[!is.na(get(sos_col))]

  green_dat<-green_dat[!is.na(get(sos_col)),n_seasons:=.N,by=.(pixel,season)]
  (max_seasons<-green_dat[,max(n_seasons,na.rm=T)])

  # Summarize & Map
  green_dat_summary <- green_dat[
    qscore >= summary_params$q_score_min[i],
    .(
      sos_doy_med   = as.integer(median(lubridate::yday(get(sos_col)), na.rm = TRUE)),
      sos_doy_mean  = circ_mean(  lubridate::yday(get(sos_col))          ),
      sos_doy_wmean = circ_mean(  lubridate::yday(get(sos_col)), w = R2  ),
      sos_doy_sd    = circ_sd(    lubridate::yday(get(sos_col))          ),

      eos_doy_med   = as.integer(median(lubridate::yday(get(eos_col)), na.rm = TRUE)),
      eos_doy_mean  = circ_mean(  lubridate::yday(get(eos_col))          ),
      eos_doy_wmean = circ_mean(  lubridate::yday(get(eos_col)), w = R2  ),
      eos_doy_sd    = circ_sd(    lubridate::yday(get(eos_col))          ),

      len_mean      = mean(season_len(get(sos_col), get(eos_col)), na.rm = TRUE),
      len_sd        =  sd(season_len(get(sos_col), get(eos_col)), na.rm = TRUE),

      useable_seasons = .N
    ),
    by = .(pixel, season,n_seasons)
  ][, sos_doy_cv := sos_doy_sd / sos_doy_mean]

  # add lat/lon
  green_dat_summary<-merge(green_dat_summary,coords_index,all.x=T)

  green_dat_summary[,failed_seasons:=round(1-useable_seasons/n_seasons,3)]

  # create rasters
  tmpl <- terra::rast(list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)[1])

  metrics <- setdiff(names(green_dat_summary),c("pixel","season","x","y","n_seasons","useable_seasons"))

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

  rm(green_dat_summary,green_dat,green_dat_seasons)
  gc()
}

# Explore results
if(F){
r <- rast(ncol=100, nrow=100)
values(r) <- sample(1:365, ncell(r), replace=TRUE)
# Create a circular palette (e.g. HSV hue from 0 to 360 degrees)
require(viridisLite)
circular_colors <- rainbow(365, start = 0, end = 1)  # Full color wheel
plot(stacked_seasons[[grep("sos_doy_med",names(stacked_seasons),value=T)]],col = circular_colors)
plot(stacked_seasons[[grep("sos_doy_wmean",names(stacked_seasons),value=T)]],col = circular_colors)
plot(stacked_seasons[[grep("len_mean",names(stacked_seasons),value=T)]])
plot(stacked_seasons[[grep("failed_seasons",names(stacked_seasons),value=T)]])
}

# 7) Save Metadata & Future Enhancements ####
# - Store greenup variability (SD)
# - Store per-season cluster
# - Cross-validate against rainfall onset
# - Compare duration between green-up and green-down as proxy for season length

# 8) Create country subsets ####

save_dir<-file.path(dirs$nvdi_phenology,"countries")
if(!dir.exists(save_dir)){dir.create(save_dir)}

file<-list.files(dirs$nvdi_phenology,"pixel_index",full.names = T)
pixel_index<-read_parquet(file)

pixel_map <- rast(
  x = pixel_index,
  type = "xyz",     # indicates x, y, value columns
  crs = "EPSG:4326" # adjust if needed
)
#  Assign values using the pixel column
values(pixel_map) <- pixel_index$pixel

# Admin 0 map
africa0<-list.files(dirs$boundaries,"_a0_",full.names = T)
africa0<-read_parquet(africa0)

africa0_sf <- st_as_sf(africa0)
africa0_vect <- vect(africa0_sf)


# Admin 1 map
africa1<-list.files(dirs$boundaries,"_a1_",full.names = T)
africa1<-read_parquet(africa1)

africa1_sf <- st_as_sf(africa1)
africa1_vect <- vect(africa1_sf)

countries<-sort(unique(africa0_vect$iso3))

for(country in countries){
  cat("Processing",country,"         /r")
  save_file<-file.path(save_dir,paste0(country,"_seasonal-phenology.parquet"))

  country0<-africa0_vect[africa0_vect$iso3==country,]
  country1<-africa1_vect[africa1_vect$iso3==country,c("admin1_name")]
  country1 <- aggregate(country1, by = "admin1_name")

  # Subset pixel map to country
  pixel_map_cr<-mask(crop(pixel_map,country0),country0)

  # Extract pixels by county
  px_by_county <- data.table(extract(pixel_map_cr, country1))
  country1_dat<-as.data.table(country1)[,ID:=.I][,agg_n:=NULL]
  px_by_county <-merge(px_by_county,country1_dat,all.x=T)

  # Create pheno files for country
  pheno_files<-data.table(file=file.path(dirs$nvdi_phenology,"temp",paste0(px_by_county$pixel,".RData")))
  pheno_files[,exists:=file.exists(file)]

  missing<-pheno_files[exists==F,.N]
  perc<-round(100*missing/nrow(pheno_files),2)

  if(missing>0){
    cat("Warning: ",missing," missing pixels.",perc,"%")
  }

  pheno_files<-pheno_files[exists==T]

  # Load pheno data for country
  country_dat<-rbindlist(pblapply(1:nrow(pheno_files),function(i){
    x<-load.Rdata2(basename(pheno_files$file[i]),path=dirname(pheno_files$file[i]))
    return(x)
  }))

  # Merge admin ID
  country_dat<-merge(country_dat,px_by_county[,!"ID"],all.x=T)

  # Save result
  write_parquet(country_dat,save_file)
}

# 8.1) Merge rainfall data with country subsets ####
chirps_dir<-file.path(dirname(dirs$chirps_v3),paste0(basename(dirs$chirps_v3),"_countries"))

for(i in 1:length(countries)){
  country<-countries[i]
  cat("Processing",country,i,"/",length(countries),"         \r")
  nvdi_file<-file.path(save_dir,paste0(country,"_seasonal-phenology.parquet"))
  chirps_file<-file.path(chirps_dir,paste0(country,".parquet"))
  save_file<-file.path(save_dir,paste0(country,"_seasonal-phenology_plus-rain.parquet"))

  if(!file.exists(save_file)){
   nvdi_dat<-read_parquet(nvdi_file)
   chirps_dat<-read_parquet(chirps_file)

   setDT(chirps_dat)
   setDT(nvdi_dat)

   # Convert to compact date type
   date_cols <- names(which(sapply(nvdi_dat, function(x) class(x)[1] == "Date")))
   nvdi_dat[, (date_cols) := lapply(.SD, function(x) as.IDate(as.character(x))),
            .SDcols = date_cols]

   chirps_dat[, date2 := as.IDate(date[1], format = "%Y-%m-%d"),by=date]
   chirps_dat[, date3 := as.IDate(date2[1]),by=date2]
   chirps_dat[,c("date","date2"):=NULL]
   setnames(chirps_dat,"date3","date")

   # Helpful index for speed on the big table
   setkey(chirps_dat, pixel, date)

   # Assign the season flag to CHIRPS rows within the TRS2 interval (inclusive)
   chirps_dat[
     nvdi_dat,
     on = .(pixel, date >= TRS2.sos, date <= TRS2.eos),
     flag := i.flag
   ]

   chirps_trs2<-chirps_dat[!is.na(flag),.(rain_trs2=sum(value,na.rm=T)),by=.(x,y,pixel,flag)]

   chirps_dat[,flag:=NULL
              ][nvdi_dat,
                on = .(pixel, date >= TRS5.sos, date <= TRS5.eos),
                flag := i.flag
                ]

   chirps_trs5<-chirps_dat[!is.na(flag),.(rain_trs5=sum(value,na.rm=T)),by=.(x,y,pixel,flag)]

   chirps_dat[,flag:=NULL
   ][nvdi_dat,
     on = .(pixel, date >= DER.sos, date <= DER.eos),
     flag := i.flag
   ]

   chirps_der<-chirps_dat[!is.na(flag),.(rain_der=sum(value,na.rm=T)),by=.(x,y,pixel,flag)]


   chirps_dat[,flag:=NULL
   ][nvdi_dat[!is.na(Greenup) & !is.na(Senescence)],
     on = .(pixel, date >= Greenup, date <= Senescence),
     flag := i.flag
   ]

   chirps_gs<-chirps_dat[!is.na(flag),.(rain_gs=sum(value,na.rm=T)),by=.(x,y,pixel,flag)]

   # Merge back rainfall
   nvdi_dat[chirps_trs2,
     on = .(pixel, flag),
     rain_trs2:= i.rain_trs2]

   nvdi_dat[chirps_trs5,
            on = .(pixel, flag),
            rain_trs5 := i.rain_trs5]

   nvdi_dat[chirps_gs,
            on = .(pixel, flag),
            rain_greenup_scenescence := i.rain_gs]

   nvdi_dat[chirps_der,
            on = .(pixel, flag),
            rain_der := i.rain_der]

   write_parquet(nvdi_dat,save_file)

  }

}

  ## 8.2) Upload to S3 ####
# Includes functions to upload data to S3 bucket
source("https://raw.githubusercontent.com/AdaptationAtlas/hazards_prototype/refs/heads/main/R/haz_functions.R")

s3_bucket<-"s3://digital-atlas/domain=phenology/type=sos_eos/source=glass_nvdi/region=africa/processing=raw/level=adm1/"
folder_local<-file.path(dirs$nvdi_phenology,"countries")

files<-list.files(folder_local,full.names = T,recursive=F,include.dirs = F)

upload_files_to_s3(files = files,
                   selected_bucket=s3_bucket,
                   max_attempts = 3,
                   overwrite=F,
                   mode="public-read")


# End of Script ####
