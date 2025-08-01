# load the module
source("R/00_download_chc_daily.R")
source("R/00_get_agera5_v2.R")

# Load packages
if (!requireNamespace("pacman", quietly=TRUE)) install.packages("pacman")
pacman::p_load(sf,elevatr,terra,data.table,terra,geodata,remotes,s3)

if(!require(MODIStsp)){
  remotes::install_github("lbusett/MODIStsp")
  require(MODIStsp)
}

# 1) Download chirps v3 ####
get_chirps_v3(
  start_year = 1981,
  end_year   = 2024,
  out_dir    = dirs$chirps_v3,
  bbox = af_bbox,
  round_digits=1,
  workers    = parallel::detectCores()-2)

# 2) Download chirts era5 - tmax ####

# Note funcion dev needs to output missing files/ingrity check option,
get_chirts_era5_tmax(start_year=1981,
                     end_year=2024,
                     out_dir=dirs$chirts_era5_max,
                     bbox=af_bbox,
                     round_digits=1,
                     workers= parallel::detectCores()-2,
                     timeout=60)

# 2) Download chirts era5 - tmin ####
get_chirts_era5_tmin(start_year=1981,
                     end_year=2024,
                     out_dir=dirs$chirts_era5_min,
                     bbox=af_bbox,
                     round_digits=1,
                     workers= parallel::detectCores()-2,
                     timeout=60)

# 2) Download chirts era5 - vpd - NOT YET AVAILBLE####
if(F){
get_chirts_era5_vpd(start_year=1981,
                     end_year=2024,
                     out_dir=dirs$chirts_era5_min,
                     bbox=af_bbox,
                     round_digits=1,
                     workers= parallel::detectCores()-2,
                     timeout=60)
}

# 4) Download PET ####
get_hobbins_refet(start_year=1983,
                         end_year=2024,
                         out_dir=dirs$hobbins_ref_et,
                         bbox=af_bbox,
                         round_digits=1,
                         timeout=60,
                         workers = parallel::detectCores()%/%2)

# 5) Download AgERA5 ####

user<-Sys.getenv("CDSAPI_UID", unset = NA)

ecmwfr::wf_set_key(user=user,
                   key=Sys.getenv("CDSAPI_KEY", unset = NA))

parameter_sets <- data.frame(
  variable = c("solar_radiation_flux", "2m_relative_humidity", "10m_wind_speed"),
  statistics = c(NA, NA, "24_hour_mean"),
  time = c(NA, "12_00", NA),
  stringsAsFactors = FALSE
)

years  <- as.character(1995:2024)
months <- sprintf("%02d", 1:12)

for (i in seq_len(nrow(parameter_sets))) {
  cat("Running parameter set",i,"/",nrow(parameter_sets),"\n")
  getAgERA5_v2(variable=parameter_sets$variable[i],
               stat=parameter_sets$statistics[i],
               time=parameter_sets$time[i],
               years=years,
               months=months,
               user = user,
               datadir=dirs$agera5_v2,
               bbox=af_bbox,
               workers=5,
               parallel_retries = 3)
}

# 6) Download Elevation ####

# Define study‐area bbox
corners <- matrix(
  c(af_bbox["xmin"], af_bbox["ymin"],
    af_bbox["xmin"], af_bbox["ymax"],
    af_bbox["xmax"], af_bbox["ymax"],
    af_bbox["xmax"], af_bbox["ymin"],
    af_bbox["xmin"], af_bbox["ymin"]),
  ncol = 2, byrow = TRUE
)
poly_sf <- st_polygon(list(corners)) |> st_sfc(crs = 4326)
poly_sf_feature <- st_sf( geometry = poly_sf )

# In most agro-ecological or growing season analyses, what matters is broad elevation gradients (low vs. highland)
# and large‐scale slope/aspect. A 1 km DEM captures those gradients just fine, and won’t artificially introduce noise
# that your climate layers can’t “see.”

dem_30m <- get_elev_raster(
  locations = poly_sf_feature,
  z       = 7,          # ~1.2km
  clip    = "bbox",      # crop exactly to the bbox you gave
  prj     = "EPSG:4326", # output in geographic
  src     = "aws"        # pull “SRTM 1 arc-second” from AWS terrain‐tiles
)

# inspect / save
plot(dem_30m)
terra::writeRaster(dem_30m,
                   filename = file.path(dirs$srtm,"DEM_SRTM_30m_Africa.tif"),
                   overwrite=TRUE)

# 7) Download soilgrids ####
depths<-c(20,50)
vars<-c("OC","sand","clay")
vd<-expand.grid(var=vars,depth=depths,stringsAsFactors = F)

lapply(1:nrow(vd),FUN=function(i){
  cat(paste(vd[i,],collapse="/"),"\n")
  soil_data <- geodata::soil_af_isda(
    var   = vd$var[i],       # organic carbon stock
    depth = vd$depth[i],     # topsoil
    path  = dirs$isda
  )
})

pblapply(tolower(vars),FUN=function(var){
  files<-list.files(dirs$isda,paste0(var,".*tif"),full.names = T,recursive=T)
  depth<-unlist(tstrsplit(basename(files),"_",keep=3))
  depths<-strsplit(gsub("cm","",depth),"-")
  new_depth<-paste0(min(as.numeric(unlist(depths))),"-",max(as.numeric(unlist(depths))),"cm")
  new_file<-file.path(dirs$isda,gsub(depth[1],new_depth,basename(files[1])))
  weights<-unlist(lapply(depths,FUN=function(d){mean(as.numeric(d))}))
  weights<-weights/sum(weights)
  rast_dat<-terra::rast(files)
  rast_dat<-terra::crop(rast_dat,af_v)
  rast_dat<-rast_dat * weights
  r_weighted <- sum(rast_dat)
  terra::writeRaster(r_weighted,new_file)
})

# 8) Download GLASS NDVI ####
  ## 8.1) Download hdf filee ####
  # PDF https://www.glass.hku.hk/pdf/algorithms/NDVI/Improved%20global%20250%20m%208-day%20NDVI%20and%20EVI%20products%20from%202000%E2%80%932021%20using%20the%20LSTM%20model.pdf

  for(year in 2000:2024){
    cat("Downloading glass nvdi year =",year,"\n")
    download_glass_ndvi(year = year,
                        out_dir = dirs$glass_ndvi,
                        overwrite = FALSE,
                        workers = future::availableCores() - 2)
  }

  ## 8.2) Convert to tif ####
  # We need to convert the data to a more friendly format, this is more easily done in the terminal than R
  # See ~/Documents/rprojects/seasonality/scripts/convert_glass_ndvi.sh


  ## 8.3) Delete bad files ####
  # Check the file conversion log to find any errors in the downloaded data
  log_path<-"scripts/gdal_errors.log"
  log_lines <- readLines(log_path)
  error_lines <- grep("❌", log_lines, value = TRUE)
  failed_files <- sub(".*ERROR processing ", "", error_lines)
  files<-list.files(dirs$glass_ndvi,paste(failed_files,collapse="|"),full.name=T)
  # If bad files, delete them files and then rerun 8.1
  unlink(files)
  if(length(files)>0){
    stop("Bad nvdi files present rerun section 8.1 and terminal tif conversion script")
  }


# 9) Download boundaries ####

  #   Create an S3FileSystem object for anonymous read access ####
  s3 <- s3fs::S3FileSystem$new(anonymous = TRUE)

  s3_files<-c(
    admin_0_file="s3://digital-atlas/domain=boundaries/type=admin/source=gaul2024/region=africa/processing=simplified/level=adm0/atlas_gaul24_a0_africa_simple-highres.parquet",
    admin_1_file="s3://digital-atlas/domain=boundaries/type=admin/source=gaul2024/region=africa/processing=simplified/level=adm1/atlas_gaul24_a1_africa_simple-highres.parquet"
  )

  local_files<-file.path(dirs$boundaries,basename(s3_files))

  lapply(1:length(local_files),FUN=function(i){
    file<-local_files[i]
    if(!file.exists(file)){
      s3$file_download(s3_files[i],file)
    }
  })
