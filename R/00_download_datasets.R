# load the module
source("R/00_download_chc_daily.R")
source("R/00_get_agera5_v2.R")
source("R/00_setup_folders.R")

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

dem <- get_elev_raster(
  locations = poly_sf_feature,
  z       = 7,          # ~1.2km
  clip    = "bbox",      # crop exactly to the bbox you gave
  prj     = "EPSG:4326", # output in geographic
  src     = "aws"        # pull “SRTM 1 arc-second” from AWS terrain‐tiles
)

# inspect / save
plot(dem)
terra::writeRaster(dem,
                   filename = file.path(dirs$srtm,"DEM_SRTM_Africa.tif"),
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
  ## 8.1) Download hdf file ####
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

# 10) Download LULC ####
  # You will need to download manually from:
    # https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=download
    # year = 2022
    # version = v2.1.1

# 11) Download Aridity Index (AI) ####
# Source: Global Aridity Index & ET0 Database v3 (Zomer et al. 2022) on Figshare
# (annual AI raster, 30 arc-sec). See dataset/readme for variable naming.   [oai_citation:3‡figshare.com](https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448?utm_source=chatgpt.com)

download_figshare_zip_by_doi <- function(doi, zip_name, dest_zip, overwrite = FALSE) {
  article_id <- sub("^.*figshare\\.(\\d+)\\..*$", "\\1", doi)
  if (!grepl("^[0-9]+$", article_id)) stop("Could not parse Figshare article id from DOI: ", doi)

  api  <- sprintf("https://api.figshare.com/v2/articles/%s", article_id)
  meta <- jsonlite::fromJSON(api)
  files <- meta$files

  hit <- files[files$name == zip_name, ]
  if (!nrow(hit)) {
    stop("ZIP not found: ", zip_name, "\nAvailable:\n- ", paste(files$name, collapse = "\n- "))
  }

  if (file.exists(dest_zip) && !overwrite) return(normalizePath(dest_zip))
  dir.create(dirname(dest_zip), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest_zip) && overwrite) file.remove(dest_zip)

  # ---- robust download for big files (curl resume + retries) ----
  url <- hit$download_url[1]

  cmd <- sprintf(
    "curl -L --fail --retry 10 --retry-delay 5 --connect-timeout 30 --max-time 0 -C - -o %s %s",
    shQuote(dest_zip), shQuote(url)
  )
  status <- system(cmd)
  if (status != 0) stop("curl download failed (status ", status, ").")

  normalizePath(dest_zip)
}

extract_and_find_tif <- function(zip_path, out_dir, tif_pattern = "\\.tif$") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(zip_path, exdir = out_dir)

  tifs <- list.files(out_dir, pattern = tif_pattern, recursive = TRUE, full.names = TRUE)
  if (!length(tifs)) stop("No GeoTIFFs found after unzip in: ", out_dir)

  tifs
}

ai_doi <- "10.6084/m9.figshare.7504448.v4"

zip_name <- "Global-AI_ET0__annual_v3_1.zip"
zip_path <- file.path(dirs$aridity, zip_name)

zip_path <- download_figshare_zip_by_doi(
  doi = ai_doi,
  zip_name = zip_name,
  dest_zip = zip_path,
  overwrite = TRUE
)

unzip_dir <- file.path(dirs$aridity, "Global-AI_ET0__annual_v3_1")
tifs <- extract_and_find_tif(zip_path, unzip_dir)
unlink(zip_path)  # clean up zip

ai_tif <- tifs[grepl("AI", basename(tifs), ignore.case = TRUE) &
              !grepl("ET0|PET", basename(tifs), ignore.case = TRUE)][1]

if (is.na(ai_tif) || !nzchar(ai_tif)) {
  stop("Could not uniquely identify AI GeoTIFF. Candidate files:\n- ", paste(basename(tifs), collapse="\n- "))
}

ai_r <- terra::rast(ai_tif)
ai_r_crop<-terra::crop(ai_r,af_v)
terra::writeRaster(ai_r_crop,
            filename = file.path(dirs$aridity,"Global_Aridity_Index_Africa.tif"),
            overwrite=TRUE)

unlink(dirname(ai_tif), recursive = TRUE)  # clean up unzipped files
plot(ai_r_crop)
rm(ai_r, ai_r_crop)


