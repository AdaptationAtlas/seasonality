#!/usr/bin/env Rscript
# -----------------------------------------------------------
# 00_setup_folders.R   (run once after you mounted the SSD)
# -----------------------------------------------------------
Sys.setenv(ANALOGUE_DATA_ROOT = "/Volumes/clim_dat")
root <- Sys.getenv("ANALOGUE_DATA_ROOT")

Sys.setenv(PATH = paste("/opt/homebrew/bin", Sys.getenv("PATH"), sep = ":"))
Sys.which("gdal_translate")


if (!dir.exists(root)) {
  stop("⚠️  ANAlOGUE_DATA_ROOT does not exist.")
}

dirs <- list(
  chirps_v3="climate_raw/chirps/chirps_v3_cog",
  chirts_era5_max="climate_raw/chirts_era5/tmax_cog",
  chirts_era5_min="climate_raw/chirts_era5/tmin_cog",
  agera5_v2="climate_raw/agera5_v2",
  hobbins_ref_et="climate_raw/hobbins_ref_et",
  boundaries="static_raw/boundaries",
  srtm="static_raw/strm",
  soilgrids="static_raw/soilgrids",
  isda="static_raw/isda",
  glass_ndvi="climate_raw/glass_ndvi",
  glass_ndvi_tif="climate_raw/glass_ndvi_tif",
  onset="climate_derived/onset",
  nvdi_phenology="climate_derived/glass_phenology",
  seasonal="climate_derived/seasonal",
  similarity="climate_derived/similarity",
  rf_forests="models/rf_forests",
  validation_reports="models/validation_reports",
  points="points",
  scratch="scratch"
)

dir_names<-names(dirs)
dirs<-lapply(dirs,FUN=function(d){
  p <- file.path(root, unlist(d))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
})
names(dirs)<-dir_names

# drop a meta file
meta <- list(
  created  = as.character(Sys.time()),
  hostname = Sys.info()[["nodename"]],
  root     = root
)
jsonlite::write_json(meta, file.path(root, "climate_derived", "setup_meta.json"), auto_unbox = TRUE, pretty = TRUE)

cat("✅ Folder structure created under", root, "\n")

list_dirs <- function(path, depth = 2, indent = 0) {
  cat(strrep(" ", indent), basename(path), "\n", sep = "")
  if (depth == 0) return()
  dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  for (d in dirs) list_dirs(d, depth - 1, indent + 2)
}

print(list_dirs(root, depth = 3))

# Make africa bbox
af_bbox<-c(
xmin  =  -26,   # 26 ° W  (just west of Cabo Verde)
xmax  =   64,   # 64 ° E  (just east of Rodrigues at 63.4 ° E)
ymin  =  -35,  # 35 ° S  (a bit south of Cape Agulhas at 34.8 ° S)
ymax  =   38   # 38 ° N  (north of Morocco’s northern tip at 35.9 ° N)
)

af_v <- terra::vect(
  terra::ext(af_bbox["xmin"], af_bbox["xmax"], af_bbox["ymin"], af_bbox["ymax"]),
  crs = "EPSG:4326"
)
