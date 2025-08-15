
# Load packages ####
pacman::p_load(
  arrow,
  geoarrow,
  dbscan,
  data.table,    # fast data manipulation
  terra,         # raster operations
  sf,
  pbapply
)

# Set directories ####
save_dir<-file.path(dirname(dirs$chirps_v3),paste0(basename(dirs$chirps_v3),"_countries"))
if(!dir.exists(save_dir)){dir.create(save_dir)}


# Load NVDI pixel map ###
file<-list.files(dirs$nvdi_phenology,"pixel_index",full.names = T)
pixel_index<-read_parquet(file)

pixel_map <- rast(
  x = pixel_index,
  type = "xyz",     # indicates x, y, value columns
  crs = "EPSG:4326" # adjust if needed
)
#  Assign values using the pixel column
values(pixel_map) <- pixel_index$pixel

# Load admin vector ####
# Admin 0 map
africa0<-list.files(dirs$boundaries,"_a0_",full.names = T)
africa0<-read_parquet(africa0)

africa0_sf <- st_as_sf(africa0)
africa0_vect <- vect(africa0_sf)[,"iso3"]
africa0_vect <- aggregate(africa0_vect, by = "iso3")
africa0_vect$area<-expanse(africa0_vect)
africa0_vect<-africa0_vect[order(africa0_vect$area),]
countries<-africa0_vect$iso3

# List rain files ####
rain_files<-list.files(dirs$chirps_v3,".tif$",recursive=T,full.names = T)

# Check for any files that have not been cropped and fix
if(F){
ext_check<-pblapply(rain_files,FUN=function(file){
  terra::ext(terra::rast(file))[1:4]
})

ext_check<-do.call("rbind",ext_check)
ext_check<-data.table(ext_check)
ext_check[,ext:=paste(xmin,xmax,ymin,ymax)
          ][,factor:=as.numeric(as.factor(ext))
            ][,file:=rain_files]

if(ext_check[,length(unique(factor))]>0){
bad_files<-ext_check[factor==1,file]
good_rast<-rast(ext_check[factor==2,file][1])

for(file in bad_files){
  cat(file,"                  \r")
  bad_rast<-rast(file)+0
  fixed_rast<-crop(bad_rast,good_rast)
  stopifnot(ext(fixed_rast)==ext(good_rast))
  writeRaster(fixed_rast,file,overwrite=T,
              filetype = "COG",
              gdal = c("COMPRESS=ZSTD", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512"))
}

}
}

rain_dat<-rast(rain_files)

# This loop is converts a raster stack of daily rainfall for a country in a long-format data.table where each row
# is a day x pixel rainfall value. Pixel values are added to each row so we can merge with the glass nvdi data.

for(i in 1:length(countries)){
  country<-countries[i]
  cat("Processing",country,i,"/",length(countries),"         \r")
  save_file<-file.path(save_dir,paste0(country,".parquet"))

  if(!file.exists(save_file)){

  # Get vector for the selected country
  country0<-africa0_vect[africa0_vect$iso3==country,]

  # Crop and mask pixel map to country vector
  pixel_map_cr<-mask(crop(pixel_map,country0),country0)

  # Crop rainfall to country vector
  rain_crop<-crop(rain_dat,country0)

  # Resample rainfall to pixel map
  rain_crop_rs<-resample(rain_crop,pixel_map_cr)
  rm(rain_crop)
  gc()

  # Mask rainfall to country vector
  rain_crop_rs<-mask(rain_crop_rs,country0)
  rain_crop_rs<-c(pixel_map_cr,rain_crop_rs)

  # Get raster values as a wide data.table
  rain_dt<-data.table(as.data.frame(rain_crop_rs,xy = TRUE, na.rm = FALSE))
  rm(rain_crop_rs)
  gc()

  # Melt into long form
  rain_dt_melt<-melt(rain_dt,id.vars=c("x","y","pixel"))
  rm(rain_dt)
  gc()

  # Remove na values to reduce size
  rain_dt_melt<-rain_dt_melt[!is.na(pixel)]
  # Convert variable (file) name to date format
  rain_dt_melt[,variable:=gsub("chirps-v3.0.","",variable[1]),by=variable]
  rain_dt_melt[,variable:=as.Date(variable[1],format = "%Y.%m.%d"),by=variable]
  setnames(rain_dt_melt,"variable","date")

  # Round values to reduce file size
  rain_dt_melt[,value:=round(value,2)]

  # Save result in parquet format
  write_parquet(rain_dt_melt,save_file)
  rm(rain_dt_melt)
  gc()

  }
}


