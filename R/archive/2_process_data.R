devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/sos_functions.R")


chirps_save_dir<-paste0(chirps_dir,"/intermediate/array_countries_dekad/temp")
chirps_save_dir_dekad<-paste0(chirps_dir,"/intermediate/array_countries_dekad")

if(!dir.exists(chirps_save_dir)){
  dir.create(chirps_save_dir,recursive = T)
}

if(!dir.exists(chirps_save_dir_dekad)){
  dir.create(chirps_save_dir_dekad,recursive = T)
}

# List CHIRPS files
FILES<-data.frame(File=list.files(path = chirps_dir,'.tif'))
FILES<-droplevels(FILES)

Dates<-unlist(lapply(strsplit(FILES$File,"v2.0."),"[[",2))
Dates<-unlist(lapply(strsplit(Dates,".tif"),"[[",1))
Dates<-as.Date(Dates,format = "%Y.%m.%d")

FILES$Year<-format(Dates,"%Y")
FILES$Month<-format(Dates,"%m")
FILES$Day<-format(Dates,"%d")
FILES$Date<-Dates
FILES$Yday<-format(Dates,"%j")
FILES$Code<-paste(FILES$Year,FILES$Yday,sep="-")

# Load a CHIRPS raster to get geometries
CHIRPSrast<-terra::rast(paste0(chirps_dir,"/",FILES$File[1]))
