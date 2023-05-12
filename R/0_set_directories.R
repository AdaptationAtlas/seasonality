if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}

packages<-c()

pacman::p_load(packages)


if(F){
# CGlabs
data_dir<-"/home/jovyan/common_data"
output_dir<-"/home/jovyan/common_data/atlas_sos"

}else{

# PS local machine
data_dir<-"D:/datasets"
output_dir<-"C:/datasets/seasonality"
}

chirps_dir<-paste0(data_dir,"/chirps/raw")
hobbins_dir<-paste0(data_dir,"/hobbins_ref_et/raw")
chirts_dir<-paste0(data_dir,"/chirts/raw")


if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = T)
}
