if(!require(pacman)){
  install.packages("pacman")
  library(pacman)
}

packages<-c("devtools")

pacman::p_load(packages)


if(T){
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

# Starting conditions
# Is chirps dataset downloaded and complete?
chirps_dl<-T
# Is hobbins_ref_et dataset downloaded and complete?
hobbins_dl<-T
# Is chirts dataset downloaded and complete?
chirts_dl<-F

# If you need to set up github token
usethis::create_github_token()
gitcreds::gitcreds_set()
