# Download chirps
if(!chirps_dl){
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_chirps.R")

  options(timeout = 150000)
  download_chirps(start_year=1981,end_year = 2022,end_day=365,save_dir = chirps_dir,quiet=T)

}

# Download hobbins
if(!hobbinss_dl){
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_hobbins.R")

  options(timeout = 150000)
  download_hobbins(start_year=1981,end_year = 2022,end_day=365,save_dir = chirps_dir,quiet=T)

}
