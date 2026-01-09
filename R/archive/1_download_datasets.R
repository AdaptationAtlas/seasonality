# Download chirps
if(!chirps_dl){
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_chirps.R")

  options(timeout = 150000)
  download_chirps(start_year=1980,end_year = 2022,end_day=365,save_dir = chirps_dir,quiet=T)

}

# Download hobbins
if(!hobbins_dl){
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_hobbins.R")
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_and_unzip.R")

  options(timeout = 150000)
  download_hobbins(start_year=1980,end_year = 2022,end_dekad=1,save_dir = hobbins_dir,quiet=T)

}

# Download chirts
if(!chirts_dl){
  devtools::source_url("https://raw.githubusercontent.com/AdaptationAtlas/seasonality/main/R/functions/download_chirts.R")
  download_chirts(start_year=1983,end_year=2016,save_dir=chirps_dir,quiet=F)
}
