# load the module
source("R/00_download_chc_daily.R")

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
                         end_year=1983,
                         out_dir=dirs$hobbins_ref_et,
                         bbox=af_bbox,
                         round_digits=1,
                         timeout=60,
                         workers = parallel::detectCores()%/%2)




