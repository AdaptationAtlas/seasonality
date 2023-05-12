#' Download CHIRPS precipitation data from UCSB website
#'
#' Download CHIRPS precipitation data from University of California Santa Barbara (UCSB) Climate Hazard Group website. The data are daily values for Africa, in GeoTIFF format.
#'
#' @param start_year The year to start downloading data from. Default is 1980.
#' @param end_year The year to stop downloading data at. Default is 2022.
#' @param end_day The day to stop downloading data at for end_year. Default is 365.
#' @param save_dir The directory to save the downloaded data to. Default is chirps_dir.
#' @param quiet A logical parameter to suppress messages or not during download. Default is FALSE.
#' @return This function downloads CHIRPS dataset for the specified time period in the directory specified by save_dir variable.
#' @export
download_chirps<-function(start_year=1980,end_year=2022,end_day=365,save_dir=chirps_dir,quiet=F){

  URLmaster<-"https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_daily/tifs/p05/" # define the URL

  # Create save_dir directory if it does not exist
  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive=T)
  }

  # loop through the years
  for(YEAR in start_year:end_year){

    # Set the correct number of days in the year
    if(YEAR==end_year){ # if the year is the end year
      end_day<-end_day # set the end day to the end day
    }else{
      end_day<-as.numeric(format(as.Date(paste0(YEAR,"-12-31")),"%j")) # set the end day to the last day of the year
    }

    # Set the correct number of days in leap years
    if(end_day==365 & YEAR %in% seq(1978,2100,4)){ # if the end day is 365 and the year is a leap year
      end_day<-366 # set the end day to 366
    }

    # loop through the days of each year
    for(DAY in 1:end_day){

      # Set the date variables
      DATE<-as.Date(paste0(YEAR,"-",DAY),format="%Y-%j") # set the date
      DAY<-format(DATE,"%d") # set the day
      MONTH<-format(DATE,"%m") # set the month

      # Set the file variables
      FILE<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif.gz") # set the file name
      FILE2<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif") # set the file name

      # Set the URL and destination file variables
      URL<-paste0(URLmaster,YEAR,"/",FILE) # set the URL
      URL2<-paste0(URLmaster,YEAR,"/",FILE2) # set the URL
      destfile<-paste0(save_dir,"/",FILE) # set the destination file
      destfile2<-paste0(save_dir,"/",FILE2) # set the destination file

      # Download the file if it does not exist already
      if(!(file.exists(destfile)|file.exists(destfile2))){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",FILE2))
        flush.console()
        try(download.file(URL2, destfile2,quiet=quiet))
      }
    }
  }
}
