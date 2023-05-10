#' Download CHIRPS precipitation data from UCSB website
#'
#' Download CHIRPS precipitation data from University of California Santa Barbara (UCSB) Climate Hazard Group website. The data are daily values for Africa, in GeoTIFF format.
#'
#' @param StartYear The year to start downloading data from. Default is 1980.
#' @param EndYear The year to stop downloading data at. Default is 2022.
#' @param EndDay The day to stop downloading data at for EndYear. Default is 365.
#' @param SaveDir The directory to save the downloaded data to. Default is chirps_dir.
#' @param quiet A logical parameter to suppress messages or not during download. Default is FALSE.
#' @return This function downloads CHIRPS dataset for the specified time period in the directory specified by SaveDir variable.
#' @export
DownloadCHIRPS<-function(StartYear=1980,EndYear=2022,EndDay=365,SaveDir=chirps_dir,quiet){

  URLmaster<-"https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_daily/tifs/p05/" # define the URL

  # Create SaveDir directory if it does not exist
  if(!dir.exists(SaveDir)){
    dir.create(SaveDir,recursive=T)
  }

  # loop through the years
  for(YEAR in StartYear:EndYear){

    # Set the correct number of days in the year
    if(YEAR==EndYear){ # if the year is the end year
      ENDDAY<-EndDay # set the end day to the end day
    }else{
      ENDDAY<-as.numeric(format(as.Date(paste0(YEAR,"-12-31")),"%j")) # set the end day to the last day of the year
    }

    # Set the correct number of days in leap years
    if(ENDDAY==365 & YEAR %in% seq(1978,2100,4)){ # if the end day is 365 and the year is a leap year
      ENDDAY<-366 # set the end day to 366
    }

    # loop through the days of each year
    for(DAY in 1:ENDDAY){

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
      destfile<-paste0(SaveDir,"/",FILE) # set the destination file
      destfile2<-paste0(SaveDir,"/",FILE2) # set the destination file

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

options(timeout = 150000)
DownloadCHIRPS(StartYear=2022,EndYear = 2022,EndDay=335,SaveDir = SaveDir,quiet=T)

