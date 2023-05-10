


DownloadCHIRPS<-function(StartYear=1980,EndYear=2022,EndDay=365,SaveDir=chirps_dir,quiet){

  URLmaster<-"https://data.chc.ucsb.edu/products/CHIRPS-2.0/africa_daily/tifs/p05/" # define the URL


  if(!dir.exists(SaveDir)){ # if the directory does not exist
    dir.create(SaveDir,recursive=T) # create the directory
  }

  for(YEAR in StartYear:EndYear){ # loop through the years

    if(YEAR==EndYear){ # if the year is the end year
      ENDDAY<-EndDay # set the end day to the end day
    }else{
      ENDDAY<-as.numeric(format(as.Date(paste0(YEAR,"-12-31")),"%j")) # set the end day to the last day of the year
    }

    # deal with leap years
    if(ENDDAY==365 & YEAR %in% seq(1978,2100,4)){ # if the end day is 365 and the year is a leap year
      ENDDAY<-366 # set the end day to 366
    }

    for(DAY in 1:ENDDAY){ # loop through the days


      DATE<-as.Date(paste0(YEAR,"-",DAY),format="%Y-%j") # set the date

      DAY<-format(DATE,"%d") # set the day
      MONTH<-format(DATE,"%m") # set the month

      FILE<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif.gz") # set the file name
      FILE2<-paste0("chirps-v2.0.",YEAR,".",MONTH,".",DAY,".tif") # set the file name

      URL<-paste0(URLmaster,YEAR,"/",FILE) # set the URL
      URL2<-paste0(URLmaster,YEAR,"/",FILE2) # set the URL
      destfile<-paste0(SaveDir,"/",FILE) # set the destination file
      destfile2<-paste0(SaveDir,"/",FILE2) # set the destination file





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

