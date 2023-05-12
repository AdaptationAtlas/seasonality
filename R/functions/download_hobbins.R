#' Download Hobbins' Reference Evapotranspiration data
#'
#' Downloads Hobbins' Reference Evapotranspiration (PET) tif.gz files from https:/data.chc.ucsb.edu/products/Hobbins_RefET/ETos_p10_dekad_africa/tifs_gz/ and unzips them into the given destination directory.
#'
#' @param start_year starting year of data to download (defaults to 1980)
#' @param end_year ending year of data to download (defaults to 2022)
#' @param end_dekad ending dekad of data to download (defaults to 36)
#' @param save_dir destination directory to save the downloaded files
#' @param quiet whether or not to suppress the output of the download progress (defaults to TRUE)
#' @param unzip whether or not to unzip the downloaded files (defaults to TRUE)
#' @param remove whether or not to remove the original downloaded compressed files (defaults to TRUE)
#'
#' @return A message acknowledging the completion of the file download and unzip process.
#'
#' @import download_and_unzip
#'
#' @examples
#' download_hobbins(start_year = 1990, end_year = 1995, save_dir = "./data/HobbinsRefET/")
#'
#' @export
download_hobbins<-function(start_year=1980,end_year=2022,end_dekad=36,save_dir,quiet=T,unzip=T,remove=T){


  if(!is.na(save_dir) & substr(save_dir,nchar(save_dir),nchar(save_dir))!="/"){
    save_dir<-paste0(save_dir,"/")
  }

  URLmaster<-"https://data.chc.ucsb.edu/products/Hobbins_RefET/ETos_p10_dekad_africa/tifs_gz/"


  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive=T)
  }


  for(YEAR in start_year:end_year){ # MinYear:MaxYear (Min = 1983 Max = Present)

    if(YEAR==end_year){
      end_dekad1<-end_dekad
    }else{
      end_dekad1<-36
    }

    for(DEKAD in 1:end_dekad1){

      if(DEKAD<10){
        DEKAD<-paste0("0",as.character(DEKAD))
      }

      # Set the file variables
      FILE<-paste0("pet_",YEAR,DEKAD,".tif.gz")
      FILE2<-paste0("pet_",YEAR,DEKAD,".tif")

      # Set the URL and destination file variables
      URL<-paste0(URLmaster,"/",FILE) # set the URL
      URL2<-paste0(URLmaster,"/",FILE2) # set the URL
      destfile<-paste0(save_dir,"/",FILE) # set the destination file
      destfile2<-paste0(save_dir,"/",FILE2) # set the destination file


      # Download the file if it does not exist already
      if(!(file.exists(destfile)|file.exists(destfile2))){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",FILE))
        flush.console()
        try(download_and_unzip(URL, filename=destfile2,destname=destfile,quiet=quiet,unzip=unzip,remove=remove))
      }else{
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("File present: ",FILE))
        flush.console()
      }

    }
  }
  cat("Download and unzip complete.")
}
