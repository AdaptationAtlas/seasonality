#' Download CHIRTS Daily Climate Data
#'
#' Downloads daily climate data from CHIRTS for Africa with a spatial resolution of p05
#'
#' @param start_year The year to start downloading data (1983 by default).
#' @param end_year The year to stop downloading data (2016 by default).
#' @param save_dir The directory where the data will be saved (chirps_dir by default).
#' @param quiet If TRUE, progress will not be printed to the console. FALSE by default.
#'
#' @importFrom dir.create Create the save_dir directory if it does not exist
#' @importFrom file.exists Check if the destination file exists
#' @importFrom download.file Download a file from a URL.
#'
#' @examples
#' \dontrun{
#' download_chirts(start_year = 1990, end_year = 1992)
#' }
#'
#' @export
download_chirts<-function(start_year=1983,end_year=2016,save_dir=chirps_dir,quiet=F){

  URLmaster<-"https://data.chc.ucsb.edu/products/CHIRTSdaily/v1.0/africa_netcdf_p05/" # define the URL

  # Create save_dir directory if it does not exist
  if(!dir.exists(save_dir)){
    dir.create(save_dir,recursive=T)
  }

  # loop through the years
  for(YEAR in start_year:end_year){

    # variables
    for(VAR in c("Tmax","Tmin")){

      # Set the file variables
      FILE<-paste0(VAR,".",YEAR,".nc") # set the file name

      # Set the URL and destination file variables
      URL<-paste0(URLmaster,"/",FILE) # set the URL
      destfile<-paste0(save_dir,"/",FILE) # set the destination file

      # Download the file if it does not exist already
      if(!(file.exists(destfile))){
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("Downloading file: ",FILE))
        flush.console()
        try(download.file(URL, destfile,quiet=quiet))
      }else{
        # Display progress
        cat('\r                                                                                                                                          ')
        cat('\r',paste0("File present: ",FILE))
        flush.console()
      }
    }
  }
}
