#' Get watershed stressor data
#'
#' @param boo_ws boolean directing CASTool to run the watershed stressor analysis
#' @param boo_helper boolean directing CASTool to download data from the CASTool helper package
#' @param out_dir directory where data outputs are written
#' @param dn_checked_sk name of checked data directory
#' @param region analysis region name
#' @param config list containing the configuration parameters
#'
#' @return A list with watershed stressor data and metadata.
#' If the user elects to not run the watershed stressor module, the function return list with data and metadata set to NULL.
#' @export
getWSStressData <- function(boo_ws = NULL,
                            boo_helper = NULL,
                            out_dir = NULL,
                            dn_checked_sk = NULL,
                            region = NULL,
                            config = NULL){

  if(is.null(config) == FALSE){
    out_dir <- config$out_dir
    dn_checked_sk <- config$dn_checked_sk
    region <- config$region
  }

  if(boo_ws == FALSE){
    data_stressorWS <- NULL
    data_stressorinfoWS <- NULL

  } else if (boo_ws == TRUE & boo_helper == FALSE) {
    data_stressorWS     <- readRDS(file.path(out_dir, dn_checked_sk, "data_stressorWS.rds"))
    data_stressorinfoWS <- readRDS(file.path(out_dir, dn_checked_sk,
                                             "data_stressorinfoWS.rds"))
    if(exists("data_stressorWS") & exists("data_stressorinfoWS")){
      message("Watershed stressor data provided by user loaded.")
    } else{
      message("Watershed stressor data provided by user not found.")
    }

  } else if(boo_ws == TRUE & boo_helper == TRUE){
    conusStates <- setdiff(datasets::state.name, c("Alaska", "Hawaii"))
    regionAvailable <- region %in% conusStates

    if(regionAvailable == TRUE){
      message("Downloading watershed stressor data from helper package.")

      data_stressorWS <- CASToolHelperPckg::getWSStressorData(region)
      data_stressorinfoWS <- CASToolHelperPckg::getWSStressorInfo()

      if("comid" %in% names(data_stressorWS)){
        data_stressorWS <- data_stressorWS |> dplyr::rename("COMID" = "comid")
      }

    } else{
      message("Watershed stressor data not available for the specified region. ")
    }
  }

  return(list(data_stressorWS = data_stressorWS,
              data_stressorinfoWS = data_stressorinfoWS))
}
