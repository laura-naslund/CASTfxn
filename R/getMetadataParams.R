#' Get metadata parameters
#'
#' @param out_dir path to results directory
#' @param dn_checked_sk name of checked data directory
#' @param config list containing the configuration parameters
#'
#' @return list with metadata parameters
#' @export
getMetadataParams <- function(out_dir = NULL,
                              dn_checked_sk = NULL,
                              config = NULL){

  if(is.null(config) == FALSE){
    out_dir <- config$out_dir
    dn_checked_sk <- config$dn_checked_sk
  }

  ## Load CASTool_Metadata ####
  data_CASTmeta <- readRDS(file.path(out_dir, dn_checked_sk, "CASTmetadata.rds"))

  ## Set correct classes
  out <- setNames(as.list(data_CASTmeta$Value), data_CASTmeta$Variable)

  out[c("DOlim",
        "pHlimLow",
        "pHlimHigh",
        "samplim",
        "r2_cutoff",
        "p.val_cutoff")] <- lapply(out[c("DOlim",
                                         "pHlimLow",
                                         "pHlimHigh",
                                         "samplim",
                                         "r2_cutoff",
                                         "p.val_cutoff")], as.numeric)

  out[c("helperImport",
        "exploreWSStressor",
        "calcRelAbund",
        "targetSampleLabels",
        "removeOutliers",
        "useAllCompReaches")] <- lapply(out[c("helperImport",
                                              "exploreWSStressor",
                                              "calcRelAbund",
                                              "targetSampleLabels",
                                              "removeOutliers",
                                              "useAllCompReaches")], as.logical)

  ## Make comma separated strings into lists
  out$biocommlist <- data_CASTmeta |>
    dplyr::filter(Variable == "biocommlist") |>
    dplyr::pull(Value) |>
    tolower() |>
    stringr::str_split(", |,") |>
    unlist()

  out$lagdays <- data_CASTmeta |>
    dplyr::filter(Variable == "lagdays") |>
    dplyr::pull(Value) |>
    stringr::str_split(",") |>
    unlist() |>
    stringr::str_trim() |>
    as.integer()

  ## Read loaded.rds ####
  ### Describes which input files have been loaded
  data_loaded <- readRDS(file.path(out_dir, dn_checked_sk, "loaded.rds"))
  out$loaded <- as.character(data_loaded$Object)

  ## Describe whether measured and/or modeled data loaded
  if(stringr::str_detect(out$loaded, "chem") |> any()){
    out$boo_meas <- TRUE
  } else{
    out$boo_meas <- FALSE
  }

  if(stringr::str_detect(out$loaded, "model") |> any()){
    out$boo_model <- TRUE
  } else{
    out$boo_model <- FALSE
  }

  return(out)
}

