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

  ## Describe whether timestamped or undated data loaded
  if(stringr::str_detect(out$loaded, "tsInfo") |> any()){
    out$boo_ts <- TRUE
  } else{
    out$boo_ts <- FALSE
  }

  if(stringr::str_detect(out$loaded, "undInfo") |> any()){
    out$boo_und <- TRUE
  } else{
    out$boo_und <- FALSE
  }

  return(out)
}

