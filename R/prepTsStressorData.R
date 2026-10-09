#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#
#  R version 4.4.3
#
#' Prepare timestamped stressor data
#'
#' @description Prepare timestamped stressor data for CASTool
#'
#' @details Reviews each uploaded file against both user-input data (either via the shiny app or contained in the CASToolMetadata.xlsx file) to evaluate 1) if the file exists, 2) if it contains the required columns (which requires that they be named correctly), 3) whether datatypes meet requirements.
#'
#' If any required files do not exist, return the list of missing files and shut down.
#' If all required files exist, proceed to internal checks of the columns and data.
#' Lastly, perform minimal joins to determine missing values.
#'
#' dir_data and dir_results should be absolute and not relative paths.
#' The function `normalizePath` can be used to convert from relative to absolute path.
#'
#' @param in_dir x
#' @param out_dir directory where all rds files are stored
#' @param fn_data x
#' @param fn_meta x
#' @param removeOutliers x
#' @param sub_dir x
#'
#' @return A list containing data_Sites and data_cluster to be used in the CASTool.
#'
#' @examples
#' # None at this time
#' @export
prepTsStressorData <- function(config = NULL,
                                 # meta = NULL,
                                 in_dir = NULL,
                                 out_dir = NULL,
                                 fn_data,
                                 fn_meta,
                                 removeOutliers = NULL,
                                 sub_dir = "_Histograms") {
  # Global Bindings
  UseInStressorID <- StdParamName <- StationID <- StressSampleID <-
    StressSampleDate <- TransfResult <- IQRmethod <- SDmethod <- Outlier <- NULL

  # Debug
  boo.debug = TRUE
  if (boo.debug) {
    sub_dir = "_Histograms"
  }

  if(is.null(config) == FALSE){
    in_dir <- file.path(config$out_dir, config$dn_checked_sk)
    out_dir <- config$out_dir
  }

  # if(is.null(meta) == FALSE){
  #   removeOutliers <- meta$removeOutliers
  # }

  # Create output folder
  out.folders <- c(out_dir, sub_dir)

  for (i in 1:length(out.folders)) {
    if (i == 1) {
      dir.path <- file.path(out.folders[i])
    } else {
      dir.path <- file.path(dir.path, out.folders[i])
    }
    if (dir.exists(dir.path) == FALSE) {
      dir.create(dir.path)
    }
  }

  out_dir <- dir.path


  # Load RDS files
  data_tsInfo   <- readRDS(file.path(in_dir, fn_meta))
  data_tsAll    <- readRDS(file.path(in_dir, fn_data))

  ## Get metadata for all tsured stressors
  data_tsInfo   <- data_tsInfo |>
    dplyr::filter(UseInStressorID == 1)

  ## Get tsured stressor values to use in the stressor id
  params2use      <- as.character(data_tsInfo$StdParamName)
  data_tsAll    <- dplyr::filter(data_tsAll, StdParamName %in% params2use)

  ## getOutliers ----
  ## returns a dataframe with SampleID, StdParamName, TransfResult,
  ## IQRmethod, SDmethod, Outlier
  data_tsOutliers <- getOutliers(df_data   = data_tsAll,
                                   df_meta   = data_tsInfo,
                                   dir_plots = out_dir)
  ## Merge outlier flags with raw data by sample ID
  data_tsAll <- merge(data_tsAll, data_tsOutliers,
                        by.x = c("StationID", "StressSampleID", "StressSampleDate",
                                 "StdParamName", "ResultValue"),
                        by.y = c("StationID", "StressSampleID", "StressSampleDate",
                                 "StdParamName", "ResultValue"),
                        all.x = TRUE)
  data_tsAll <- data_tsAll |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, StdParamName,
                  TransfResult, IQRmethod, SDmethod, Outlier)
  # Clean up
  rm(data_tsOutliers, params2use)

  ## Average duplicate data -- remove outliers first, if desired
  if (!is.na(removeOutliers) & removeOutliers == TRUE) {
    data_tsoutliers <- data_tsAll |>
      dplyr::filter(Outlier == "Outlier")
    data_tsAll <- data_tsAll |>
      dplyr::filter(Outlier %in% c("Good", "NE"))
  } else {
    data_tsoutliers <- NULL
  }

  data_tsAll <- data_tsAll |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, StdParamName,
                  TransfResult) |>
    dplyr::mutate(StressSampleDate = lubridate::parse_date_time(StressSampleDate,
                                              orders = c("ymd", "mdy", "dmy")) |>
                    lubridate::date(),
                  # StressSampleID = stringr::str_replace_all(StressSampleID, "[:punct:]", "_"),
                  # StationID = stringr::str_replace_all(StationID, "[:punct:]", "_")
                  ) |>
    dplyr::group_by(StationID, StressSampleID, StressSampleDate, StdParamName) |>
    dplyr::mutate(TransfResult = mean(TransfResult, na.rm = TRUE)) |>
    dplyr::filter(!is.na(TransfResult))
  data_tsAll_dups <- data_tsAll[duplicated(data_tsAll),]
  data_tsAll <- unique(data_tsAll) # should be unique, long-form sample/analyte

  # Duplicate all pH values, one to use for high pH as a stressor and one for low pH as a stressor
  data_pHHigh <- data_tsAll |>
    dplyr::filter(stringr::str_detect(tolower(StdParamName), "^ph$") == TRUE) |>
    dplyr::mutate(StdParamName = "pH_high")
  data_tsAll <- data_tsAll |>
    dplyr::mutate(StdParamName = ifelse(stringr::str_detect(tolower(StdParamName), "^ph$") == TRUE, "pH_low",
                                        StdParamName))
  data_tsRaw <- rbind(data_tsAll, data_pHHigh)
  rm(data_pHHigh)

  # Replace metadata pH with pH_low and pH_high
  pH_meta <- data_tsInfo |>
    dplyr::filter(stringr::str_detect(tolower(StdParamName), "^ph$"))

  if(nrow(pH_meta)>0){
    pH_metaHigh <- pH_meta |> dplyr::mutate(StdParamName = "pH_high", Label = "pH (high)", DirIncStress = "Inc")
    pH_metaLow <- pH_meta |> dplyr::mutate(StdParamName = "pH_low", Label = "pH (low)", DirIncStress = "Dec")

    data_tsInfo <- data_tsInfo |>
      dplyr::filter(stringr::str_detect(tolower(StdParamName), "^ph$") == FALSE) |>
      dplyr::bind_rows(pH_metaHigh) |>
      dplyr::bind_rows(pH_metaLow)
  }

  # Adjust dates for samples collected on multiple dates with the same ID
  analytes <- unique(as.character(data_tsRaw$StdParamName))
  data_tsRaw <- as.data.frame(data_tsRaw) |>
    dplyr::filter(!is.na(TransfResult)) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate,
                  StdParamName, TransfResult) |>
    tidyr::pivot_wider(names_from = StdParamName,
                       values_from = TransfResult)
  data_tsRaw <- data_tsRaw |>
    dplyr::group_by(StationID, StressSampleID) |>
    dplyr::mutate(StressSampleDate = min(StressSampleDate)) |>
    dplyr::ungroup()
  data_tsRaw <- data_tsRaw |>
    tidyr::pivot_longer(cols = dplyr::all_of(analytes),
                        names_to = "StdParamName",
                        values_to = "TransfResult")
  data_tsRaw <- dplyr::filter(data_tsRaw, !is.na(TransfResult))

  mytsData <- list(data_tsInfo = data_tsInfo,
                     data_tsRaw = data_tsRaw,
                     data_tsoutliers = data_tsoutliers)

  return(mytsData)

}
