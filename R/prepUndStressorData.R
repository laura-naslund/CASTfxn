#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#
#  R version 4.4.3
#
#' @description Prepare undated stressor data for CASTool
#'
#' @details Prepare undated stressor data for CASTool. Filter dataset to parameters selected by the user. Identify and (if selected) remove outliers.
#'
#' @param in_dir directory containing checked CASTool data inputs
#' @param out_dir directory where data outputs are written
#' @param fn_data name of timestamped stressor rds data file
#' @param fn_meta name of timestamped stressor rds metadata file
#' @param removeOutliers boolean indicating whether to remove observations identified as outliers
#' @param sub_dir subdirectory where stressor histograms are written
#'
#' @return A list containing formatted undated stressor data (data_undRaw), metadata (data_undInfo), and identified outliers (data_undoutliers).
#'
#' @export
prepUndStressorData <- function(in_dir,
                                out_dir,
                                fn_data,
                                fn_meta,
                                config = NULL,
                                removeOutliers,
                                sub_dir = "_Histograms") {

  # Debug
  boo.debug = TRUE
  if (boo.debug) {
    sub_dir = "_Histograms"
  }

  if(is.null(config) == FALSE){
    in_dir <- file.path(config$out_dir, config$dn_checked_sk)
    out_dir <- config$out_dir
  }

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
  data_undInfo <- readRDS(file.path(in_dir, fn_meta))
  data_undAll  <- readRDS(file.path(in_dir, fn_data))

  ## Get metadata for all measured stressors
  data_undInfo   <- data_undInfo  |>
    dplyr::filter(UseInStressorID == 1)

  ## Get measured stressor values to use in the stressor id
  params2use    <- as.character(data_undInfo$StdParamName)
  data_undAll <- dplyr::filter(data_undAll, StdParamName %in% params2use)

  ## getOutliers returns a dataframe with undSampleID, StdParamName, TransfResult,
  ## IQRmethod, SDmethod, Outlier
  data_undOutliers <- getOutliers(df_data   = data_undAll,
                                    df_meta   = data_undInfo,
                                    dir_plots = out_dir)
  ## Merge outlier flags with raw data by sample ID
  data_undAll <- merge(data_undAll, data_measOutliers,
                        by.x = c("StationID", "StressSampleID", "StressSampleDate",
                                 "StdParamName", "ResultValue"),
                        by.y = c("StationID", "StressSampleID", "StressSampleDate",
                                 "StdParamName", "ResultValue"),
                        all.x = TRUE)
  data_undAll <- data_undAll  |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, StdParamName,
                  TransfResult, IQRmethod, SDmethod, Outlier)
  # Clean up
  rm(data_undOutliers, params2use)

  ## Average duplicate data -- remove outliers first, if desired
  if (removeOutliers) {
    data_undoutliers <- data_undAll  |>
      dplyr::filter(Outlier == "Outlier")
    data_undAll <- data_undAll  |>
      dplyr::filter(Outlier %in% c("Good", "NE"))
  } else {
    data_undoutliers <- NULL
  }

  data_undAll <- data_undAll  |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, StdParamName,
                  TransfResult)  |>
    dplyr::group_by(StationID, StressSampleID, StressSampleDate, StdParamName)  |>
    dplyr::mutate(TransfResult = mean(TransfResult, na.rm = TRUE))  |>
    dplyr::filter(!is.na(TransfResult))

  data_undAll <- unique(data_undAll) # should be unique, long-form sample/analyte

  myundData <- list(data_undInfo    = data_undInfo,
                     data_undRaw      = data_undAll |>
                       dplyr::mutate(StressSampleDate = NA) |>
                       dplyr::select(StationID, StressSampleID, StressSampleDate, TransfValue),
                     data_undoutliers = data_undoutliers)

  return(myundData)

}
