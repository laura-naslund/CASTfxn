#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.2
#
#' @title Write outlier gap statements
#'
#' @description Writes gap statements for identified outliers for target site, inside-the-case sites,
#' and outside-the-case sites.
#'
#' Uses the library dplyr.
#'
#' @param TargetSiteID Target site identifier
#' @param df_outliers Dataframe containing any stressor values identified outliers.
#' @param df_stressInfo Dataframe containing stressor metadata.
#' @param df_Sites Dataframe containing site information
#' @param siteDetectsAll All stressors ever detected in any target site samples.
#' @param compSites Vector containing comparator site identifiers (necessary if useBC = TRUE).
#' @param useBC boolean to use biological similarity filter
#' @param allSites Vector containing all "outside the case" identifiers (necessary if useBC = TRUE)
#' @param dir_results directory to write outliers description to
#'
#' @keywords internal
#' @examples
#' # None at this time
#' @export
writeOutliers <- function(TargetSiteID,
                          df_outliers,
                          df_stressInfo,
                          df_Sites,
                          siteDetectsAll,
                          compSites = NULL,
                          useBC = FALSE,
                          allSites = NULL,
                          dir_results = file.path(getwd(), "Results")) {

  # Global Bindings
  data_outliers <- data_stressInfo <- list.CompSites <-
    StationID <- TransfResult <- Outlier <- StdParamName <- IncaseCol <-  NULL

  # Debug
  boo.debug = FALSE
  #
  if (boo.debug) {
    df_outliers = data_outliers
    df_stressInfo = data_stressInfo
    TargetSiteID = TargetSiteID
    siteDetects = siteDetectsAll
    compSites = list.CompSites$comp.sites
    allSites = list.CompSites$all.sites
    dir_results = dir_results
  }

  # Initialize gaps df
  df_gap <- data.frame(fxnname = character(), condition = character(), result = character(), comment = character())

  if(useBC == TRUE){
    compSites <- compSites
    allSites <- allSites
  } else{
    TargetSiteCluster <- df_Sites |>
      dplyr::filter(StationID == TargetSiteID) |>
      dplyr::pull(IncaseCol)

    compSites <- df_Sites |>
      dplyr::filter(IncaseCol == TargetSiteCluster) |>
      dplyr::pull(StationID)

    allSites <- df_Sites |>
      dplyr::filter(is.na(IncaseCol) | IncaseCol != TargetSiteCluster) |>
      dplyr::pull(StationID)
  }

  # add case identifier to outliers dataframe and add stressor parameter label
  df_outliers <- df_outliers |>
    dplyr::mutate(StationIDType = dplyr::case_when(
      StationID == TargetSiteID ~ "Target site",
      StationID %in% compSites ~ "Inside-the-case",
      TRUE ~ "Outside-the-case")) |>
    dplyr::left_join(df_stressInfo |> dplyr::select("StdParamName", "Label"), by = "StdParamName") |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, StationIDType, StdParamName, Label, TransfResult,
                  IQRmethod, SDmethod)

  file.outliers <- paste0(TargetSiteID, "_outliers.csv")
  fn.outliers <- file.path(dir_results, TargetSiteID, file.outliers)

  write.csv(df_outliers, fn.outliers, row.names=FALSE)

  n_target_outliers <- df_outliers |>
    dplyr::filter(StationIDType == "Target site") |>
    nrow()

  n_inside_outliers <- df_outliers |>
    dplyr::filter(StationIDType == "Inside-the-case") |>
    nrow()

  n_outside_outliers <- df_outliers |>
    dplyr::filter(StationIDType == "Outside-the-case") |>
    nrow()

  gap.statement.target <- data.frame(
    fxnname = "writeOutliers",
    condition = "Target site outliers",
    result = as.character(n_target_outliers),
    comment = paste0("See ", file.outliers, " in the ", TargetSiteID, " results folder for detailed description of outliers.")
  )

  gap.statement.inside <- data.frame(
    fxnname = "writeOutliers",
    condition = "Inside-the-case outliers",
    result = as.character(n_inside_outliers),
    comment = paste0("See ", file.outliers, " in the ", TargetSiteID, " results folder for detailed description of outliers.")
  )

  gap.statement.outside <- data.frame(
    fxnname = "writeOutliers",
    condition = "Outside-the-case outliers",
    result = as.character(n_outside_outliers),
    comment = paste0("See ", file.outliers, " in the ", TargetSiteID, " results folder for detailed description of outliers.")
  )

  df_gap <- df_gap |>
    dplyr::bind_rows(gap.statement.target) |>
    dplyr::bind_rows(gap.statement.inside) |>
    dplyr::bind_rows(gap.statement.outside)

  return(list(df_gap = df_gap))
}
