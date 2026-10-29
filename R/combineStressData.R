#' Combine stress data
#'
#' @param boo_ts boolean indicating the presence of timestamped stressor data
#' @param boo_und boolean indicating the presence of undated stressor data
#' @param boo_outliers boolean indicating whether to remove outliers
#' @param ts_stress timestamped stressor data
#' @param und_stress undated stressor data
#'
#' @return A list with combined stressor data (data_stress), metadata (data_stressInfo), and outliers table (data_stressOutliers)
#' @export

combineStressData <- function(boo_ts = NULL,
                              boo_und = NULL,
                              boo_outliers = NULL,
                              ts_stress,
                              und_stress){

  if (isTRUE(boo_ts) && isTRUE(boo_und)) {
    # Combine metadata
    data_stressInfo <- rbind(ts_stress[["data_tsInfo"]], und_stress[["data_undInfo"]])

    # Combine data
    data_stress     <- rbind(ts_stress[["data_tsRaw"]], und_stress[["data_undRaw"]])

    # Combine outliers
    if (isTRUE(boo_outliers)) {
      data_stressOutliers <- rbind(ts_stress[["data_tsoutliers"]], und_stress[["data_undoutliers"]])
    }

  } else if (isTRUE(boo_ts)) {
    data_stressInfo <- ts_stress[["data_tsInfo"]]
    data_stress     <- ts_stress[["data_tsRaw"]]

    if (isTRUE(boo_outliers)) {
      data_stressOutliers <- ts_stress[["data_tsoutliers"]]
    } else {
      data_stressOutliers <- NULL
    }

  } else {
    data_stressInfo <- und_stress[["data_undInfo"]]
    data_stress     <- und_stress[["data_undRaw"]]

    if (isTRUE(boo_outliers)) {
      data_stressOutliers <- und_stress[["data_undoutliers"]]
    } else {
      data_stressOutliers <- NULL
    }
  }

  return(list(data_stressInfo = data_stressInfo,
              data_stress = data_stress,
              data_stressOutliers = data_stressOutliers))
}
