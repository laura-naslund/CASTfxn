#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Identifies Stressor Outliers
#'
#' @description Flags stressor outliers using two methods: 3 times IQR and
#'              6 times sd. If either suggest a given value as an outlier, the
#'              value is identified as an outlier.
#'
#' @details In addition to flagging likely outliers (by both 3x IQR and 6X sd),
#'          this function writes histograms of all untransformed stressors and
#'          for transformed stressors when transformation is specified in the
#'          stressor metadata. Writes histograms with resolution of 100 rather than
#'          default set in config to optimize the size of the results folder.
#'
#' Uses the dplyr library
#'
#' @param df_data dataframe containing all stressor data values
#' @param df_meta dataframe containing stressor metadata
#' @param dir_plots directory to which to write stressor histograms
#'
#' @return Dataframe containing sample ID, stressor name, stressor value,
#'         IQR method flag, SD method flag, Outlier flag
#'
# @keywords
#' @examples
#' # None at this time
#' @export
getOutliers <- function(df_data,
                        df_meta,
                        dir_plots) {##FUNCTION.START

  # Debug
  boo_DEBUG <- FALSE
  if (boo_DEBUG == TRUE) {
    df_data = data_Stress
    df_meta = data_stressInfo
    dir_plots = file.path(out.dir, "_Histograms")
  }

  # Ensure uniqueness of df_meta
  df_meta <- df_meta |>
    dplyr::select(StdParamName, LogTransf) |>
    dplyr::group_by(StdParamName) |>
    dplyr::summarise(LogTransf = max(LogTransf), .groups = "drop_last")

  # Merge data with metadata (note that some metadata may not be in the raw data)
  df_data <- merge(df_data, df_meta[, c("StdParamName", "LogTransf")],
                   by.x = "StdParamName", by.y = "StdParamName",
                   all.x = TRUE)
  df_data <- df_data |>
    dplyr::mutate(ResultValue = ifelse(ResultValue < 0 & LogTransf == 1,
                                       NA, ResultValue)) |>
    dplyr::filter(!is.na(ResultValue)) |>
    dplyr::mutate(TransfResult = ifelse(LogTransf == 1,
                                        suppressWarnings(log1p(ResultValue)),
                                        ResultValue))

  params <- unique(as.character(df_data$StdParamName))

  for (a in seq_along(params)) { # Iterate over parameters

    paramName <- params[a]
    msg <- paste0("Evaluating ", paramName)
    message(msg)

    df_sub <- df_data |>  # Subset for just one parameter at a time
      dplyr::filter(StdParamName == paramName)

    fn <- paste0(paramName, ".png")
    fn2 <- paste0(paramName, "_transf.png")

    # plot histogram of data
    png(file.path(dir_plots, fn), width = 8, height = 6, units = "in", res = 100)
    hist(df_sub$ResultValue, main = paste0("Histogram of ", paramName, " observations"), xlab = paramName)
    dev.off()


    if (unique(df_sub$LogTransf) == 1) {
      png(file.path(dir_plots, fn2), width = 8, height = 6, units = "in", res = 100)
      hist(df_sub$TransfResult, main = paste0("Histogram of Log1p transformed ", paramName, " observations"), xlab = paramName)
      dev.off()
    }

    # 3*IQR method for identifying outliers
    iqr <- stats::quantile(df_sub$TransfResult, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr1.5 <- 1.5 * stats::IQR(df_sub$TransfResult, na.rm = TRUE)
    outlowlim <- iqr[1] - iqr1.5
    outhilim <- iqr[2] + iqr1.5

    df_sub <- df_sub |>
      dplyr::mutate(IQRmethod = dplyr::case_when(dplyr::between(TransfResult,
                                                          outlowlim, outhilim) ~ "Good",
                                                 TransfResult < outlowlim ~ "Outlier low",
                                                 TransfResult > outhilim ~ "Outlier high",
                                                 TRUE ~ "NE"))

    # 6*sd method for identifying outliers
    paramMean <- mean(df_sub$TransfResult[is.finite(df_sub$TransfResult)],
                      na.rm = TRUE)
    paramSD <- stats::sd(df_sub$TransfResult[is.finite(df_sub$TransfResult)],
                         na.rm = TRUE)
    df_sub <- df_sub |>
      dplyr::mutate(SDmethod = dplyr::case_when((abs(TransfResult - paramMean) >
                                                    (6 * paramSD)) ~ "Outlier",
                                                 paramSD == 0 ~ "NE",
                                                 TRUE ~ "Good"))

    # Combine the findings
    # If both IQR & SD methods agree on "Good", then mark as "Good"
    # If either IQR method OR SD method aren't evaluated, then mark as "NE"
    # If either IQR method or SD method suggest an outlier, then mark as "Outlier"
    df_sub <- df_sub |>
      dplyr::mutate(IQRmethod = ifelse(!is.na(IQRmethod), IQRmethod, "NE"),
                    SDmethod = ifelse(!is.na(SDmethod), SDmethod, "NE")) |>
      dplyr::mutate(Outlier = dplyr::case_when((IQRmethod == "Good") & (SDmethod == "Good") ~ "Good",
                                               (IQRmethod == "NE") | (SDmethod == "NE") ~ "NE",
                                               grepl("Outlier", IQRmethod) | (SDmethod == "Outlier") ~
                                                 "Outlier",
                                               TRUE ~ "Good"))

    if (a == 1) {
      df_temp <- df_sub
    } else {
      df_temp <- rbind(df_temp, df_sub)
    }
  } # End parameter iteration

  myOutliers <- df_temp

  return(myOutliers)

}
