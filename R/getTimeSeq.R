#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Time Sequence Line of Evidence (Graphics only)
#'
#' @description Graph time-specific stressor-response values.
#'
#' @details Generates faceted time sequence graphics (stressor/response, one
#' atop the other). All stressor/response data are graphed.
#' Entries for the line of evidence table only apply to paired data.
#'
#' Uses the packages dplyr, ggplot2, lubridate, tidyr
#'
#' @param TargetSiteID target site identifier
#' @param biocomm biological response community (bmi, algae, or fish)
#' @param config list containing configuration parameters
#' @param df_resp dataframe containing the target site biological response samples'
#' index and metric values
#' @param df_respinfo dataframe containing response metadata
#' @param df_stress dataframe containing all stressor data
#' @param df_stressinfo dataframe containing stressor metadata filtered to stressors that passed the co-occurrence screening
#' @param df_paired paired stressor-response data used only to generate scores table
#' @param plotdpi Standardized dpi value for graphics. Default = plot_dpi. 600
#' @param plotH Standardized height value for graphics. Default = plot_H. 6
#' @param plotW Standardized width value for graphics. Default = plot_W. 8
#' @param plotunits Standardized units for graphics. Default = plot_units. "in"
#' @param dir_results Directory containing results
#' @param dir_sub Subdirectory for weight of evidence outputs. Default = "_WoE"
#' @param boo_plot Boolean to save plots. Default = TRUE.
#'
#' @return One or more pngs in SiteID/Biocomm/Stressor subfolder of the
#'        "Results" folder. No scores are currently generated.
#'
#' @keywords internal
#' @export
getTimeSeq <- function(TargetSiteID,
                       biocomm,
                       config = NULL,
                       colBio,
                       df_stress,
                       df_resp,
                       df_respinfo,
                       df_stressinfo,
                       df_paired,
                       plotdpi = 600,
                       plotH = 6,
                       plotW = 8,
                       plotunits = "in",
                       dir_results = file.path(getwd(),"Results"),
                       dir_sub = "_WoE",
                       boo_plot = TRUE) {##FUNCTION.START

  # Debug
  boo_DEBUG <- FALSE

  if (boo_DEBUG == TRUE) {
    TargetSiteID
    biocomm = bioComm
    df_resp = bioMetricData
    df_respinfo = data_bmiMetricsInfo
    df_stress = list.Stress$data_stress
    df_stressinfo = df_stressorMetadata
    df_paired = df_PairedStressResp
    plotdpi = 600
    plotH = 6
    plotW = 8
    plotunits = "in"
    dir_results = dir_results
    dir_sub = "_WoE"
  }

  if(is.null(config) == FALSE){
    plotdpi       <- config$plot_dpi
    plotH         <- config$plot_H
    plotW         <- config$plot_W
    plotunits     <- config$plot_units
    boo_plot      <- config$boo_plot
  }

  biocomm <- toupper(biocomm)

  # Write results directory ----
  out.folders <- c(dir_results, TargetSiteID, biocomm, dir_sub)

  for (i in 1:length(out.folders)) {
    if (i == 1) {
      path <- file.path(out.folders[i])
    } else {
      path <- file.path(path, out.folders[i])
    }
    if (dir.exists(path) == FALSE) {
      dir.create(path)
    }
  }

  path <- file.path(dir_results, TargetSiteID, biocomm, dir_sub)

  # Get vector of target site stressors
  stressors <- df_stressinfo$Stressor
  metrics <- df_respinfo$MetricName

  # Subset stress and response dataframes to TargetSite
  df_stress <- df_stress |>
    dplyr::filter(StationID == TargetSiteID) |>
    dplyr::filter(is.na(TransfResult)==FALSE)

  df_resp <- df_resp |>
    dplyr::filter(StationID == TargetSiteID)

  count = 1
  totplots <- length(stressors) * length(metrics)

  # Loop over stressors
  msgNum <- 0
  for (s in seq_along(stressors)) {
    stressname <- stressors[s]
    stressLabel <- df_stressinfo$Label[df_stressinfo$Stressor == stressname]

    # Write graphics directory ----
    out.dir <- dirname(dir_results)
    out.folders <- c(out.dir, basename(dir_results), TargetSiteID, biocomm, stressname)

    for (i in 1:length(out.folders)) {
      if (i == 1) {
        dir_path_stress <- file.path(out.folders[i])
      } else {
        dir_path_stress <- file.path(dir_path_stress, out.folders[i])
      }
      if (dir.exists(dir_path_stress) == FALSE) {
        dir.create(dir_path_stress)
      }
    }

    # Loop over responses
    for (r in seq_along(metrics)) {
      metricname <- metrics[r]
      metricLabel <- df_respinfo$MetricLabel[df_respinfo$MetricName == metricname]
      msgNum <- msgNum + 1

      # Create filename for graphic
      fn = paste0(TargetSiteID, "_", make.names(stressname), "_",
                  biocomm, "_", make.names(metricname), "_TS.png")
      fpath = file.path(dir_path_stress, fn)

      # Subset dataframe for current stressor and response
      df.plotstress <- df_stress |>
        dplyr::filter(StdParamName == stressname) |>
        dplyr::rename(SampleDate = StressSampleDate, Value = TransfResult) |>
        dplyr::select(StationID, SampleDate, Value) |>
        dplyr::mutate(type = "Stressor", Value = signif(Value, digits = 3))

      df.plotresp <- df_resp |>
        dplyr::select(StationID, RespSampleDate, dplyr::all_of(metricname)) |>
        dplyr::rename(SampleDate = RespSampleDate, Value = {{metricname}}) |>
        dplyr::mutate(type = "Response", Value = signif(Value, digits = 3)) |>
        dplyr::filter(!is.na(Value))

      df.plot <- rbind(df.plotstress, df.plotresp) |>
        dplyr::mutate(type = factor(type, levels = c("Stressor", "Response"),
                                    labels = c(stressLabel, metricLabel)))

      # Get min/max date range
      minDate <- as.Date(min(df.plot$SampleDate) - 30)
      minYear <- lubridate::year(minDate)
      minDate <- lubridate::ymd(paste0(minYear, "/01/01"))
      maxDate <- as.Date(max(df.plot$SampleDate) + 30)
      maxYear <- lubridate::year(maxDate)
      maxDate <- lubridate::ymd(paste0(maxYear, "/12/31"))

      df.plot <- df.plot |>
        tidyr::complete(tidyr::nesting(type),
                        SampleDate = seq(minDate, maxDate, by = "day"))

      msg <- paste0("Plotting bar graphs (", msgNum, "/", totplots, ") ",
                    stressname, " and ", metricname)
      message(msg)

      p_ts <- ggplot2::ggplot(df.plot, ggplot2::aes(x = SampleDate, y = as.numeric(Value))) +
        ggplot2::geom_col(col = "black", fill = "black", linewidth = 0.2, alpha = 0.5)+
        ggplot2::geom_label(ggplot2::aes(x = SampleDate, y = as.numeric(Value), label = Value, vjust = 1))+
        ggplot2::facet_wrap(~type, ncol = 1, scales = "free_y")+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
                       panel.grid.minor = ggplot2::element_blank())+
        ggplot2::scale_x_date(limits = c(minDate, maxDate), date_labels = "%Y-%m-%d", date_breaks = "3 months")+
        ggplot2::labs(title = paste(TargetSiteID, "Time sequence line of evidence"),
                                   subtitle = "Does elevated stress precede the biological response?",
                                   x = "Sample Date", y = "Value")

      if(boo_plot){
        ggplot2::ggsave(filename = fpath, p_ts, dpi = plotdpi, width = plotW,
                        height = plotH, units = plotunits)
      }## IF ~ boo_plot ~ END

    } # END loop over responses & graphics

  } # END loop over stressors

  # Use paired data to prepare scores table
  df_TS_scores <-   df_paired |>
    dplyr::filter(StationID == TargetSiteID) |>
    dplyr::select(StationID, StressSampleDate, RespSampleDate, StressSampleID,
                  RespSampleID, Quality, dplyr::all_of(stressors)) |>
    dplyr::left_join(df_resp |> dplyr::select(StationID, RespSampleID,
                                              RespSampleDate, dplyr::all_of(metrics)),
                     by = c("StationID", "RespSampleID", "RespSampleDate")) |>
    tidyr::pivot_longer(cols = dplyr::all_of(stressors), names_to = "Stressor",
                        values_to = "StressorValue") |>
    tidyr::pivot_longer(cols = dplyr::all_of(metrics), names_to = "bioIndexName",
                        values_to = "bioIndex") |>
    dplyr::filter(bioIndexName == colBio, !is.na(StressorValue)) |>
    dplyr::mutate(bioComm = biocomm, LoE = "TS", Score = "NE") |>
    dplyr::left_join(df_stressinfo |> dplyr::select(Stressor, Label),
                     by = "Stressor") |>
    dplyr::select(-Stressor) |>
    dplyr::rename(Stressor = Label)

  return(df_TS_scores)

}
