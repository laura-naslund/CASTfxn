#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Weight-of-Evidence summary
#'
#' @description Summarize weight of evidence scores
#'
#' Uses the packages dplyr and tidyr.
#'
#' @param TargetSiteID target side identifier
#' @param biocomm current biological community
#' @param dfLoE dataframe containing the lines of evidence scores
#' @param dfStress metadata for stressors passing the co-occurrence screen
#' @param dir_results directory to save results
#' @param dir_WoE subdirectory for weight of evidence outputs
#' @param plotdpi plot resolution
#' @param plotH plot height
#' @param plotW plot width
#' @param plotunits plot units
#' @param boo_plot boolean for generating plots
#'
#' @return No returns. Writes a weight of evidence table with scores for each line of evidence
#' and a summary weight of evidence table which counts the number of lines of evidence
#' supporting, refuting, indeterminate and not evaluated for each target site
#' stressor-response pair. Also writes a bar chart depicting the summary weight of evidence table.
#'
#' @keywords internal
#' @export
getWoE <- function(TargetSiteID,
                   biocomm,
                   dfLoE,
                   dfStress,
                   dir_results,
                   dir_WoE = "_WoE",
                   config = NULL,
                   plotdpi,
                   plotH,
                   plotW,
                   plotunits,
                   boo_plot = TRUE) {##FUNCTION.START

  # QC data
  boo_DEBUG <- FALSE

  if (boo_DEBUG == TRUE) {
    TargetSiteID = TargetSiteID
    biocomm = bioComm
    dfLoE = df_LoE
    dfStress = df_stressorMetadata
    dir_results = dir_results
    dir_WoE = "_WoE"
  }

  if(is.null(config) == FALSE){
    plotdpi <- config$plot_dpi
    plotH <- config$plot_H
    plotW <- config$plot_W
    plotunits <- config$plot_units
    boo_plot <- config$boo_plot
  }

  # define pipe
  biocomm <- toupper(biocomm)

  # Write results directory ----
  out.dir <- dirname(dir_results)
  out.folders <- c(out.dir, basename(dir_results), TargetSiteID, biocomm,
                   dir_WoE)

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

  # Create summary line of evidence table, rolling up supporting, refuting, indeterminate, and not evaluated
  dfLoE_summary <- dfLoE |>
    dplyr::mutate(Score = dplyr::if_else(Score == "NE", NA, Score),
                  Score = suppressWarnings(as.numeric(Score)),
                  StressorValue = round(StressorValue, 3),
                  bioComm = toupper(bioComm)) |>
    dplyr::group_by(StationID, StressSampleID, StressSampleDate,
                    RespSampleID, RespSampleDate, bioComm, Stressor,
                    StressorValue) |>
    dplyr::summarize(NumSupport = sum(Score > 0, na.rm = TRUE),
                     NumRefute = sum(Score < 0, na.rm = TRUE),
                     NumIndeterminate = sum(Score == 0, na.rm = TRUE),
                     NumNotEvaluated = sum(is.na(Score)),
                     .groups = "drop") |>
    dplyr::arrange(StationID, Stressor, StressSampleDate)

  write.csv(dfLoE_summary,
              file.path(dir.path, paste0(TargetSiteID, "_LoESummary.csv")), row.names = FALSE)

  # Generate summary line of evidence plot
  dfLoESummary_pivot <- dfLoE_summary |>
    dplyr::select(Stressor, StressSampleID, StressSampleDate, NumSupport,
                  NumRefute, NumIndeterminate) |>
    dplyr::rename("Support" = "NumSupport", "Refute" = "NumRefute", "Indeterminate" = "NumIndeterminate") |>
    tidyr::pivot_longer(cols = c("Support", "Refute", "Indeterminate"), names_to = "Evidence", values_to = "Count") |>
    dplyr::mutate(Evidence = forcats::fct_relevel(as.factor(Evidence), "Support", "Refute", "Indeterminate")) |>
    dplyr::mutate(
      StressSampleID = as.character(StressSampleID),
      Stressor = stringr::str_wrap(Stressor, width = 40),
      StressSampleID = dplyr::if_else(nchar(StressSampleID)>20,
                                      paste0(stringr::str_sub(StressSampleID, 1, ceiling(nchar(StressSampleID)/2)), "\n",
                                             stringr::str_sub(StressSampleID, ceiling(nchar(StressSampleID)/2) + 1, nchar(StressSampleID)), "\n"),
                                      StressSampleID))

  stressor_order <- dfLoESummary_pivot |>
    dplyr::filter(Evidence == "Support") |>
    dplyr::group_by(Stressor) |>
    dplyr::summarize(max_count = max(Count)) |>
    dplyr::arrange(max_count) |>
    dplyr::pull(Stressor)

  sample_order <- dfLoESummary_pivot |>
    dplyr::distinct(StressSampleID, StressSampleDate) |>
    dplyr::arrange(dplyr::desc(StressSampleDate)) |>
    dplyr::pull(StressSampleID)


  dfLoESummary_pivot <- dfLoESummary_pivot |>
    dplyr::mutate(
      Stressor = as.factor(Stressor),
      Stressor = forcats::fct_relevel(Stressor, stressor_order),
      StressSampleID = as.factor(StressSampleID),
      StressSampleID = forcats::fct_relevel(StressSampleID, sample_order)
      )

  p1 <- ggplot2::ggplot(dfLoESummary_pivot,
                        ggplot2::aes(x = Stressor, y = Count, fill = StressSampleID))+
    ggplot2::geom_bar(width = 0.5, stat = "identity", position = "dodge")+
    ggplot2::scale_fill_brewer(palette = "Paired")+
    ggplot2::coord_flip()+
    ggplot2::facet_wrap(~Evidence, ncol = 3)+
    ggplot2::theme_bw()+
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))+
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 10),
                   strip.text = ggplot2::element_text(size = 11), axis.title = ggplot2::element_text(size = 12))

  fn_png_p1 <- paste0(TargetSiteID, "_",
                      biocomm, "_LoESummaryFig.png")

  if ((boo_plot) == TRUE) {
    ggplot2::ggsave(filename = file.path(dir.path, fn_png_p1),
                    plot = p1, dpi = plotdpi, width = plotW, height = plotH,
                    units = plotunits)
  }

  dfLoE <- dfLoE |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, bioComm, Stressor, StressorValue, LoE, Score) |>
    dplyr::mutate(StressorValue = round(StressorValue, 3),
                  bioComm = toupper(bioComm)) |>
    tidyr::pivot_wider(id_cols = c(StationID, StressSampleID, StressSampleDate,
                                   RespSampleID, RespSampleDate, bioComm,
                                   Stressor, StressorValue), names_from = LoE,
                       values_from = Score) |>
    dplyr::arrange(StationID, Stressor, StressSampleDate)

  write.csv(dfLoE, file.path(dir.path, paste0(TargetSiteID, "_LoEs.csv")), row.names = FALSE)

}

