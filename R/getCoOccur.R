#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Co-Occurrence Line of Evidence
#'
#' @description Compares target sample values to unimpaired, comparator samples.
#' Generates line of evidence scores and boxplots.
#' Identifies stressors to remove from further consideration.
#'
#' @param TargetSiteID target site identifier
#' @param df_data data frame with paired stressor/response data, where stressor
#'                data are transformed.
#' @param detects stressors detected at the target site with sufficient comparator observations
#' @param df_stressinfo dataframe containing stressor metadata
#' @param biocomm current biological community
#' @param config list containing configuration parameters
#' @param colBio df_data column with biological index numeric value.
#' @param pHlimLow The lower limit of pH considered to be supportive of a
#'                 biological community. Defaults to pH of 6.5
#' @param pHlimHigh The upper limit of pH considered to be supportive of a
#'                  biological community. Defaults to pH of 9.
#' @param DOlim The lower limit of DO in mg/L considered to be supportive of a
#'              biological community. Defaults to 7 mg/L
#' @param plotvars Colors, fills, shapes, transparencies for each type (target,
#'                 not degraded, degraded, inside-the-case, outside-the-case).
#' @param plotdpi DPI for plots for standardization. Default = plot_dpi.  600
#' @param plotH Plot height for standardization. Default = plot_H. 6
#' @param plotW Plot width for standardization. Default = plot_W. 8
#' @param plotunits Plot units for standardization. Default = plot_units. "in"
#' @param dir_plots Directory to save plots.  Default = working directory and Results.
#' @param dir_sub Subdirectory for weight of evidence tables. Default = "_WoE"
#' @param boo_plot Boolean value to save plots.  Default = TRUE.
#' @param incaseLabel name of in-case designation (e.g., ClusterID)
#' @param targetSampleLabels boolean to include target sample response labels on figures
#'
#' @return Writes individual plots as pngs. Writes co-occurrence scores as csv.
#' Returns a list with filtered stressor metadata, stressors that did not pass the
#' co-occurence line of evidence screen, dataframe of co-occurrence scores, and gap dataframe.
#'
#' @export
getCoOccur <- function(TargetSiteID,
                       df_data,
                       detects,
                       df_stressinfo,
                       biocomm,
                       config = NULL,
                       colBio,
                       pHlimLow = 6.5,
                       pHlimHigh = 9,
                       DOlim = 7,
                       plotvars,
                       plotdpi = 600,
                       plotH = 6,
                       plotW = 8,
                       plotunits = "in",
                       dir_plots,
                       dir_sub = "_WoE",
                       boo_plot = TRUE,
                       incaseLabel = NULL,
                       targetSampleLabels = FALSE) {

  boo_DEBUG <- FALSE

  if (boo_DEBUG==TRUE) {
    TargetSiteID = TargetSiteID
    df_data = df_PairedStressResp
    detects = siteDetectsAll
    df_stressinfo = data_stressInfo
    biocomm = bioComm
    colBio = bioIndex
    pHlimLow = pHlimLow
    pHlimHigh = pHlimHigh
    DOlim = DOlim
    plotvars = plotvars
    plotdpi = 600
    plotH = 6
    plotW = 8
    plotunits = "in"
    dir_plots = dir_plots
    dir_sub = "_WoE"
    boo_plot = TRUE
  }

  if(is.null(config) == FALSE){
    plotvars      <- config$data_plotvars
    plotdpi       <- config$plot_dpi
    plotH         <- config$plot_H
    plotW         <- config$plot_W
    plotunits     <- config$plot_units
    boo_plot      <- config$boo_plot
  }

  biocomm <- toupper(biocomm)

  # Write results directory ----
  out.dir <- dirname(dir_plots)
  out.folders <- c(out.dir, basename(dir_plots), TargetSiteID, biocomm, dir_sub)

  for (i in 1:length(out.folders)) {
    if (i == 1) {
      dir_path <- file.path(out.folders[i])
    } else {
      dir_path <- file.path(dir_path, out.folders[i])
    }
    if (dir.exists(dir_path) == FALSE) {
      dir.create(dir_path)
    }
  }

  # Initialize gaps df
  df_gap <- data.frame(fxnname = character(), condition = character(), result = character(), comment = character())

  ## Plot, Variables
  # Define generic plot variables ----
  plotvars  <- plotvars |>
    dplyr::filter(Type %in% c("target", "insideND", "insideD"))
  bio_col     <- unlist(plotvars$Fill)
  bio_fill    <- unlist(plotvars$Fill)
  bio_shape   <- unlist(plotvars$Shape)
  bio_size    <- unlist(plotvars$Size)
  bio_alpha   <- unlist(plotvars$Alpha)

  ## Plot, Variables, Target Site Line
  targ_line_col <- bio_col[1]
  targ_line_lty <- 2
  targ_line_lwd <- 0.5

  # arrow labels
  aLabPos <- "1"
  aLabZero <- "0"
  aLabNeg <- "-1"

  ## Limit line types (dotted)
  lim_lty <- 3

  # Comparator site, not degraded data
  df.comp <- df_data |>
    dplyr::filter(IncaseYN == 1,
                  Quality == "Not degraded",
                  StationID != TargetSiteID)

  # Calculate comparator data statistics
  df.stats <- df.comp |>
    dplyr::select(dplyr::all_of(detects)) |>
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Stressor",
                        values_to = "StressorValue") |>
    dplyr::mutate(NotNA = ifelse(!is.na(StressorValue), 1, 0)) |>
    dplyr::group_by(Stressor) |>
    dplyr::summarize(n = sum(NotNA, na.rm = TRUE),
                     q25 = stats::quantile(StressorValue, probs = 0.25, na.rm = TRUE),
                     q50 = stats::quantile(StressorValue, probs = 0.50, na.rm = TRUE),
                     q75 = stats::quantile(StressorValue, probs = 0.75, na.rm = TRUE),
                     minVal = min(StressorValue, na.rm = TRUE),
                     maxVal = max(StressorValue, na.rm = TRUE))

  df_stressinfo <- df_stressinfo |>
    dplyr::rename("Stressor" = "StdParamName") |>
    dplyr::mutate(Label = dplyr::if_else(LogTransf == 1, paste("Log1p", Label), Label))

  df.stats <- df.stats |>
    dplyr::left_join(df_stressinfo |>
                       dplyr::select(Stressor, LogTransf, DirIncStress, Label),
                     by = "Stressor")

  # Score target site samples
  df.scores <- df_data |>
    dplyr::filter(StationID == TargetSiteID) |>
    tidyr::pivot_longer(cols = all_of(detects), names_to = "Stressor", values_to = "StressorValue") |>
    dplyr::left_join(df.stats, by = "Stressor") |>
    dplyr::mutate(Sc_Box = dplyr::case_when(
      (Stressor == "pH_low" & StressorValue < pHlimLow) ~ 1,
      (Stressor == "pH_high" & StressorValue > pHlimHigh) ~ 1,
      (stringr::str_detect(Stressor |> tolower() |> gsub("\\s", "", x = _), "^do$|dissolvedoxygen") & StressorValue < DOlim) ~ 1,
      (DirIncStress == "Dec") ~ dplyr::case_when(StressorValue > q50 ~ -1,
                                                 StressorValue < q25 ~ 1,
                                                 TRUE ~0),
      (DirIncStress == "Inc") ~ dplyr::case_when(StressorValue > q75 ~ 1,
                                                 StressorValue < q50 ~ -1,
                                                 TRUE ~ 0)),
      bioComm = biocomm) |>
    dplyr::select(-dplyr::any_of(c("OutcaseCol", "BioComm", "RefSiteFlag", "IncaseYN", "OutcaseYN", "BetterThan", "LogTransf", "DirIncStress")))

  i.Group <- unique(df.scores$IncaseCol)
  j.len <- length(detects)

  # Loop over detected stressors ####
  for (j in seq_along(detects)) {##FOR.j.START

    stressname <- detects[j]

    # Write graphics directory ----
    out.dir <- dirname(dir_plots)
    out.folders <- c(out.dir, basename(dir_plots), TargetSiteID, biocomm, stressname)

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

    message(paste0("Processing stressor (", j, "/", j.len, ") ",
                   stressname, ".\n"))
    utils::flush.console()

    # Subset dataframes by current stressor
    df.scores.stress <- df.scores |>
      dplyr::filter(Stressor == stressname) |>
      dplyr::mutate(case = "Target sample(s)")

    df.comp.stress <- df.comp |>
      dplyr::select(dplyr::all_of(stressname), Quality) |>
      dplyr::rename(Stressor = stressname) |>
      tidyr::drop_na()

    df.stats.stress <- df.stats |>
      dplyr::filter(Stressor == stressname)


    # File Name
    fn_png_p1 <- paste0(TargetSiteID, "_", make.names(stressname), "_",
                        biocomm, "_", colBio, "_CO.png")

    # Plot labels
    scores <- df.scores.stress |>
      dplyr::pull(Sc_Box)

    scores.text <- ifelse(is.na(scores), "NE",
                          ifelse(scores == 1, "1 (Supporting)",
                                 ifelse(scores == 0, "0 (Indeterminate)",
                                        ifelse(scores == -1, "-1 (Refuting)", "NA"))))

    lab.Score <- paste0("Score = ", paste0(scores.text, collapse = ", "))
    lab.N <- paste0("n = ", nrow(df.comp.stress))
    lab_comp <- paste0("Inside-the-case samples selected from ", incaseLabel,
                       " = ", i.Group)

    if(df.stats.stress$DirIncStress == "Inc"){
      subtitleCO <- "Is the target sample stressor value elevated compared to unimpaired, comparator samples?"
    } else if(df.stats.stress$DirIncStress == "Dec"){
      subtitleCO <- "Is the target sample stressor value depressed compared to unimpaired, comparator samples?"
    }

    lab.sub <- paste0("Paired stressor/response samples considered not degraded ",
                      "from inside the case (", lab.N, ").\n", lab.Score, ".")

    # scoring lines
    maxVal <- df.stats.stress |> dplyr::pull(maxVal)
    minVal <- df.stats.stress |> dplyr::pull(minVal)

    if (df.stats.stress$DirIncStress == "Dec") {
      box_qHI <- df.stats.stress |> dplyr::pull(q50)
      box_qLO <- df.stats.stress |> dplyr::pull(q25)
      segNeg <- ((maxVal - box_qHI) / 2) + box_qHI
      segZero <- ((box_qHI - box_qLO) / 2) + box_qLO
      segPos <- ((box_qLO - minVal) / 2) + minVal
    } else {
      box_qHI <- df.stats.stress |> dplyr::pull(q75)
      box_qLO <- df.stats.stress |> dplyr::pull(q50)
      segPos <- ((maxVal - box_qHI) / 2) + box_qHI
      segZero <- ((box_qHI - box_qLO) / 2) + box_qLO
      segNeg <- ((box_qLO - minVal) / 2) + minVal
    }

    # plot1, ggplot ####
    p1 <- ggplot2::ggplot(df.comp.stress, ggplot2::aes(y = Stressor, x = 1, group = 1)) +
      ggplot2::geom_boxplot(outliers = FALSE,
                            staplewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(data = df.scores.stress, ggplot2::aes(yintercept = StressorValue, lty = case), color = targ_line_col,
                          lwd = targ_line_lwd, na.rm = TRUE, show.legend = FALSE)+
      ggplot2::geom_vline(data = data.frame(xintercept = 1 - 0.5, lab = "Target sample(s)"),
                          ggplot2::aes(xintercept = xintercept, linetype = lab), inherit.aes = FALSE, alpha = 0) +
      ggplot2::guides(linetype = ggplot2::guide_legend(order = 1,
                                                       override.aes = list( alpha = 1, colour = targ_line_col, linewidth = targ_line_lwd ))) +
      ggplot2::scale_linetype_manual(name = "", values = targ_line_lty)+
      ggplot2::geom_jitter(ggplot2::aes(color = Quality, shape = Quality,
                                        fill = Quality), alpha = 0.5,
                           na.rm = TRUE, width = 0.25, height = 0.01) +
      ggplot2::annotate("segment", y = minVal, yend = box_qLO, x = c(1.5, 1.5),
                        color = "orange", linewidth = 0.7, alpha = 0.6,
                        arrow = grid::arrow(ends = "both", type = "open",
                                            length = grid::unit(0.2, "cm"))) +
      ggplot2::annotate("segment", y = box_qLO, yend = box_qHI, x = c(1.5, 1.5),
                        color = "orange", linewidth = 0.7, alpha = 0.6,
                        arrow = grid::arrow(ends = "both", type = "open",
                                            length = grid::unit(0.2, "cm"))) +
      ggplot2::annotate("segment", y = box_qHI, yend = maxVal, x = c(1.5, 1.5),
                        color = "orange", linewidth = 0.7, alpha = 0.6,
                        arrow = grid::arrow(ends = "both", type = "open",
                                            length = grid::unit(0.2, "cm"))) +
      ggplot2::annotate("text", x = 1.53, y = c(segNeg, segZero, segPos),
                        label = c(aLabNeg, aLabZero, aLabPos), color = "orange") +
      ggplot2::scale_color_manual(name = 'Samples',
                                  breaks = c("Degraded", "Not degraded"),
                                  values = bio_fill, drop = TRUE) +
      ggplot2::scale_fill_manual(name = "Samples",
                                 breaks = c("Degraded", "Not degraded"),
                                 values = bio_fill, drop = TRUE) +
      ggplot2::scale_shape_manual(name = "Samples",
                                  breaks = c("Degraded", "Not degraded"),
                                  values = bio_shape, drop = TRUE) +
      ggplot2::labs(title = paste0(TargetSiteID, ": Co-occurrence line of evidence"), subtitle = subtitleCO,
                    caption = lab.sub, y = df.stats.stress$Label, x = lab_comp) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())

    # Add target sample labels if requested
    if(targetSampleLabels == TRUE){
      p1 <- p1 +
        ggrepel::geom_label_repel(mapping = ggplot2::aes(x = 0.5, y = StressorValue, label = RespSampleID, group = NA), data = df.scores.stress, color = targ_line_col,size = 8/ggplot2::.pt, inherit.aes = FALSE)
    }

    ggplot2::ggsave(filename = file.path(dir_path_stress, fn_png_p1), plot = p1,
                    dpi = plotdpi, width = plotW, height = plotH, units = plotunits)

  }##FOR.j.END

  # Write scores file ####
  fn.scores <- file.path(dir_path, paste0(TargetSiteID, "_", biocomm, "_CO_Scores.csv"))
  write.csv(df.scores, file = fn.scores, row.names = FALSE)

  # Write gaps
  ## Defining notstressors in two steps is necessary to retain stressors with at least one score 0 or 1
  stressors <- unique(df.scores$Stressor[df.scores$Sc_Box != -1])
  notstressors <- unique(df.scores$Stressor[df.scores$Sc_Box == -1])
  notstressors <- setdiff(notstressors, stressors)

  if (length(notstressors) > 0) {
    for (s in seq_along(notstressors)) {
      msg <- paste0(notstressors[s], " identified as not a candidate cause for ",
                    TargetSiteID," for the ", bioComm, " community.")
      message(msg)

      # No identified stressors may be a data gap, but may not be, either
      gap.statement <- data.frame(
        fxnname = "getCoOccur",
        condition = paste0(notstressors[s], " score for ", bioComm,
                           " refutes"),
        result = as.character(-1),
        comment = msg)

      df_gap <- df_gap |>
        dplyr::bind_rows(gap.statement)
    }
  }

  # Prep df.scores for export to include in df_LoEs
  df.scores <- df.scores |>
    dplyr::filter(Stressor %in% stressors) |>
    dplyr::select(!Stressor) |>
    dplyr::rename(bioIndex = dplyr::all_of(colBio), Score = Sc_Box,
                  Stressor = Label) |>
    dplyr::mutate(LoE = "CO", bioIndexName = colBio) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, bioComm, bioIndexName, bioIndex, Quality,
                  Stressor, StressorValue, LoE, Score)

  # Filtered stressor metadata
  sstv.name <- paste0("SSTVname.", tolower(biocomm))
  sens.max <- paste0("SensMax.", tolower(biocomm))
  sens.min <- paste0("SensMin.", tolower(biocomm))

  df.stressorMetadata <- df_stressinfo |>
    dplyr::filter(Stressor %in% stressors) |>
    dplyr::select(Stressor, LogTransf, DirIncStress, Label, SSIndex,
                  dplyr::all_of(sstv.name), dplyr::all_of(sens.max),
                  dplyr::all_of(sens.min))

  return(list(df_stressorMetadata = df.stressorMetadata,
              notEvaluated = notstressors,
              df_COscores = df.scores,
              df_gap = df_gap))

}##FUNCTION.END
