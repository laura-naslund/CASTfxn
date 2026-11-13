#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Verified Predictions
#'
#' @description Conduct stressor specific index analysis
#'
#' @param TargetSiteID target site identifier
#' @param stressors.ssi vector of stressors identified as candidate causes with SSIs
#' @param info.stress.ssi stressor metadata filtered to those with SSIs
#' @param info.ssi metrics metadata filtered to SSIs
#' @param df_paired paired stressor biological response data
#' @param biocomm current biological community
#' @param df_bioMetricData metric data for the current biological community
#' @param config list of configuration parameters
#' @param plotvars plotting parameters
#' @param plotdpi standardized plot dpi. Default = 600
#' @param plotH standardized plot height. Default = 6
#' @param plotW standardized plot width. Default = 8
#' @param plotunits units for plot height and width. Default = "in"
#' @param dir_plots results directory
#' @param dir_sub name of directory to write scores to. Default = "_WoE"
#' @param boo_plot = TRUE
#' @param incaseLabel label designating a site inside the case
#' @param targetSampleLabels boolean to label target response values with sample ID
#'
#' @return Writes boxplot and logistic regression scores to the _WoE folder.
#' Generates associated figures in the corresponding stressor folders.
#' Returns scores dataframe.

#' @export
getVPSSI <- function(TargetSiteID,
                     stressors.ssi,
                     info.stress.ssi,
                     info.ssi,
                     df_paired,
                     biocomm,
                     df_bioMetricData,
                     config = cfg,
                     plotvars,
                     plotdpi = 600,
                     plotH = 6,
                     plotW = 8,
                     plotunits = "in",
                     dir_plots,
                     dir_sub = "_WoE",
                     boo_plot = TRUE,
                     incaseLabel,
                     targetSampleLabels = targetSampleLabels) {

  # Debugging
  boo.DEBUG <- FALSE
  #
  if (boo.DEBUG == TRUE) {##IF.boo.DEBUG.START
    TargetSiteID = TargetSiteID
    stressors.ssi = stressors.ssi
    df_stressinfo = df_stressorMetadata
    df_paired = df_PairedStressResp
    biocomm = bioComm
    df_bioMetricData = bioMetricData
    plotvars = plotvars
    plotdpi = 600
    plotH = 6
    plotW = 8
    plotunits = "in"
    dir_plots = dir_results
    dir_sub = "_WoE"
    boo_plot = boo_plot_user
  }##IF.boo.DEBUG.END

  if(is.null(config) == FALSE){
    plotvars         <- config$data_plotvars
    plotdpi          <- config$plot_dpi
    plotH            <- config$plot_H
    plotW            <- config$plot_W
    plotunits        <- config$plot_units
  }

  # define pipe; standardize biocomm
  biocomm_upper <- toupper(biocomm)

  # Initialize gap df
  df_gap <- data.frame(fxnname = character(), condition = character(), result = character(), comment = character())

  # Prepare generic plot variables ----
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

  ## Plot, arrow labels
  aLabPos <- "1"
  aLabZero <- "0"
  aLabNeg <- "-1"

  # Pull ssi
  ssi <- info.ssi |> dplyr::distinct(MetricName) |> dplyr::pull(MetricName)

  ## SSI loop
  df.box.scores <- data.frame()
  df.log.scores <- data.frame()

  for(i in seq_along(ssi)){

    # Current ssi and associated stressors
    current_ssi <- ssi[i]
    current_stressors <- info.stress.ssi |>
      dplyr::filter(SSIndex == current_ssi) |>
      dplyr::pull(Stressor)

    # Data gaps: check that ssi is in the biological metric data
    if (current_ssi %in% colnames(df_bioMetricData) == FALSE) {
      gap.statement <- data.frame(
        fxnname = "getVPSSI",
        condition = "No stressor-specific index data available",
        result = ssi.name,
        comment = paste0(ssi.name, " not in ", biocomm, " metric data.")
      )

      df_gap <- df_gap |> dplyr::bind_rows(gap.statement)
      next
    }

    # Add ssi to stressor data and add ssi metadata
    df.data <- df_paired |>
      dplyr::filter(IncaseYN == 1) |>
      dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                    RespSampleDate, IncaseYN, BioComm, Quality, dplyr::all_of(current_stressors))

    df.ssi.metrics <- df_bioMetricData |>
      dplyr:: select(StationID, RespSampleID, RespSampleDate, dplyr::all_of(current_ssi))

    df.data <- df.data |>
      dplyr::left_join(df.ssi.metrics, by = c("StationID", "RespSampleID", "RespSampleDate")) |>
      tidyr::pivot_longer(cols = dplyr::all_of(current_ssi), names_to = "SSIndex",
                          values_to = "SSIValue") |>
      tidyr::pivot_longer(cols = dplyr::all_of(current_stressors), names_to = "Stressor",
                          values_to = "StressorValue")|>
      dplyr::left_join(info.ssi, by = c("SSIndex" = "MetricName"))

    # Co-occurrence data
    ## Filter to not degraded
    df.comp <- df.data |>
      dplyr::filter(Quality == "Not degraded",
                    StationID != TargetSiteID)

    ## Calculate statistics
    df.stats <- df.comp |>
      dplyr::mutate(NotNA = ifelse(!is.na(SSIValue), 1, 0)) |>
      dplyr::group_by(SSIndex) |>
      dplyr::summarize(n = sum(NotNA, na.rm = TRUE),
                       q25 = stats::quantile(SSIValue, probs = 0.25, na.rm = TRUE),
                       q50 = stats::quantile(SSIValue, probs = 0.50, na.rm = TRUE),
                       q75 = stats::quantile(SSIValue, probs = 0.75, na.rm = TRUE),
                       minVal = min(SSIValue, na.rm = TRUE),
                       maxVal = max(SSIValue, na.rm = TRUE))

    current.box.scores <- df.data |>
      dplyr::filter(StationID == TargetSiteID) |>
      dplyr::filter(is.na(SSIValue) == FALSE) |>
      dplyr::left_join(df.stats, by = "SSIndex") |>
      dplyr::mutate(Sc_Box = dplyr::case_when(
        (TrendWIncStress == "Dec") ~ dplyr::case_when(SSIValue > q50 ~ -1,
                                                      SSIValue < q25 ~ 1,
                                                      TRUE ~ 0),
        (TrendWIncStress == "Inc") ~ dplyr::case_when(SSIValue > q75 ~ 1,
                                                      SSIValue < q50 ~ -1,
                                                      TRUE ~ 0)),
        bioComm = biocomm) |>
      dplyr::mutate(Score_text = dplyr::case_when(
        Sc_Box == 1 ~ "1 (Supporting)",
        Sc_Box == 0 ~ "0 (Indeterminate)",
        Sc_Box == -1 ~ "-1 (Refuting)"
      ))


    # Boxplot aesthetics
    if(current.box.scores$TrendWIncStress[1] == "Dec"){
      box_qHI <- current.box.scores$q50[1]
      box_qLO <- current.box.scores$q25[1]

      str_subtitle <- paste0("Is the target sample stressor-specific index depressed compared to those from unimpaired, comparator samples?") |>
        stringr::str_wrap(100)

    } else if(current.box.scores$TrendWIncStress[1] == "Inc"){
      box_qHI <- current.box.scores$q75[1]
      box_qLO <- current.box.scores$q50[1]

      str_subtitle <- paste0("Is the target sample stressor-specific index elevated compared to those from unimpaired, comparator samples?")|>
        stringr::str_wrap(100)
    }

    maxVal <- current.box.scores$maxVal[1]
    minVal <- current.box.scores$minVal[1]

    segNeg <- ((maxVal - box_qHI) / 2) + box_qHI
    segZero <- ((box_qHI - box_qLO) / 2) + box_qLO
    segPos <- ((box_qLO - minVal) / 2) + minVal

    str_title <- paste0(TargetSiteID, ": Stressor-specific index - Co-occurrence line of evidence ",
                        "for ", current.box.scores$MetricLabel[1]) |>  stringr::str_wrap(75)

    lab.N <- paste0("n = ", current.box.scores$n[1])
    lab.Score <- paste0("Score = ", paste0(current.box.scores$Score_text, collapse = ", "))

    str_caption <- paste0("Inside-the-case samples with paired stressor/response ",
                          "samples (", lab.N, ").\n", lab.Score, ".\n*Sample ",
                          "quality rated based on overall biological index.")

    str_ylab <- paste0("Inside-the-case samples selected from = ", incaseLabel)

    msg <- paste0("printing boxplot for ", current_ssi)
    message(msg)

    temp_df <- data.frame(xval = current.box.scores$SSIValue, yval = 0.9, labelval = current.box.scores$RespSampleID, Case = "Target sample(s)")

    # plot boxplot
    p1 <- ggplot2::ggplot(df.comp, ggplot2::aes(y = SSIValue,
                                                        x = 1,
                                                        group = 1)) +
      ggplot2::geom_boxplot(outliers = FALSE,
                            na.rm = TRUE, staplewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(data = temp_df, ggplot2::aes(yintercept = xval, lty = Case), color = targ_line_col,
                          lwd = targ_line_lwd, na.rm = TRUE, show.legend = FALSE) +
      # dummy line to get orientation correct for legend
      ggplot2::geom_vline(data = data.frame(xintercept = 1 - 0.5, lab = "Target sample(s)"), ggplot2::aes(xintercept = xintercept, linetype = lab), inherit.aes = FALSE, alpha = 0) +
      ggplot2::guides(linetype = ggplot2::guide_legend(order = 1, override.aes = list( alpha = 1, colour = targ_line_col, linewidth = targ_line_lwd )))+
      ggplot2::scale_linetype_manual(name = "", values = targ_line_lty)+
      ggplot2::geom_jitter(ggplot2::aes(color = "black", shape = Quality,
                                        fill = Quality), alpha = 0.5,
                           na.rm = TRUE, width = 0.25, height = 0) +
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
      ggplot2::annotate("text", x = 1.52, y = c(segNeg, segZero, segPos),
                        label = c(aLabNeg, aLabZero, aLabPos), color = "orange") +
      ggplot2::scale_color_manual(name = "Samples*",
                                  breaks = c("Not degraded", "Degraded"),
                                  values = bio_col[2:3], drop = TRUE) +
      ggplot2::scale_fill_manual(name = "Samples*",
                                 breaks = c("Not degraded", "Degraded"),
                                 values = bio_fill[2:3], drop = TRUE) +
      ggplot2::scale_shape_manual(name = "Samples*",
                                  breaks = c("Not degraded", "Degraded"),
                                  values = bio_shape[2:3], drop = TRUE) +
      ggplot2::labs(title = str_title,
                    subtitle = str_subtitle,
                    caption = str_caption,
                    y = current.box.scores$MetricLabel[1],
                    x = str_ylab) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
                     plot.caption = ggplot2::element_text(size = 8)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())

    if(targetSampleLabels == TRUE){
      p1 <- p1+
        ggrepel::geom_label_repel(data = temp_df, ggplot2::aes(x = 0.5, y = xval, label = labelval, group = NA), color = targ_line_col,
                                  size = 8/ggplot2::.pt)
    }

    for(k in seq_along(current_stressors)){
      stress <- current_stressors[k]
      slab <- info.stress.ssi |> dplyr::filter(Stressor == stress) |> dplyr::pull(Label)
      dir_path_stress <- file.path(dirname(dir_plots), basename(dir_plots), TargetSiteID, biocomm_upper, stress)

      # Sufficiency LoE
      if(is.na(current.box.scores$CutoffValue[1]) == FALSE){

        # Logistic regression scoring
        current_df.suff <- df.data |>
          dplyr::filter(Stressor == stress) |>
          dplyr::rowwise() |>
          dplyr::mutate(SSIQuality = qualityAssign(metric_ind = InclusiveIndicator,
                                                   metric_dir = TrendWIncStress,
                                                   metric_cutval = CutoffValue,
                                                   metric = SSIValue),
                        y.name = ifelse(SSIQuality == "Degraded", 1, 0)) |>
          dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID, RespSampleDate, Stressor, StressorValue, SSIndex, SSIValue, SSIQuality, y.name) |>
          tidyr::drop_na()

        fit <- stats::glm(y.name ~ StressorValue, data = current_df.suff, family = stats::binomial)

        current.log.scores <- current_df.suff |>
          dplyr::filter(StationID == TargetSiteID)

        newdat <- data.frame(StressorValue = seq(min(current_df.suff$StressorValue, na.rm = TRUE),
                                                 max(current_df.suff$StressorValue, na.rm = TRUE), len = 100))
        newdat$y.name <- stats::predict(fit, newdata = newdat, type = "response")

        current.log.scores$VPpred_Deg <- stats::predict(fit, newdata = data.frame(StressorValue = current.log.scores$StressorValue), type = "response")

        current.log.scores <- current.log.scores |>
          dplyr::mutate(Sc_VPlog = dplyr::case_when(
            VPpred_Deg <= 0.2 ~ -1,
            (VPpred_Deg > 0.2 & VPpred_Deg <= 0.5) ~ 0,
            VPpred_Deg > 0.5 ~ 1,
            TRUE ~ NA
          )) |>
          dplyr::mutate(Score_text = dplyr::case_when(
            Sc_VPlog == -1 ~ "-1 (Refuting)",
            Sc_VPlog == 0 ~ "0 (Indeterminate)",
            Sc_VPlog == 1 ~ "1 (Supporting)"
          )) |>
          dplyr::distinct() |>
          dplyr::mutate(n = nrow(current_df.suff)) |>
          dplyr::select(-y.name)

        # define filename for logistic regression plot
        fn_png_p2 <- paste0(TargetSiteID, stress, "_", biocomm_upper, "_", current_ssi,
                            "_VIreg.png")

        # Identify plot vars for logistic plots
        negStart <- 0
        negEnd   <- 0.2  # same as zeroStart
        zeroEnd  <- 0.5  # same as posStart
        posEnd   <- 1
        midNeg   <- ((negEnd - negStart) / 2) + negStart
        midZero  <- ((zeroEnd - negEnd) / 2) + negEnd
        midPos   <- ((posEnd - zeroEnd) / 2) + zeroEnd

        # plot logistic regressions ####
        ylabel <- paste0("Relative probability of degraded condition (", ssi, ")")
        maintitleSR <- paste0(TargetSiteID,
                              ": Stressor-specific index - Sufficiency line of evidence")
        subtitleSR <-"Is the target sample stressor-specific index value associated with a higher probability of biological impairment?" |>
          stringr::str_wrap(75)

        captionSR <- paste(paste0("All inside-the-case samples (n = ", unique(current.log.scores$n), ")."),
                           paste0("Score = ", paste(current.log.scores$Score_text, collapse = ", "), "."),
                           sep = "\n")

        xmin <- min(current_df.suff$StressorValue, na.rm = TRUE)
        xmax <- max(current_df.suff$StressorValue, na.rm = TRUE)
        xseg <- xmax + (0.02 * xmax)

        msg <- paste0("printing logistic regression for ", stress, " against ", current_ssi)
        message(msg)


        temp_df_suff <- data.frame(xval = current.log.scores$StressorValue, yval = 0.9, labelval = current.log.scores$RespSampleID, Case = "Target sample(s)")

        p2 <- ggplot2::ggplot(current_df.suff, ggplot2::aes(x = StressorValue, y = y.name)) +
          ggplot2::geom_point(ggplot2::aes(color = "black", shape = SSIQuality,
                                           fill = SSIQuality, size = SSIQuality),
                              alpha = 0.5, size = 2, na.rm = TRUE) +
          ggplot2::geom_line(ggplot2::aes(y = y.name, x = StressorValue), data = newdat
                             , color = "black", lwd = 1, na.rm = TRUE) +
          ggplot2::scale_fill_manual(name = "Samples",
                                     breaks = c("Not degraded", "Degraded"),
                                     values = bio_col[2:3], drop = FALSE) +
          ggplot2::scale_color_manual(name = "Samples",
                                      breaks = c("Not degraded", "Degraded"),
                                      values = bio_col[2:3], drop = FALSE) +
          ggplot2::scale_shape_manual(name = "Samples",
                                      breaks = c("Not degraded", "Degraded"),
                                      values = bio_shape[2:3], drop = FALSE) +
          ggplot2::scale_size_manual(name = "Samples",
                                     breaks = c("Not degraded", "Degraded"),
                                     values = bio_size[2:3], drop = FALSE) +
          ggplot2::geom_vline(data = temp_df_suff, ggplot2::aes(xintercept = xval, lty = Case), color = targ_line_col,
                              lwd = targ_line_lwd,
                              na.rm = TRUE) +
          ggplot2::scale_linetype_manual(name = "", values =  targ_line_lty)+
          ggplot2::geom_hline(yintercept = c(0.2, 0.5), color = "black"
                              , lty = 2, na.rm = TRUE) +
          ggplot2::annotate("segment", y = negStart, yend = negEnd, x = xseg,
                            color = "orange", linewidth = 0.7, alpha = 0.6,
                            arrow = grid::arrow(ends = "both", type = "open",
                                                length = grid::unit(0.2, "cm"))) +
          ggplot2::annotate("segment", y = negEnd, yend = zeroEnd, x = xseg,
                            color = "orange", linewidth = 0.7, alpha = 0.6,
                            arrow = grid::arrow(ends = "both", type = "open",
                                                length = grid::unit(0.2, "cm"))) +
          ggplot2::annotate("segment", y = zeroEnd, yend = posEnd, x = xseg,
                            color = "orange", linewidth = 0.7, alpha = 0.6,
                            arrow = grid::arrow(ends = "both", type = "open",
                                                length = grid::unit(0.2, "cm"))) +
          ggplot2::annotate("text", x = xmax, y = c(midNeg, midZero, midPos),
                            label = c(aLabNeg, aLabZero, aLabPos), color = "orange") +
          ggplot2::labs(y = ylabel, x = slab) +
          ggplot2::theme_bw() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)) +
          ggplot2::labs(title = maintitleSR, subtitle = subtitleSR,
                        caption = captionSR)

        if(targetSampleLabels == TRUE){
          p2 <- p2 +
            ggrepel::geom_label_repel(data = temp_df_suff, ggplot2::aes(x = xval, y = yval, label = labelval, group = NA), color = targ_line_col,
                                      size = 8/ggplot2::.pt)
        }

        # Save logistic regression
        ggplot2::ggsave(filename = file.path(dir_path_stress, fn_png_p2),
                        plot = p2, dpi = plotdpi, width = plotW,
                        height = plotH, units = plotunits)

        df.log.scores <- dplyr::bind_rows(df.log.scores, current.log.scores)
      } # sufficiency if

      fn_png_p1 <- paste0(TargetSiteID, "_", stress, "_", biocomm_upper, "_", current_ssi,
                          "_VIbox.png")

      ggplot2::ggsave(filename = file.path(dir_path_stress, fn_png_p1),
                      plot = p1, dpi = plotdpi, width = plotW, height = plotH,
                      units = plotunits)

    } # Seq along current stressors

    df.box.scores <- dplyr::bind_rows(df.box.scores, current.box.scores)
  } # Seq along ssi

  dir.path <- file.path(dirname(dir_plots), basename(dir_plots), TargetSiteID, biocomm_upper, dir_sub)

  df.box.scores_ret <- df.box.scores |>
    dplyr::left_join(info.stress.ssi |> dplyr::select(Stressor, Label), by = "Stressor") |>
    dplyr::select(-Stressor) |>
    dplyr::rename(Stressor = Label, bioIndexName = SSIndex, Label = MetricLabel,
                  bioIndex = SSIValue, Score = Sc_Box) |>
    dplyr::mutate(LoE = "VP_SSIbox", Quality = NA) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, bioComm, bioIndexName, bioIndex, Quality,
                  Stressor, StressorValue, LoE, Score)


  fn.scores.box <- paste(TargetSiteID, biocomm_upper, "VPSSIBox", "Scores.csv", sep = "_")
  utils::write.csv(df.box.scores, file.path(dir.path, fn.scores.box), row.names = FALSE)

  df.log.scores.ret <- data.frame()
  if (nrow(df.log.scores) > 0) {
    fn.scores.log <- paste(TargetSiteID, biocomm_upper, "VPSSILog", "Scores.csv", sep = "_")
    utils::write.csv(df.log.scores, file.path(dir.path, fn.scores.log), row.names = FALSE)

    df.log.scores.ret <- df.log.scores |>
      dplyr::left_join(info.stress.ssi |> dplyr::select(Stressor, Label), by = "Stressor") |>
      dplyr::select(-Stressor) |>
      dplyr::rename(Score = Sc_VPlog,
                    bioIndex = SSIValue, bioIndexName = SSIndex,
                    Quality = SSIQuality, Stressor = Label) |>
      dplyr::mutate(LoE = "VP_SSIlog", bioComm = biocomm) |>
      dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                    RespSampleDate, bioComm, bioIndexName, bioIndex, Quality,
                    Stressor, StressorValue, LoE, Score)
  }

  return(list(df.scores = dplyr::bind_rows(df.box.scores_ret, df.log.scores.ret),
              df_gap = df_gap))

}
