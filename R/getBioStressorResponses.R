#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Biological Stressor-Response Gradient Lines of Evidence
#'
#' @description Use linear regression to evaluate stressor-response gradients
#'              inside the case and outside the case.
#'
#' Requires packages dplyr, ggplot2, stringr, tidyr
#'
#' @param TargetSiteID target site identifier
#' @param df_stressinfo metadata for stressors that passed the co-occurrence screen
#' @param df_respinfo biological response metric metadata for the current biological community
#' @param df_respdata biological response metric data for the current biological community
#' @param df_datapaired paired stressor-response data for the current biological community
#' @param biocomm current biological community
#' @param bioindex name of index for current biological community
#' @param min_cases minimum number of paired samples to advance stressor through analysis
#' @param p.val_cutoff p-value cutoff above which the slope is not considered significant
#'                     Default = 0.05
#' @param r2_cutoff r2 value below which the relationship has too much variance
#'                  Default = 0.1
#' @param plotvars colors, shapes, fills, and transparencies for each type (target,
#'                 not degraded, degraded, inside-the-case, and outside-the-case)
#' @param refOutline color of the reference sites outline
#' @param plotdpi standardized dpi for all plots
#' @param plotH standardized height for all plots
#' @param plotW standardized width for all plots
#' @param plotunits standardized units for all plots
#' @param dir_sub subdirectory for weight of evidence outputs from this function. Default = "_WoE"
#' @param dir_plots directory to save plots
#' @param boo_plot Boolean value to save plots. Default = TRUE.
#'
#' @return Writes a stressor-response correlation table and table of scores.
#' Generates linear regression figures for inside and outside of the case
#' as well as a visual depiction of the correlation table. Returns a list containing
#' scores for the bioindex and datagaps table.
#'
#' @export
getBioStressorResponses <- function(TargetSiteID,
                                    df_stressinfo,
                                    df_respinfo,
                                    df_respdata,
                                    df_datapaired,
                                    config = NULL,
                                    biocomm,
                                    bioindex,
                                    min_cases = samplim,
                                    p.val_cutoff = 0.05,
                                    r2_cutoff = 0.1,
                                    plotvars = data_plotvars,
                                    refOutline = refOutline_col,
                                    plotdpi,
                                    plotH,
                                    plotW,
                                    plotunits,
                                    dir_plots,
                                    dir_sub = "_WoE",
                                    boo_plot = TRUE,
                                    targetSampleLabels = targetSampleLabels) {##FUNCTION.START

  boo.DEBUG <- FALSE

  if (boo.DEBUG == TRUE) {
    TargetSiteID = TargetSiteID
    df_stressinfo = df_stressorMetadata
    df_respinfo = bioMetricInfo
    df_respdata = bioMetricData
    df_datapaired = df_PairedStressResp
    biocomm = bioComm
    bioindex = bioIndex
    min_cases = samplim
    p.val_cutoff = 0.05
    r2_cutoff = 0.2
    plotvars = data_plotvars
    refOutline = refOutline_col
    plotdpi = plot_dpi
    plotH = plot_H
    plotW = plot_W
    plotunits = plot_units
    dir_plots = dir_plots
    dir_sub = "_WoE"
    boo_plot = TRUE
    targetSampleLabels = TRUE
  }

  if(is.null(config) == FALSE){
    plotvars      = config$data_plotvars
    refOutline    = config$refOutline_col
    plotdpi       = config$plot_dpi
    plotH         = config$plot_H
    plotW         = config$plot_W
    plotunits     = config$plot_units
    boo_plot      = config$boo_plot
  }

  # Definitions ####
  biocomm <- toupper(biocomm)

  # Write results directory ----
  out.folders <- c(dir_plots, TargetSiteID, biocomm, dir_sub)

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

  # Plot vars, inside-the-case
  plotvarsIn  <- plotvars |>
    dplyr::filter(Type %in% c("target", "insideND", "insideD"))
  bio_fill_inside    <- rev(unlist(plotvarsIn$Fill)) # Degraded, Not degraded, Target
  bio_shape_inside   <- rev(unlist(plotvarsIn$Shape)) # down triangle, circle, triangle
  bio_size_inside    <- rev(unlist(plotvarsIn$Size)) # Degraded, Not degraded, Target
  bio_alpha_inside   <- rev(unlist(plotvarsIn$Alpha)) # Degraded, Not degraded, Target

  # Plot vars, outside-the-case
  plotvarsOut  <- plotvars |>
    dplyr::filter(Type %in% c("target", "outsideND", "outsideD"))
  bio_fill_outside    <- rev(unlist(plotvarsOut$Fill)) # Degraded, Not degraded, Target
  bio_shape_outside   <- rev(unlist(plotvarsOut$Shape)) # down triangle, circle, triangle
  bio_size_outside    <- rev(unlist(plotvarsOut$Size)) # Degraded, Not degraded, Target
  bio_alpha_outside   <- rev(unlist(plotvarsOut$Alpha)) # Degraded, Not degraded, Target

  # Merge other metrics into paired dataset with transformed stressor values
  df_datapaired <- merge(df_datapaired, df_respdata,
                         by = c("StationID", "RespSampleID", "RespSampleDate",
                                bioindex, "Quality"),
                         all.x = TRUE)
  df_datapaired <- df_datapaired |> dplyr::select_if(not_all_na)

  df_SiteData <- df_datapaired[df_datapaired$StationID == TargetSiteID, ]
  df_CompData <- df_datapaired[df_datapaired$IncaseYN == 1, ]
  df_NotCompData <- df_datapaired[df_datapaired$OutcaseYN == 1 & df_datapaired$IncaseYN == 0, ]

  # Get stressors and responses for looping over
  stressors <- df_stressinfo$Stressor
  responses <- df_respinfo$MetricName

  p.len <- length(stressors)
  q.len <- length(responses)

  # Scores filepath
  fn_scores <- paste0(TargetSiteID, "_", biocomm, "_BG_Scores.csv")
  fp_scores <- file.path(dir_path, fn_scores)
  df.sc <- data.frame()

  # Initialize df.CorrTable
  fn_corr <- paste0(TargetSiteID, "_", biocomm, "_BG_InsideCorrs.csv")
  fp_corr <- file.path(dir_path, fn_corr)
  df.CorrTable <- data.frame()

  data_gaps <- data.frame()

  # FOR.p ####
  for (p in seq_along(stressors)) {

    stressName <- stressors[p]
    stressLabel <- df_stressinfo$Label[df_stressinfo$Stressor == stressName]

    # Write graphics directory ----
    out.dir <- dirname(dir_plots)
    out.folders <- c(out.dir, basename(dir_plots), TargetSiteID, biocomm, stressName)

    for (f in seq_along(out.folders)) {
      if (f == 1) {
        dir_path_stress <- file.path(out.folders[f])
      } else {
        dir_path_stress <- file.path(dir_path_stress, out.folders[f])
      }
      if (dir.exists(dir_path_stress) == FALSE) {
        dir.create(dir_path_stress)
      }
    }

    # FOR.q ####
    for (q in seq_along(responses)) {

      respName <- responses[q]
      respLabel <- as.character(df_respinfo$MetricLabel[df_respinfo$MetricName == respName])
      pq <- q.len * (p - 1) + q
      pq.len <- p.len * q.len

      # Determine expected direction of slope
      dirIncStress <- unique(df_stressinfo$DirIncStress[df_stressinfo$Stressor == stressName])
      dirIncStress <- tolower(dirIncStress)
      dirRespBad <- unique(df_respinfo$TrendWIncStress[df_respinfo$MetricName == respName])
      dirRespBad <- tolower(dirRespBad)
      if (dirIncStress == "inc" & dirRespBad == "dec") {
        # increasing stressor decreases biological integrity when bad values are
        # lower than good values; slope should be negative
        exp.dir <- -1
      } else if (dirIncStress == "inc" & dirRespBad == "inc") {
        # increasing stressor decreases biological integrity when bad values are
        # higher than good values; slope should be positive
        exp.dir <- 1
      } else if (dirIncStress == "dec" & dirRespBad == "dec") {
        # decreasing stressor decreases biological integrity when bad values are
        # lower than good values; slope should be positive
        exp.dir <- 1
      } else { # (dirIncStress == "dec" & dirRespBad == "inc")
        # decreasing stressor decreases biological integrity when bad values are
        # higher than good values; slope should be negative
        exp.dir <- -1
      }

      # Subset dataframes to current stressor and response
      df_SiteData_sr <- df_SiteData |>
        dplyr::select(RespSampleID, dplyr::all_of(respName), dplyr::all_of(stressName)) |>
        dplyr::rename(Response = respName, Stressor = stressName) |>
        dplyr::mutate(Case = "Target sample(s)") |>
        tidyr::drop_na()

      df_CompData_sr <- df_CompData |>
        dplyr::select(all_of(respName), all_of(stressName), Quality) |>
        dplyr::rename(Response = respName, Stressor = stressName) |>
        tidyr::drop_na()

      df_NotCompData_sr <- df_NotCompData |>
        dplyr::select(all_of(respName), all_of(stressName), Quality) |>
        dplyr::rename(Response = respName, Stressor = stressName) |>
        tidyr::drop_na()

      # Write sample number data gaps
      gaps_sr <- bioStrRspGaps(df_CompData_sr, df_NotCompData_sr, stressName, respName)
      data_gaps <- dplyr::bind_rows(data_gaps, gaps_sr$gap)

      # Inside figure and score
      if(gaps_sr$inside_ne == FALSE){

        # Scoring
        score_sr_inside <- bioStrRespScore(df_CompData_sr, "inside",
                                           p.val_cutoff, r2_cutoff, exp.dir)

        sr_inside_mod <- score_sr_inside$model_summary

        SRLin_Score_inside <- as.character(sr_inside_mod$score)

        # Plotting
        inside_plot <- bioStrRspPlot(score_sr_inside$data, df_SiteData_sr, "inside", stressLabel, respLabel,
                                     targetSampleLabels, TargetSiteID, sr_inside_mod$caption,
                                     stressName, respName, plotW, plotH, plotunits, plotdpi)

        # Correlation
        corr_gaps <- bioStrRespCorrGaps(df_CompData_sr, stressName, respName)

        if(nrow(corr_gaps)==0){
          df.corr.sr <- data.frame(StationID = TargetSiteID,
                                   biocomm = biocomm,
                                   stressName = stressName,
                                   stressLabel = stressLabel,
                                   respName = respName,
                                   respLabel = respLabel,
                                   n = nrow(df_CompData_sr),
                                   statistic = sr_inside_mod$tval,
                                   p.value = sr_inside_mod$pvalue,
                                   estimate = cor(score_sr_inside$data$Stressor, score_sr_inside$data$Response) |> round(2),
                                   r2 = sr_inside_mod$r2)

          df.CorrTable <- dplyr::bind_rows(df.CorrTable, df.corr.sr)
        } else{
          data_gaps <- dplyr::bind_rows(data_gaps, corr_gaps)
        }

      } else{
        SRLin_Score_inside <- "NE"
      }

      # Outside figure and score
      if(gaps_sr$outside_ne == FALSE){

        # Scoring
        score_sr_outside <- bioStrRespScore(df_NotCompData_sr, "outside",
                                            p.val_cutoff, r2_cutoff, exp.dir)

        SRLin_Score_outside <- as.character(score_sr_outside$model_summary$score)

        # Plotting
        outside_plot <- bioStrRspPlot(score_sr_outside$data, df_SiteData_sr, "outside", stressLabel, respLabel,
                                      targetSampleLabels, TargetSiteID, score_sr_outside$model_summary$caption,
                                      stressName, respName, plotW, plotH, plotunits, plotdpi)

      } else{
        SRLin_Score_outside <- "NE"
      }

      msg.status <- paste0("Item (", pq, "/", pq.len, "), ", stressName,
                           " (", p, "/", p.len, "), ", respName, " (", q,
                           "/", q.len, "); score = ", SRLin_Score_inside, " (inside), ", SRLin_Score_outside, " (outside)")
      message(msg.status)

      df.sc.sr <- data.frame(StationID = TargetSiteID,
                             biocomm = biocomm,
                             stressName = stressName,
                             stressLabel = stressLabel,
                             respName = respName,
                             respLabel = respLabel,
                             n_site = nrow(df_SiteData_sr),
                             n_comp = nrow(df_CompData_sr),
                             SRLin_Score_inside = SRLin_Score_inside,
                             n_out = nrow(df_NotCompData_sr),
                             SRLin_Score_outside = SRLin_Score_outside)

      df.sc <- dplyr::bind_rows(df.sc, df.sc.sr)
    } # response loop end
  } # stressor loop end

  bioStrCorrPlot(TargetSiteID, df.CorrTable, dir_path, plotW, plotH, plotunits, plotdpi)

  write.csv(df.CorrTable, fp_corr, row.names = FALSE)
  write.csv(df.sc, fp_scores, row.names = FALSE)

  df_SiteData_join <- df_SiteData |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, BioComm, Quality, dplyr::all_of(stressors), dplyr::all_of(bioindex)) |>
    tidyr::pivot_longer(cols = all_of(stressors), names_to = "Stressor", values_to = "StressorValue") |>
    tidyr::pivot_longer(cols = all_of(bioindex), names_to = "bioIndexName", values_to = "bioIndex")

  df.scores <- df.sc |>
    dplyr::filter(respName == bioindex) |>
    dplyr::select(StationID, stressName, stressLabel, SRLin_Score_inside, SRLin_Score_outside) |>
    dplyr::left_join(df_SiteData_join, by = c("StationID", "stressName"="Stressor")) |>
    dplyr::select(-stressName) |>
    dplyr::rename("Stressor" = "stressLabel", "bioComm" = "BioComm") |>
    tidyr::pivot_longer(cols = c(SRLin_Score_inside, SRLin_Score_outside), names_to = "LoE", values_to = "Score")

  return(list(df.scores = df.scores,
              df_gap = data_gaps))

} # function end

#' Generate scores for the biological gradient line of evidence
#'
#' @param data paired stressor response data filtered to current stressor, response, and case
#' @param type case (inside or outside)
#' @param p.val_cutoff p value cutoff
#' @param r2_cutoff r2 cutoff
#' @param exp.dir expected slope direction
#'
#' @return data with 75th percentile prediction interval and summary of regression model
#' @export
bioStrRespScore <- function(data, type = c("inside", "outside"),  p.val_cutoff, r2_cutoff, exp.dir){
  model <- stats::lm(Response ~ Stressor, data = data,
                     na.action = stats::na.exclude)

  model_pred <- stats::predict(model, data, interval = "prediction", level = 0.75)

  data <- cbind(data, model_pred)

  model_summary <- data.frame(intercept = model$coefficients[1] |> unname() |> round(1),
                              slope = model$coefficients[2] |> unname() |> round(1),
                              r2 = summary(model)$r.squared |> round(3),
                              pvalue = summary(model)$coefficients[2,4] |> round(3),
                              n = nrow(data),
                              tval = summary(model)$coefficients[2,3] |> round(3))

  if(model_summary$r2 > r2_cutoff & model_summary$pvalue < p.val_cutoff){

    slope.dir <- ifelse(model_summary$slope < 0, -1, 1)

    if(slope.dir == exp.dir){
      model_summary$score <- 1
      model_summary$score_txt <- paste0("1 (Supporting): R2 > ", r2_cutoff, ", p-value < ", p.val_cutoff)
    } else{
      model_summary$score <- -1
      model_summary$score_text <- paste0("-1 (Refuting): R2 > ", r2_cutoff, ", p-value < ", p.val_cutoff)
    }

  } else{
    model_summary$score <- 0
    model_summary$score_txt <- paste0("0 (Indeterminate): R2 < ", r2_cutoff, ", p-value > ", p.val_cutoff)
  }

  model_summary$caption <- paste0("Regression (", type, "-the-case samples): y = ",
                                  model_summary$slope, " x + ", model_summary$intercept,
                                  " ~ R2 = ", model_summary$r2, " ~ p-value = ",
                                  model_summary$pvalue, " ~ n = ", model_summary$n,
                                  "\nscore = ", model_summary$score_txt)

  return(list(data = data,
              model_summary = model_summary))
}

#' Generate biological gradient (linear regression) plot
#'
#' @param data paired stressor response data filtered to current stressor, response, and case
#' @param site_data target site paired stressor response data filtered to current stressor, response, and case
#' @param type case (inside or outside)
#' @param stressLabel current stressor label
#' @param respLabel current response label
#' @param targetSampleLabels boolean to label target site samples with response ID
#' @param TargetSiteID target site identifier
#' @param caption caption associated with scoring
#' @param stressName current stressor name
#' @param respName current response name
#' @param width plot width
#' @param height plot height
#' @param units plot units
#' @param dpi plot resolution
#'
#' @return Saves plots. No return.
#' @export
bioStrRspPlot <- function(data, site_data, type = c("inside", "outside"), stressLabel, respLabel, targetSampleLabels, TargetSiteID, caption, stressName, respName, width, height, units, dpi){
  plot1 <- ggplot2::ggplot() +
    ggplot2::geom_smooth(data = data,
                         ggplot2::aes(x = Stressor, y = Response),
                         method = stats::lm,
                         formula = 'y ~ x',
                         color = "black",
                         fill = "black",
                         alpha = 0.2,
                         show.legend = FALSE,
                         na.rm = TRUE) +
    ggplot2::geom_line(data = data,
                       ggplot2::aes(x = Stressor, y = lwr),
                       color = "black",
                       linetype = "dashed",
                       show.legend = FALSE,
                       na.rm = TRUE)+
    ggplot2::geom_line(data = data,
                       ggplot2::aes(x = Stressor, y = upr),
                       color = "black",
                       linetype = "dashed",
                       show.legend = FALSE,
                       na.rm = TRUE) +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(x = Stressor, y = Response,
                                     color = Quality, fill = Quality,
                                     shape = Quality), alpha = 0.5,
                        na.rm = TRUE) +
    ggplot2::scale_fill_manual(name = "Quality",
                               breaks = c("Degraded", "Not degraded"),
                               values = get(paste0("bio_fill_", type)), drop = FALSE) +
    ggplot2::scale_color_manual(name = "Quality",
                                breaks = c("Degraded", "Not degraded"),
                                values = get(paste0("bio_fill_", type)), drop = FALSE) +
    ggplot2::scale_shape_manual(name = "Quality",
                                breaks = c("Degraded", "Not degraded"),
                                values = get(paste0("bio_shape_", type)), drop = FALSE) +
    ggplot2::geom_point(data = site_data ,
                        ggplot2::aes(x = Stressor, y = Response, size = Case),
                        color = "black", shape = get(paste0("bio_shape_", type))[3],
                        fill = get(paste0("bio_fill_", type))[3],
                        na.rm = TRUE) +
    ggplot2::scale_size_manual(name = "", values = get(paste0("bio_size_", type))[3]*1.5)+
    ggplot2::labs(x = stressLabel,
                  y = respLabel,
                  title = paste0(TargetSiteID, ": Biological gradient (linear regression) line of evidence"),
                  subtitle = paste("Is there a linear relationship between stressor and biological response from", type,
                                   "the case?\nLinear regression with 75th percentile prediction interval"),
                  caption = caption)+
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
                   plot.caption = ggplot2::element_text(size = 8),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 10))

  if(targetSampleLabels == TRUE){
    plot1 <- plot1 +
      ggrepel::geom_label_repel(data = site_data,
                                ggplot2::aes(x = Stressor, y = Response, label = RespSampleID),
                                color = get(paste0("bio_fill_", type))[3],
                                size = 8/ggplot2::.pt)
  }

  file_append <- ifelse(type == "inside", "GI.png", "GO.png")

  fn_png_out <- paste(TargetSiteID, stressName, biocomm, respName, file_append, sep = "_")
  fn_png_out <- file.path(dir_path_stress, fn_png_out)

  if (boo_plot) {
    ggplot2::ggsave(fn_png_out, plot1, width = plotW, height = plotH,
                    units = plotunits, dpi = plotdpi)
  }

}

#' Write gaps for running biological gradient line of evidence
#'
#' @param df_CompData_sr inside the case paired stressor response data filtered to current stressor and response
#' @param df_NotCompData_sr outside the case paired stressor response data filtered to current stressor and response
#' @param stressName current stressor name
#' @param respName current response name
#'
#' @return
#' @export Returns list with gap statement dataframe and boolean indicators of
#' to evaluate biological gradient line of evidence for inside and outside the case.
bioStrRspGaps <- function(df_CompData_sr, df_NotCompData_sr, stressName, respName){
  n_inside <- nrow(df_CompData_sr)
  n_outside <- nrow(df_NotCompData_sr)

  inside_ne <- FALSE
  outside_ne <- FALSE

  inside.gap.statement <- data.frame()
  outside.gap.statement <- data.frame()

  inside.condition = paste0("Number of inside-the-case samples; stressor: ", stressName, "; response: ", respName)
  outside.condition = paste0("Number of outside-the-case samples; stressor: ", stressName, "; response: ", respName)

  # Inside the case gaps
  if(n_inside < 2){ # only possible if min_cases set to 1
    inside.gap.statement <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = inside.condition,
      result = as.character(n_inside),
      comment =  "< 2 samples")

    inside_ne <- TRUE
  } else if(n_inside < min_cases){ # should not happen due to earlier screening
    inside.gap.statement <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = inside.condition,
      result = as.character(n_inside),
      comment =  paste0("< ", min_cases, " samples"))

    inside_ne <- TRUE
  }

  # Outside the case gaps
  if(n_outside == 0){
    outside.gap.statement <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = outside.condition,
      result = as.character(n_outside),
      comment =  "0 samples")

    outside_ne <- TRUE
  } else if(n_outside < 2){
    outside.gap.statement <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = outside.condition,
      result = as.character(n_outside),
      comment =  "< 2 samples")

    outside_ne <- TRUE
  } else if(n_outside < min_cases){
    outside.gap.statement <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = outside.condition,
      result = as.character(n_outside),
      comment =  paste0("< ", min_cases, " samples"))

    outside_ne <- TRUE
  }

  gaps.list <- list(gap = dplyr::bind_rows(inside.gap.statement, outside.gap.statement),
                    inside_ne = inside_ne,
                    outside_ne = outside_ne)

  return(gaps.list)
}

#' Generate stressor response correlation plot
#'
#' @param TargetSiteID target site identifier
#' @param df_corr dataframe containing correlation coefficients for stressor response pairs
#' @param dir_path path to write plot
#' @param width plot width
#' @param height plot height
#' @param units plot units
#' @param dpi plot resoltuion
#'
#' @return Saves correlation plot. No return
#' @export
bioStrCorrPlot <- function(TargetSiteID, df_corr, dir_path, width, height, units, dpi){

  wrap_length <- 28

  p_cp <- ggplot2::ggplot(df_corr, ggplot2::aes(x = stringr::str_wrap(respLabel, wrap_length),
                                                y = stringr::str_wrap(stressLabel, wrap_length),
                                                fill = estimate)) +
    ggplot2::geom_tile(color = "black", lwd = 1, linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = estimate), color = "black", size = 2.5) +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0,
                                  limits = c(-1, 1), guide = "colorbar") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.3)),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(1)),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
                   axis.title = ggplot2::element_text(size = ggplot2::rel(1.3)),
                   axis.text.x = ggplot2::element_text(color = "black", size = ggplot2::rel(1), angle = 45,
                                                       hjust = 1),
                   axis.text.y = ggplot2::element_text(color = "black", size = ggplot2::rel(1), angle = 30)) +
    ggplot2::labs(title = paste0(TargetSiteID, ": Stressor-Response Correlations"),
                  x = "Responses",
                  y = "Stressors",
                  fill = "Estimate")

  # Save correlation plot
  fn_png_cp <- file.path(dir_path,
                         paste0(TargetSiteID, "_", biocomm, "_InsideCorrPlot.png"))
  ggplot2::ggsave(fn_png_cp, p_cp, width = plotH, height = plotW,
                  units = plotunits, dpi = plotdpi)

}

#' Generate correlation gaps
#'
#' @param data inside the case paired stressor-response data
#' @param stressName current stressor name
#' @param respName current response name
#'
#' @return Correlation gaps dataframe
#' @export
bioStrRespCorrGaps <- function(data, stressName, respName){
  gap.statement.str <- data.frame()
  gap.statement.resp <- data.frame()

  if (stats::sd(data$Stressor) == 0) { # Vertical line
    gap.statement.str <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = paste(stressName, respName, sep = ":"),
      result = "0",
      comment = "Stressor data in the inside-the-case sample set have a standard deviation of zero: all values are equal"
    )
  }
  if (stats::sd(data$Response) == 0) { # Horizontal line
    gap.statement.resp <- data.frame(
      fxnname = "getBioStressorResponse",
      condition = paste(stressName, respName, sep = ":"),
      result = "0",
      comment = "Response data in the inside-the-case sample set have a standard deviation of zero: all values are equal"
    )
  }
  return(dplyr::bind_rows(gap.statement.str, gap.statement.resp))
}
