#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Stressor Sufficiency Line of Evidence
#'
#' @description Generates logistic regressions for the sufficiency line of evidence
#'
#' @param TargetSiteID target site identifier
#' @param df_data paired stressor-response data for the current biological community
#' @param df_stressinfo  metadata for stressors that passed the co-occurrence screen
#' @param df_resp biological response metric data for the current biological community
#' @param df_respinfo biological response metadata for the current biological community
#' @param biocomm current biological community
#' @param colBio df_data column name for the field with biological index value.
#' @param config list containing configuration parameters
#' @param plotvars standardized sizes, fills, shapes, and transparencies for plots.
#' @param plotdpi standardized plot dpi. Defaults to plot_dpi. 600
#' @param plotH standardized plot height. Defaults to plot_H. 6
#' @param plotW standardized plot width. Defaults to plot_W. 8
#' @param plotunits standardized plot units. Defaults to plot_units. "in"
#' @param dir_plots directory to save plots.
#' @param dir_sub subdirectory for weight of evidence outputs from this function. Default = "_WoE"
#' @param boo_plot boolean value indicating whether or not to generate the plot. Defaults to TRUE.
#' @param targetSampleLabels boolean to label target sample values
#'
#' @return Writes individual plots as pngs and all sufficiency scores to a csv. Returns sufficiency scores for the biological index.
#' @export
getSufficiency <- function(TargetSiteID,
                           df_data,
                           df_stressinfo,
                           df_resp,
                           df_respinfo,
                           biocomm,
                           colBio,
                           config = NULL,
                           plotvars,
                           plotdpi = 600,
                           plotH = 6,
                           plotW = 8,
                           plotunits = "in",
                           dir_plots,
                           dir_sub = "_WoE",
                           boo_plot = TRUE,
                           targetSampleLabels = TRUE) {##FUNCTION.START

  # Debug
  boo_DEBUG <- FALSE

  if (boo_DEBUG==TRUE) {
    TargetSiteID = TargetSiteID
    df_data = df_PairedStressResp
    compSites = list.CompSites$comp.sites
    df_stressinfo = df_stressorMetadata
    biocomm = bioComm
    colBio = bioIndex
    plotvars = plotvars
    plotdpi = 600
    plotH = 6
    plotW = 8
    plotunits = "in"
    dir_plots = dir_results
    dir_sub = "_WoE"
    boo_plot = boo_plot_user
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

  # Get stressors
  stressors <- df_stressinfo$Stressor

  # Get responses with cutoff values
  responses <- df_respinfo |>
    dplyr::filter(is.na(CutoffValue) == FALSE) |>
    dplyr::pull(MetricName)

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

  # define scoring limits
  negStart <- 0
  negEnd   <- 0.2  # same as zeroStart
  zeroEnd  <- 0.5  # same as posStart
  posEnd   <- 1
  midNeg   <- ((negEnd - negStart) / 2) + negStart
  midZero  <- ((zeroEnd - negEnd) / 2) + negEnd
  midPos   <- ((posEnd - zeroEnd) / 2) + zeroEnd
  # arrow labels
  aLabNeg  <- "-1"
  aLabZero <- "0"
  aLabPos  <- "1"

  ## Plot, Variables
  ## Prepare colors, sizes, etc  ----
  plotvarsIn  <- plotvars |>
    dplyr::filter(Type %in% c("insideND", "insideD"))
  bio_col     <- rev(unlist(plotvarsIn$Fill)) # Degraded, Not degraded
  bio_fill    <- rev(unlist(plotvarsIn$Fill)) # Degraded, Not degraded
  bio_shape   <- rev(unlist(plotvarsIn$Shape)) # down triangle, circle
  bio_size    <- rev(unlist(plotvarsIn$Size)) # Degraded, Not degraded
  bio_alpha   <- rev(unlist(plotvarsIn$Alpha)) # Degraded, Not degraded
  targ_line_col <- plotvars$Fill[plotvars$Type == "target"]
  targ_line_lty <- 2
  targ_line_lwd <- 0.5

  # Get dataset
  df_resp <- df_resp |>  # will need to edit this if df_resp isn't filtered by user selection
    dplyr::select("StationID", "RespSampleID", "RespSampleDate", dplyr::all_of(responses)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(responses), names_to = "Response", values_to = "ResponseValue")

  df_data <- df_data |>
    dplyr::filter(IncaseYN == 1) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, Quality, dplyr::all_of(stressors)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(stressors), names_to = "Stressor", values_to = "StressorValue") |>
    dplyr::left_join(df_resp, by = c("StationID", "RespSampleID", "RespSampleDate"))


  # Calculate quantiles on Comparator Sites
  # Loop, j ####
  tot.plot <- length(stressors) * length(responses)
  plot.count <- 1

  suff_score <- NULL

  for (j in seq_along(stressors)) { ##FOR.j.START
    str <- stressors[j]
    str_label <- df_stressinfo |>
      dplyr::filter(Stressor == str) |>
      dplyr::pull(Label)

    # Write graphics directory ----
    out.dir <- dirname(dir_plots)
    out.folders <- c(out.dir, basename(dir_plots), TargetSiteID, biocomm, str)

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

    for(r in seq_along(responses)){
      resp <- responses[r]
      resp_label <- df_respinfo |>
        dplyr::filter(MetricName == resp) |>
        dplyr::pull(MetricLabel)

      message(paste0("Processing item (", plot.count, "/", tot.plot, "); ", str, ", ", resp, "\n"))
      utils::flush.console()

      plot.count <- plot.count + 1

      df_data_sr <- df_data |>
        dplyr::filter(Stressor == str, Response == resp) |>
        tidyr::drop_na()

      # If the response is not the index, assign quality
      if(resp != colBio){
        resp_info <- df_respinfo |>
          dplyr::filter(MetricName == resp)
        resp_trend <- resp_info$TrendWIncStress
        resp_cutoff <- resp_info$CutoffValue
        resp_inclind <- resp_info$InclusiveIndicator

        df_data_sr <- df_data_sr |>
          dplyr::select(-Quality) |>
          dplyr::mutate(Quality = qualityAssign(metric_ind = resp_inclind,
                                                metric_dir = resp_trend,
                                                metric_cutval = resp_cutoff,
                                                metric = ResponseValue))
      }

      df_target_sr <- df_data_sr |>
        dplyr::filter(StationID == TargetSiteID) |>
        dplyr::filter(Stressor == str, Response == resp)

      df.plot <- df_data_sr |>
        dplyr::rename(y = ResponseValue, x = StressorValue) |>
        dplyr::mutate(y.name = ifelse(Quality == "Degraded", 1, 0))


      if(nrow(df.plot) > 0){
        fit <- stats::glm(y.name ~ x, data = df.plot, family = stats::binomial)

        j_values <- data.frame(x = df_target_sr$StressorValue)

        #  Stressor Response Curve
        n_cc_df_plot <- nrow(df.plot)

        # create data for curve (type "response" gives probabilities)
        newdat <- data.frame(x = seq(min(df.plot$x, na.rm = TRUE),
                                     max(df.plot$x, na.rm = TRUE), len = 100))
        newdat$y.name <- stats::predict(fit, newdata = newdat, type = "response")

        # Scoring
        j_SR_predict <- stats::predict(fit, newdata = j_values, type = "response")
        j_SR_score <- cut(j_SR_predict,
                          breaks = c(0, 0.2, 0.5, 1),
                          labels = c(-1, 0, 1))
        j_SR_score_txt <- ifelse(j_SR_score == "1", "1 (Supporting)",
                                 ifelse(j_SR_score == "0", "0 (Indeterminate)",
                                        ifelse(j_SR_score == "-1", "-1 (Refuting)", "NA")))

        # Plot
        fn_png_p1 <- paste(TargetSiteID, str, resp, biocomm, "SU.png", sep = "_")

        captionSR <- paste(paste0("All inside-the-case samples (n=", n_cc_df_plot, ")."),
                           paste0("Score = ", paste(j_SR_score_txt, collapse = ", "), "."),
                           sep = "\n")

        xmax <- max(df.plot$x, na.rm = TRUE)
        xseg <- xmax + (0.02 * xmax)

        temp_df <- data.frame(xval = df_target_sr$StressorValue, yval = 0.9, labelval = df_target_sr$RespSampleID, case = "Target sample(s)")

        # plot1, ggplot ####
        p1 <- ggplot2::ggplot(df.plot, ggplot2::aes(x = x, y = y.name)) +
          ggplot2::geom_point(ggplot2::aes(color = "black", shape = Quality,
                                           fill = Quality),
                              alpha = 0.5, size = 2, na.rm = TRUE) +
          ggplot2::scale_fill_manual(name = "Samples",
                                     breaks = c("Degraded", "Not degraded"),
                                     values = bio_col, drop = FALSE) +
          ggplot2::scale_color_manual(name = "Samples",
                                      breaks = c("Degraded", "Not degraded"),
                                      values = bio_col, drop = FALSE) +
          ggplot2::scale_shape_manual(name = "Samples",
                                      breaks = c("Degraded", "Not degraded"),
                                      values = bio_shape, drop = FALSE) +
          ggplot2::geom_vline(data = temp_df, ggplot2::aes(xintercept = xval, lty = case), color = targ_line_col,
                              lwd = targ_line_lwd,
                              na.rm = TRUE) +
          ggplot2::scale_linetype_manual(name = "", values = targ_line_lty)+
          ggplot2::geom_hline(yintercept = c(0.2, 0.5), color = "black",
                              lty = 2, na.rm = TRUE) +
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
          ggplot2::labs(y = "Relative probability of degraded condition", x = str_label) +
          ggplot2::geom_line(ggplot2::aes(y = y.name, x = x), data = newdat,
                             color = "black", lwd = 1, na.rm = TRUE) +
          ggplot2::theme_bw() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)) +
          ggplot2::labs(title = paste0(TargetSiteID, ": Sufficiency (logistic regression) line of evidence"),
                        subtitle = "Is the target sample stressor value associated with a higher probability of biological impairment?",
                        caption = captionSR)

        if(targetSampleLabels == TRUE){
          p1 <- p1 +
            ggrepel::geom_label_repel(data = temp_df, ggplot2::aes(x = xval, y = yval, label = labelval, group = NA), color = targ_line_col,
                                      size = 8/ggplot2::.pt)
        }

        if ((boo_plot) == TRUE) {
          ggplot2::ggsave(filename = file.path(dir_path_stress, fn_png_p1),
                          plot = p1, dpi = plotdpi, width = plotW, height = plotH,
                          units = plotunits)
        }

      } # if nrow(df.plot) > 0

      sr_score <- df_target_sr |>
        dplyr::mutate(n = n_cc_df_plot,
               StressLabel = str_label,
               ResponseLabel = resp_label,
               SRpred_Deg = j_SR_predict |> round(3),
               Sc_SRlog = j_SR_score)

      suff_score <- dplyr::bind_rows(suff_score, sr_score)

    } # responses end

  } # stressor end

  # Save scores file
  fn.scores <- file.path(dir_path, paste0(TargetSiteID, "_", biocomm, "_Suff_Scores.csv"))
  utils::write.csv(suff_score, file = fn.scores, row.names = FALSE)

  df.scores <- suff_score |>
    dplyr::filter(Response == colBio) |>
    dplyr::select(-Stressor) |>
    dplyr::rename(Score = Sc_SRlog,
                  Stressor = StressLabel) |>
    dplyr::mutate(LoE = "Suff",
                  bioComm = biocomm,
                  bioIndexName = Response,
                  bioIndex = ResponseValue,) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, bioComm, bioIndexName, bioIndex, Quality,
                  Stressor, StressorValue, LoE, Score)

  return(df.scores)

}##FUNCTION.END

#' Assign quality to biological metric value
#'
#' @param metric_ind inclusive indicator of quality cutoff for a biological metric (<, >, <=, >=)
#' @param metric_dir direction of stressor gradient (Inc, Dec)
#' @param metric_cutval cutoff value for biological metric quality
#' @param metric metric
#'
#' @return quality label associated with a biological metric value
#' @export
qualityAssign <- function(metric_ind, metric_dir, metric_cutval, metric){
  if (metric_ind == "\u2264" | metric_ind == "<=") {
    if (metric_dir == "Dec") {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = TRUE,
                                     labels = c("Degraded", "Not degraded"))
    } else {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = TRUE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else if (metric_ind == "\u2265" | metric_ind == ">=") {
    if (metric_dir == "Dec") {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = FALSE,
                                     labels = c("Degraded", "Not degraded"))
    } else {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = FALSE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else if (metric_ind == "<") {
    if (metric_dir == "Dec") {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = FALSE,
                                     labels = c("Degraded", "Not degraded"))
    } else {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = FALSE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else {
    if (metric_dir == "Dec") {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = TRUE,
                                     labels = c("Degraded", "Not degraded"))
    } else {
      quality <- cut(metric,
                                     breaks = c(-Inf, metric_cutval, Inf),
                                     right = TRUE,
                                     labels = c("Not degraded", "Degraded"))
    }
  }
  return(quality |> as.character())
}
