#  Copyright 2025 TetraTech. All rights reserved.
#  Use, copying, modification, or distribution of this file or any of its contents
#  is expressly prohibited without prior written permission of TetraTech.
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  R v4.4.3
#
#' @title Verified Predictions
#'
#' @description Run the stressor-specific tolerance values line of evidence
#'
#' @param TargetSiteID target site identifier
#' @param stressors.sstv vector of stressors with stressor-specific tolerance values
#' @param df_stressinfo dataframe of stressor metadata
#' @param df_paired dataframe of paired stressor-response samples
#' @param biocomm current biological community
#' @param df_bioTaxaData dataframe of relative abundance
#' @param df_MasterTaxa dataframe of master taxa with SSTV values determined for individual taxa
#' @param config list of configuation parameters
#' @param colBio name of index for the current biological community
#' @param plotvars Standardized sizes, fills, shapes, and transparencies for plots.
#' @param plotdpi Standardized plot dpi. Defaults to plot_dpi. 600
#' @param plotH Standardized plot height. Defaults to plot_H. 6
#' @param plotW Standardized plot width. Defaults to plot_W. 8
#' @param plotunits Standardized plot units. Defaults to plot_units. "in"
#' @param dir_plots directory to save plots
#' @param dir_sub name of directory to save scores. Default = "_WoE"
#' @param boo_plot boolean to save plots
#'
#' @return Writes stressor specific tolerance value scores for target site samples
#' and saves associated boxplots.

#' @export
#'
getVerifiedPredictions <- function(TargetSiteID,
                                   stressors.sstv,
                                   df_stressinfo,
                                   df_paired,
                                   biocomm,
                                   df_bioTaxaData,
                                   df_MasterTaxa,
                                   config = NULL,
                                   colBio,
                                   plotvars,
                                   plotdpi = 600,
                                   plotH = 6,
                                   plotW = 8,
                                   plotunits = "in",
                                   dir_plots,
                                   dir_sub = "_WoE",
                                   boo_plot = TRUE,
                                   targetSampleLabels = targetSampleLabels) {##FUNCTION.START

  # # Debugging
  boo.DEBUG <- FALSE
  # #
  if (boo.DEBUG == TRUE) {##IF.boo.DEBUG.START
    TargetSiteID = TargetSiteID
    stressors.sstv = stressors.sstv
    df_stressinfo = df_stressorMetadata
    df_paired = df_PairedStressResp
    biocomm = bioComm
    df_bioTaxaData = bioTaxaData
    df_MasterTaxa = bioMasterTaxa
    colBio = bioIndexGp
    plotvars = plotvars
    plotdpi = 600
    plotH = 6
    plotW = 8
    plotunits = "in"
    dir_plots = dir_results
    dir_sub = "_WoE"
    boo_plot = boo.plot.user
  }##IF.boo.DEBUG.END

  if(is.null(config) == FALSE){
    plotvars       <- config$data_plotvars
    plotdpi        <- config$plot_dpi
    plotH          <- config$plot_H
    plotW          <- config$plot_W
    plotunits      <- config$plot_units
    boo_plot       <- config$boo_plot
  }

  biocomm <- toupper(biocomm)
  biocomm_lwr <- tolower(biocomm)

  # Write results directory ----
  out.dir <- dirname(dir_plots)
  dir.path <- file.path(out.dir, basename(dir_plots), TargetSiteID, biocomm, dir_sub)
  plot.path <- file.path(out.dir, basename(dir_plots), TargetSiteID, biocomm)


  # SSTV metadata
  df_sstvinfo_abb <- df_stressinfo |>
    dplyr::filter(Stressor %in% stressors.sstv) |>
    dplyr::rename("SSTVname" = paste0("SSTVname.", biocomm_lwr),
                  "SensMax" = paste0("SensMax.", biocomm_lwr),
                  "SensMin" = paste0("SensMin.", biocomm_lwr))

  # Check that SSTVs are in MT file
  df_gap <- data.frame(fxnname = character(), condition = character(), result = character(), comment = character())
  sstv_not_MT <- setdiff(df_sstvinfo_abb$SSTVname, names(df_MasterTaxa))
  sstv_MT <- intersect(df_sstvinfo_abb$SSTVname, names(df_MasterTaxa))

  if(length(sstv_not_MT)>0){
    for(n in seq_along(sstv_not_MT)){
      gap.statement <- data.frame(
        fxnname = "getVerifiedPredictions",
        condition = sstv_not_MT[n],
        result = "0",
        comment = paste0(sstv_not_MT[n], " not available in ", biocomm, "master taxa file.")
      )

      df_gap <- df_gap |>
        dplyr::bind_rows(gap.statement)
    }

    # If no more remaining sstvs in MT file return NULL scores
    if(length(sstv_MT)==0){
      return(list(df.scores = NULL,
                  df_gap = df_gap))
    }

  }

  # Create SSTV sensitivity look up table
  df_sstvinfo <- df_sstvinfo_abb |>
    tidyr::pivot_longer(cols = c("SensMax", "SensMin")) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      result = list(assignSens(sens_input = value, sstv = SSTVname, df_MasterTaxa = df_MasterTaxa))
    ) |>
    dplyr::ungroup()|>
    dplyr::select(SSTVname, name, result) |>
    tidyr::pivot_wider(names_from = name, values_from = result)

  sstvs <- unique(df_sstvinfo$SSTVname)

  # Give labels to master taxa
  df_mt_sstv <- df_MasterTaxa |>
    dplyr::select(TaxonID, dplyr::all_of(sstvs)) |>
    dplyr::mutate(dplyr::across(everything(), as.character)) |>
    tidyr::pivot_longer(cols = !TaxonID, names_to = "SSTVname")|>
    dplyr::left_join(df_sstvinfo, by = "SSTVname") |>
    dplyr::mutate(Sensitivity = dplyr::case_when(
      value %in% (unlist(SensMin)) ~ "Most sensitive",
      value %in% (unlist(SensMax)) ~ "All sensitive",
      TRUE ~ NA
    )) |>
    dplyr::select(TaxonID, SSTVname, Sensitivity)

  # Duplicate most sensitive rows to include in all sensitive calculation
  df_mt_sstv_allsens <- df_mt_sstv |>
    dplyr::filter(Sensitivity == "Most sensitive") |>
    dplyr::mutate(Sensitivity = "All sensitive")

  df_mt_sstv <- df_mt_sstv |>
    dplyr::bind_rows(df_mt_sstv_allsens)

  # Generate target and comparators datasets
  df_bioTaxaData_sum <- df_bioTaxaData |>
    dplyr::full_join(df_mt_sstv, by = "TaxonID", relationship = "many-to-many")|>
    dplyr::group_by(StationID, RespSampleID, SSTVname, Sensitivity) |>
    dplyr::summarize(NumInds = sum(NumInd, na.rm = TRUE),
                     PctInds = sum(PctInd, na.rm = TRUE),
                     NumTaxa = dplyr::n(),
                     PctTaxa = sum(PctTaxa, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(is.na(SSTVname)==FALSE, is.na(Sensitivity) == FALSE) |>
    tidyr::pivot_longer(cols = c(NumInds, PctInds, NumTaxa, PctTaxa), names_to = "variable", values_to = "value")

  comparators <- df_paired |>
    dplyr::filter(IncaseYN == 1,
                  Quality == "Not degraded") |>
    dplyr::pull(RespSampleID)

  targets <- df_paired |>
    dplyr::filter(StationID == TargetSiteID) |>
    dplyr::pull(RespSampleID)

  df_bioTaxaData_comp <- df_bioTaxaData_sum |>
    dplyr::filter(RespSampleID %in% comparators)

  df_bioTaxaData_target <- df_bioTaxaData_sum |>
    dplyr::filter(RespSampleID %in% targets)

  # Score each metric
  df_comp_quant <- df_bioTaxaData_comp |>
    dplyr::group_by(SSTVname, Sensitivity, variable) |>
    dplyr::summarise(q25 = stats::quantile(value, probs = 0.25, na.rm = TRUE),
                     q50 = stats::quantile(value, probs = 0.50, na.rm = TRUE),
                     max = max(value, na.rm = TRUE))

  df_target_score <- df_bioTaxaData_target |>
    dplyr::left_join(df_comp_quant, by = c("SSTVname", "Sensitivity", "variable")) |>
    dplyr::mutate(Score = dplyr::case_when(value < q25 ~ 1,
                                           dplyr::between(value, q25, q50) ~ 0,
                                           value > q50 ~ -1))

  # Write summary score
  df_target_sample <- df_paired |>
    dplyr::filter(StationID == TargetSiteID) |>
    dplyr::mutate(bioIndexName = colBio) |>
    dplyr::rename(bioIndex = colBio) |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID,
                  RespSampleDate, BioComm, bioIndexName, bioIndex, Quality,
                  dplyr::all_of(stressors.sstv)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(stressors.sstv), names_to = "Stressor", values_to = "StressorValue") |>
    dplyr::left_join(df_stressinfo |> dplyr::select(Stressor, Label, "SSTVname" = all_of(paste0("SSTVname.", biocomm_lwr))), by = "Stressor")

  df_target_score_sum <- df_target_score |>
    dplyr::group_by(RespSampleID, SSTVname) |>
    dplyr::summarize(Score = mean(Score)) |>
    dplyr::left_join(df_target_sample, by = c("RespSampleID", "SSTVname")) |>
    dplyr::mutate(LoE = "VP_SSTV") |>
    dplyr::select(StationID, StressSampleID, StressSampleDate, RespSampleID, RespSampleDate,
                  BioComm, bioIndexName, bioIndex, Quality, Label, StressorValue, LoE, Score) |>
    dplyr::rename(bioComm = BioComm, Stressor = Label)

  fn_scores <-  file.path(dir.path, paste0(TargetSiteID, "_", biocomm,
                                           "_VP_SSTV_Scores.csv"))

  write.csv(df_target_score |> dplyr::select(-max), fn_scores, row.names = FALSE)

  # Plot
  for(i in seq_along(stressors.sstv)){
    stressorLabel <- df_sstvinfo_abb |> dplyr::filter(Stressor == stressors.sstv[i]) |> dplyr::pull(Label)
    sstv <- df_sstvinfo_abb |> dplyr::filter(Stressor == stressors.sstv[i]) |> dplyr::pull(SSTVname)

    sstvPlot(plotvars, TargetSiteID, stressorLabel, sstv, stressors.sstv[i],
             biocomm, colBio, df_bioTaxaData_comp, df_target_score, boo_plot,
             plot.path, plotdpi, plotW, plotH, plotunits)
  }

  return(list(df.scores = df_target_score_sum,
              df_gap = df_gap))
}


#' Assign sensitivity class
#'
#' @param sens_input sensitivity value specification
#' @param sstv applicable stressor specific tolerance value
#' @param df_MasterTaxa taxa metadata
#'
#' @return Returns SSTV values in the master taxa file belonging to a sensitivity class
#' @export
assignSens <- function(sens_input, sstv, df_MasterTaxa){
  sens_output <- unlist(stringr::str_split(sens_input, ",")) |> stringr::str_trim()

  # convert to numeric if possible
  numeric_indicator <- suppressWarnings(sens_output |> as.numeric() |> is.na() |> sum())
  if (numeric_indicator == 0) {
    sens_output <- as.numeric(sens_output)

    if (length(sens_output) > 1) {
      # Handle case where defines endpoints of a range
      sstv.gp <- df_MasterTaxa |>
        dplyr::select(dplyr::all_of(sstv)) |>
        dplyr::filter(dplyr::between(get(sstv), sens_output[1],
                                     sens_output[2])) |>
        dplyr::distinct() |>
        dplyr::pull() |>
        as.character()
    } else { # only one value
      sstv.gp <- sens_output |>
        as.character()
    }
  } else{
    sstv.gp <- sens_output
  }
  return(sstv.gp)
}

#' Stressor Specific Tolerance Value Plot
#'
#' @param plotvars plotting parameters
#' @param TargetSiteID target site identifier
#' @param stressorLabel stressor label associated with the current sstv
#' @param sstv current stressor specific tolerance value
#' @param biocomm current biological community
#' @param colBio name of the index for the current biological community
#' @param data_comp dataframe of sensitivity metric values or comparator samples
#' @param df_target_score dataframe with target sample SSTV scores
#' @param boo_plot boolean to generate figures
#' @param plot.path path to save plots
#'
#' @return No return. Saves SSTV plot to results folder
#' @export
sstvPlot <- function(plotvars, TargetSiteID, stressorLabel, sstv, stressor,
                     biocomm, colBio, data_comp, data_targ, boo_plot, plot.path,
                     plotdpi, plotW, plotH, plotunits){
  ## Plot, Variables
  ## Prepare colors, sizes, etc  ----
  plotvars  <- plotvars |>
    dplyr::filter(Type %in% c("target", "insideND", "insideD"))
  bio_fill    <- rev(unlist(plotvars$Fill)) # Degraded, Not degraded, Target
  bio_shape   <- rev(unlist(plotvars$Shape)) # down triangle, circle, and triangle
  bio_size    <- rev(unlist(plotvars$Size)) # Degraded, Not degraded, Target
  bio_alpha   <- rev(unlist(plotvars$Alpha)) # Degraded, Not degraded, Target

  # Prepare labels
  str_title <- paste0(TargetSiteID, ": Stressor-specific tolerance values line of evidence for ", stressorLabel) |> stringr::str_wrap(80)
  str_subtitle <- "Does the target sample have a lower number of individuals, percent of individuals, number of taxa, and percent of all sensitive and the most sensitive taxa compared to unimpaired, comparator samples?" |> stringr::str_wrap(100)
  legendtitle <- "Samples"
  str_xlab  <- ""

  data_targ <- data_targ |>
    dplyr::filter(SSTVname == sstv) |>
    dplyr::mutate(Quality = "Target")

  data_comp <- data_comp |>
    dplyr::filter(SSTVname == sstv) |>
    dplyr::mutate(Quality = "Not degraded")

  str_scores <-  data_targ |>
    dplyr::arrange(RespSampleID) |>
    dplyr::group_by(Sensitivity, variable, max) |>
    dplyr::summarise(Scores = toString(Score) |> stringr::str_replace_all("NA", "NE"),
                     .groups = "drop")

  str_scores <- str_scores |>
    dplyr::group_by(variable) |>
    dplyr::mutate(OverallMax = max(max),
                  Scores = paste0("Scores = ", Scores))

  fn_png_p1 <- paste0(TargetSiteID, "_", make.names(stressor), "_",
                      biocomm, "_", colBio, "_VPSSTV.png")

  p_tv <- ggplot2::ggplot(data = data_comp,
                          ggplot2::aes(x = Sensitivity, y = value, group = Sensitivity)) +
    ggplot2::geom_boxplot(data = data_comp,
                          outliers = FALSE,
                          na.rm = TRUE,
                          staplewidth = 0.5, linewidth = 0.25, width = 0.4) +
    ggplot2::geom_jitter(data = data_comp,
                         ggplot2::aes(color = Quality, shape = Quality,
                                      fill = Quality, alpha = Quality),
                         na.rm = TRUE, width = 0.15, height = 0.01,
                         size = 0.4) +
    ggplot2::geom_jitter(data = data_targ,
                         ggplot2::aes(color = Quality, shape = Quality,
                                      fill = Quality, alpha = Quality),
                         na.rm = TRUE, width = 0.2, height = 0.01,
                         size = 1.5) +
    ggplot2::geom_boxplot(data = data_comp,
                          ggplot2::aes(group = Sensitivity),
                          outliers = FALSE,
                          na.rm = TRUE,
                          staplewidth = 0.5, linewidth = 0.25, fill = NA, width = 0.4) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(. ~ variable, scales = "free") +
    ggplot2::scale_color_manual(name = legendtitle,
                                breaks = c("Degraded", "Not degraded", "Target"),
                                values = bio_fill, drop = TRUE) +
    ggplot2::scale_fill_manual(name = legendtitle,
                               breaks = c("Degraded", "Not degraded", "Target"),
                               values = bio_fill, drop = TRUE) +
    ggplot2::scale_shape_manual(name = legendtitle,
                                breaks = c("Degraded", "Not degraded", "Target"),
                                values = bio_shape, drop = TRUE) +
    ggplot2::scale_alpha_manual(name = legendtitle,
                                breaks = c("Degraded", "Not degraded", "Target"),
                                values = bio_alpha, drop = TRUE) +
    ggplot2::labs(title = str_title, subtitle = str_subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 8),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8),
                   axis.ticks.y = ggplot2::element_blank())

  p_tv <- p_tv +
    ggplot2::geom_text(data = str_scores, size = 8/2.845, hjust = 1.5, vjust = 5,
                       ggplot2::aes(x = Sensitivity, y = OverallMax * 0.75, label = Scores))

  if (boo_plot) {
    ggplot2::ggsave(filename = file.path(plot.path, stressor, fn_png_p1),
                    plot = p_tv, dpi = plotdpi, width = plotW,
                    height = plotH, units = plotunits)
  }
}
