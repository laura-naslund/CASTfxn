#' Prep biological response data
#'
#' @param file_dir directory containing checked rds files
#' @param bio biological community indicator
#' @param loaded vector of loaded files
#' @param useBC boolean to use biological similarity filter for comparator assignment
#' @param bioIndex name of index used to determine sample degraded/not-degraded status
#' @param calcRelAbund boolean to calculate relative abundance of taxa in count data
#' @param config list containing the configuration parameters
#'
#' @return list containing biological response data (data_bioMetrics), metadata (data_bioMetricsInfo), count data (if available, data_bioCounts), count metadata (if available, data_bioMasterTaxa)
#' @export
prepRespData <- function(file_dir,
                         bio,
                         loaded,
                         useBC,
                         bioIndex,
                         calcRelAbund,
                         config = NULL) {

  if(is.null(config) == FALSE){
    file_dir <- file.path(config$out_dir, config$dn_checked_sk)
    useBC    <- config$useBC
  }

  not_all_na <- function(x) {!all(is.na(x))}

  # Debug
  boo.debug <- FALSE
  if (boo.debug) {
    bio      <- "bmi"
    loaded   <- loaded
    useBC    <- useBC
    bioIndex <- bioIndexGp
  }

  bioMetrics_fp <- file.path(file_dir, paste0("data_", bio, "Metrics.rds"))
  bioMetricsInfo_fp <- file.path(file_dir, paste0("data_", bio, "MetricsInfo.rds"))

  data_bioMetrics <- readRDS(bioMetrics_fp)
  data_bioMetricsInfo <- readRDS(bioMetricsInfo_fp)

  if ((paste0("data_", bio, "MasterTaxa") %in% loaded) && (paste0("data_", bio, "Counts") %in% loaded)){
    bioMasterTaxa_fp <- file.path(file_dir, paste0("data_", bio, "MasterTaxa.rds"))
    bioCounts_fp <- file.path(file_dir, paste0("data_", bio, "Counts.rds"))

    data_bioMasterTaxa <- readRDS(bioMasterTaxa_fp)
    data_bioCounts <- readRDS(bioCounts_fp)
  }

  if (useBC == TRUE) {
    data_BCdist         <- readRDS(file.path(file_dir, "data_BCdist.rds"))
  } else {
    data_BCdist         <- NULL
  }

  # Get BIO metric data ----
  data_bioMetrics <- data_bioMetrics |>
    dplyr::select_if(not_all_na) |>
    dplyr::mutate(RespSampleDate = lubridate::parse_date_time(RespSampleDate,
                                      orders = c("ymd", "mdy", "dmy")) |>
                    lubridate::date())
  data_bioMetrics <- unique(data_bioMetrics)

  # Get BIO metric info ----
  # Add Quality (based on user-supplied data in list.bio.vars & data_bioMetricsInfo)
  data_bioMetricsInfo <- data_bioMetricsInfo |>
    dplyr::filter(UseYN == "Y") |>
    dplyr::select(MetricName, MetricLabel, IndexYN, TrendWIncStress,
                  CutoffValue, InclusiveIndicator)

  bioIndex.cutval <- data_bioMetricsInfo$CutoffValue[data_bioMetricsInfo$MetricName == bioIndex]
  bioIndex.dir <- data_bioMetricsInfo$TrendWIncStress[data_bioMetricsInfo$MetricName == bioIndex]
  bioIndex.Inc <- data_bioMetricsInfo$InclusiveIndicator[data_bioMetricsInfo$MetricName == bioIndex]

  ## Identify quality for primary index ----
  if (bioIndex.Inc == "\u2264" | bioIndex.Inc == "<=") {
    if (bioIndex.dir == "Dec") { # Lower values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = TRUE,
                                     labels = c("Degraded", "Not degraded"))
    } else { # Higher values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = TRUE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else if (bioIndex.Inc == "\u2265" | bioIndex.Inc == ">=") {
    if (bioIndex.dir == "Dec") { # Lower values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = FALSE,
                                     labels = c("Degraded", "Not degraded"))
    } else { # Higher values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = FALSE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else if (bioIndex.Inc == "<") {
    if (bioIndex.dir == "Dec") { # Lower values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = FALSE,
                                     labels = c("Degraded", "Not degraded"))
    } else { # Higher values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = FALSE,
                                     labels = c("Not degraded", "Degraded"))
    }
  } else {
    if (bioIndex.dir == "Dec") { # Lower values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = TRUE,
                                     labels = c("Degraded", "Not degraded"))
    } else { # Higher values of the index represent more stress
      data_bioMetrics$Quality <- cut(data_bioMetrics[, bioIndex],
                                     breaks = c(-Inf, bioIndex.cutval, Inf),
                                     right = TRUE,
                                     labels = c("Not degraded", "Degraded"))
    }

  } # End cut statements

  # If loaded, get bioCounts ----
  if (exists("data_bioCounts")) {
    ## Get BIO count data ----
    data_bioCounts <- data_bioCounts |>
      dplyr::mutate(RespSampleDate = lubridate::parse_date_time(RespSampleDate,
                                          orders = c("ymd", "mdy", "dmy")) |>
                      lubridate::date()
                    )

    if (calcRelAbund == TRUE) {     # Only write this column if needed
      data_bioCounts <- data_bioCounts |>
        dplyr::group_by(RespSampleID, RespSampleDate) |>
        dplyr::mutate(SampleTotAbund = sum(NumInd, na.rm = TRUE),
                         NumTaxa = dplyr::n(),
                         PctInd = round(NumInd / SampleTotAbund, 5),
                         PctTaxa = round(1 / NumTaxa, 5)) |>
        dplyr::ungroup()
    } else {
      data_bioCounts <- data_bioCounts |>
        dplyr::rename(PctInd = RelAbund) |>
        dplyr::group_by(RespSampleID, RespSampleDate) |>
        dplyr::mutate(NumTaxa = dplyr::n(),
                         PctTaxa = round(1 / SampleTotTaxa, 5)) |>
        dplyr::ungroup()
    } # end if calcRelAbund == TRUE
    data_bioCounts <- data_bioCounts |> dplyr::select(StationID, RespSampleID,
                                    RespSampleDate, TaxonID, NumInd, NumTaxa,
                                    PctInd, PctTaxa)
  }

  # Get Master Taxa file - Done above

  if(exists("data_bioMasterTaxa")==FALSE){
    data_bioMasterTaxa <- NULL
  }
  if(exists("data_bioCounts")==FALSE){
    data_bioCounts <- NULL
  }


  myRespData <- list(data_bioMetrics = data_bioMetrics,
                     data_bioMetricsInfo = data_bioMetricsInfo,
                     data_bioCounts = data_bioCounts,
                     data_bioMasterTaxa = data_bioMasterTaxa)
  return(myRespData)

}
