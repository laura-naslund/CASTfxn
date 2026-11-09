utils::globalVariables(c(
  "algIndexGp",
  "AlgSampFlag",
  "all_sites",
  "AllSamples",
  "BetterThan",
  "biocomm",
  "bioComm",
  "Biocomm",
  "BioComm",
  "biocommlist",
  "BioDeg",
  "bioIndex",
  "bioIndexGp",
  "bioIndexName",
  "bioMasterTaxa",
  "bioMetricData",
  "bioMetricInfo",
  "bioTaxaData",
  "bmiIndex",
  "bmiIndexGp",
  "boo.plot.user",
  "boo_plot_ref",
  "boo_plot_user",
  "Case",
  "case",
  "cfg",
  "char_no",
  "ClusterID",
  "COMID",
  "comp_sites",
  "Count",
  "CutoffValue",
  "data_algMetrics",
  "data_bioCoOccur",
  "data_bmiCoOccur",
  "data_bmiMetrics",
  "data_bmiMetricsInfo",
  "data_cluster",
  "DataFile",
  "data_fishMetrics",
  "data_measOutliers",
  "data_outliers",
  "data_plotvars",
  "data_respTrim",
  "data_sampSummary",
  "data_Sites",
  "data_Stress",
  "data_stressInfo",
  "data_stressorinfoWS",
  "data_stressorWS",
  "data_tsInfo",
  "data_und_all",
  "data_undInfo",
  "DateDiff",
  "DataFile",
  "df_algMetrics",
  "df_fishMetrics",
  "df_LoE",
  "df_PairedStressResp",
  "df_stressorMetadata",
  "dir.path",
  "dir_results",
  "DirIncStress",
  "errors",
  "estimate",
  "Estimate",
  "ExpectedColumns",
  "ExpectedDatatypes",
  "Evidence",
  "FileOne",
  "FilePath",
  "file_path",
  "FileTwo",
  "FinalLatitude",
  "FinalLongitude",
  "fishIndexGp",
  "FishSampFlag",
  "geometry",
  "getQualSites",
  "Group",
  "GroupName",
  "GrpNm",
  "in.dir",
  "IncaseBad",
  "IncaseCol",
  "IncaseGood",
  "incaseLabel",
  "IncaseSamples",
  "IncaseYN",
  "InclusiveIndicator",
  "Index",
  "IndexYN",
  "InsideCaseSamples",
  "IQRmethod",
  "lab",
  "label",
  "Label",
  "labelval",
  "Latitude",
  "list.bioParamsDEL",
  "list.CompSites",
  "list.SiteSummary",
  "list.Stress",
  "localdir",
  "LoE",
  "LogTransf",
  "Longitude",
  "lwr",
  "Max",
  "max_count",
  "maxIndex",
  "meanIndex",
  "MetadataFile",
  "MetricLabel",
  "MetricName",
  "Min",
  "minIndex",
  "ModeledSampleID",
  "mySiteInfo",
  "n",
  "n_comp",
  "n_site",
  "name",
  "NHD.STATE",
  "not_all_na",
  "NotNA",
  "num_char",
  "NumInd",
  "NumIndeterminate",
  "NumInds",
  "NumNotEvaluated",
  "NumRefute",
  "NumRespSamples",
  "NumStressors",
  "NumStressSamples",
  "NumSupport",
  "NumTaxa",
  "Object",
  "out.dir",
  "OutcaseBad",
  "OutcaseCol",
  "OutcaseGood",
  "outcaseLabel",
  "OutcaseSamples",
  "OutcaseYN",
  "Outlier",
  "OutsideCaseSamples",
  "OverallMax",
  "PctInd",
  "PctInds",
  "PctTaxa",
  "plot_dpi",
  "plot_H",
  "plot_units",
  "plot_W",
  "q25",
  "q50",
  "q75",
  "Quality",
  "ReachLoc",
  "ReferenceSamples",
  "refOutline_col",
  "RefSamples",
  "RefSamplesBad",
  "RefSamplesGood",
  "RefSiteFlag",
  "region",
  "RelAbund",
  "respLabel",
  "Response",
  "response_metric",
  "ResponseValue",
  "RespSampleDate",
  "RespSampleID",
  "respVal",
  "result",
  "ResultValue",
  "rowname",
  "SampleDate",
  "Samples",
  "SampleTotAbund",
  "SampleTotTaxa",
  "SampleType",
  "samplim",
  "Sc_Box",
  "Sc_SRlog",
  "Sc_VPlog",
  "Sc_VPSSI_box",
  "Score",
  "Scores",
  "SDmethod",
  "SensMax",
  "SensMin",
  "Sensitivity",
  "siteDetectsAll",
  "SourceGroup",
  "SRLin_Score_inside",
  "SRLin_Score_outside",
  "SRpred_Deg",
  "SSI",
  "SSIndex",
  "SSIqual",
  "SSIValue",
  "SSTV",
  "sstv.sensall",
  "sstv.sensmaxLabel",
  "SSTVname",
  "STATE.shp",
  "StationID",
  "StationIDType",
  "StationID.x",
  "StationID.y",
  "StationID_Master",
  "StdParamName",
  "StreamCatVar",
  "stressLabel",
  "StressLabel",
  "Stressor",
  "StressorCode",
  "StressorLabel",
  "stressor",
  "stressors",
  "StressorType",
  "StressorValue",
  "StressSampleDate",
  "StressSampleID",
  "stressVal",
  "TargetAboveMedian",
  "TargetValue",
  "TargetValue_Max",
  "TargetValue_Min",
  "TargSum",
  "TaxonID",
  "TolMax",
  "TolMin",
  "TotalSamples",
  "TotalSamplesBad",
  "TotalSamplesGood",
  "TransfResult",
  "TransfValue",
  "TrendWIncStress",
  "type",
  "Type",
  "Uploaded",
  "upr",
  "UseInStressorID",
  "UseYN",
  "V1",
  "value",
  "Value",
  "variable",
  "Variable",
  "VPpred_Deg",
  "WatershedValue",
  "WatershedValueMedian",
  "WatershedValueMedian_Min",
  "WatershedValueMedian_Max",
  "WtTot_WoE",
  "WtTotCO_boxplot",
  "WtTotSR_InCase_LinRegr",
  "WtTotSR_InCase_LogRegr",
  "WtTotSR_OutCase_LinRegr",
  "WtTotTS_barplot",
  "WtTotVP_boxplot",
  "x",
  "xintercept",
  "xval",
  "y",
  "y.name",
  "Year",
  "Year_Min",
  "Year_Max",
  "yLoc",
  "yval"
))

# generic <- stringr::str_replace_all(""," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str1 <- stringr::str_replace_all("data_Stress <- data_stressInfo <- data_respTrim <- StressSampleDate <-
#   StressSampleID <- SourceGroup <- StdParamName <- RespSampleDate <-
#   RespSampleID <- ModeledSampleID <- data_Sites <-
#   StationID <- Label <- Type <- biocomm <- SampleDate", " <- | <-\n  ", ",") |>
#   stringr::str_split(",") |> unlist()
#
# str2 <- stringr::str_replace_all("data_sampSummary <- data_Stress <- StationID <- StdParamName <-
#   TransfResult"," <- | <-\n  ", ",") |> stringr::str_split(",") |> unlist()
#
# str3 <- stringr::str_replace_all("df_stressorMetadata <-
#   bioMetricInfo <- bioMetricData <- df_PairedStressResp <- bioComm <-
#   bioIndex <- plot_dpi <- plot_H <- plot_W <- plot_units <-
#   Type <- StationID <- StressSampleID <- RespSampleID <- Quality <-
#   RefSiteFlag <- Stressor <- Response <- lwr <- upr <- StressSampleDate <-
#   RespSampleDate <- stressVal <- respVal <- n_site <- n_comp <-
#   SRLin_Score_inside <- SRLin_Score_outside <- estimate <- LoE <-
#   bioIndexName <- StressorValue <- Score <- list.SiteSummary <- samplim <-
#   data_plotvars <- refOutline_col <- Estimate"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str4 <- stringr::str_replace_all("data_Sites <- data_cluster <- data_bmiCoOccur <- bmiIndex <- StationID <-
#   OutcaseCol <- IncaseCol <- RespSampleID <- RespSampleDate <- Quality"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str5 <- stringr::str_replace_all("df_PairedStressResp <- siteDetectsAll <- data_stressInfo <-
#   list.CompSites <- bioComm <- bioIndex <- Type <- StdParamName <- Stressor <-
#   StationID <- StressSampleID <- StressSampleDate <- IncaseCol <-
#   OutcaseCol <- RespSampleID <- RespSampleDate <- RefSiteFlag <- BetterThan <-
#   Quality <- StressorValue <- NotNA <- DirIncStress <- LogTransf <- Label <-
#   incaseLabel <- Sc_Box <- bioIndexName <- LoE <- Score <- SSIndex"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str6 <- stringr::str_replace_all("data_Sites <- data_Stress <- data_bmiMetrics <- bmiIndex <-
#   StressSampleDate <- AlgSampFlag <- FishSampFlag <- StdParamName <-
#   TransfResult <- RespSampleDate <- RespSampleID <- StationID.y <-
#   StationID.x <- DateDiff <- StressSampleID <- StationID <- Quality <-
#   BioComm"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str7 <- stringr::str_replace_all("comp_sites <- all_sites <- data_bioCoOccur <- stressors <- bioMetricData <- data_stressInfo <- StationID <- StdParamName"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str8 <- stringr::str_replace_all("data_Stress <- data_stressInfo <- out.dir <- StdParamName <- LogTransf <-
#   ResultValue <- TransfResult <- IQRmethod <- SDmethod"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str9 <- stringr::str_replace_all("bioComm <- data_bioCoOccur <- bioIndex <- list.CompSites <- StationID <-
#   Quality <- IncaseCol <- OutcaseCol <- StressSampleDate <- RespSampleDate <-
#   StressSampleID <- RespSampleID <- RefSiteFlag <- IncaseYN <- OutcaseYN <-
#   BetterThan <- IncaseSamples <- IncaseGood <- IncaseBad <- OutcaseSamples <-
#   OutcaseGood <- OutcaseBad <- TotalSamples <- TotalSamplesGood <-
#   TotalSamplesBad <- RefSamples <- RefSamplesGood <- RefSamplesBad <-
#   rowname <- V1 <- Label <- getQualSites <- InsideCaseSamples <-
#   OutsideCaseSamples <- AllSamples <- ReferenceSamples <- BioComm"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str10 <- stringr::str_replace_all("list.CompSites <- data_Sites <- data_sampSummary <-
#   data_bmiMetrics <- bmiIndexGp <- data_algMetrics <- algIndexGp <-
#   data_fishMetrics <- fishIndexGp <- outcaseLabel <- incaseLabel <- in.dir <-
#   region <- StationID <- Latitude <- Longitude <- RefSiteFlag <- COMID <-
#   OutcaseCol <- IncaseCol <- SampleDate <- SampleType <- yLoc <-
#   df_algMetrics <- df_fishMetrics <- Quality <- Index <- Samples <-
#   RespSampleID <- RespSampleDate <- Score <- Case"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str11 <- stringr::str_replace_all("STATE.shp <- NHD.STATE <- data_Sites <- list.CompSites <- data_plotvars <-
#   refOutline_col <- StationID <- Longitude <- Latitude <- RefSiteFlag <-
#   COMID <- IncaseCol <- OutcaseCol <- Case <- geometry"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str12 <- stringr::str_replace_all("Variable <- StreamCatVar <- Year <- ReachLoc <- COMID <- ClusterID <- Label"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str13 <- stringr::str_replace_all("list.CompSites <- siteDetectsAll <- data_Stress <- data_stressInfo <- list.bioParamsDEL <- data_plotvars <- StationID <- StdParamName <- StressSampleID <- StressSampleDate <- ResultValue <- TransfResult <- RefSiteFlag <- GroupName <- GrpNm <- value <- Label <- boo_plot_ref <- Biocomm <- Stressor <- LogTransf <- SSTV <- SSI <- SensMin <- SensMax <- TolMin <- TolMax <- DirIncStress <- SSTVname <- SSIndex"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str14 <- stringr::str_replace_all("df_PairedStressResp <- list.CompSites <- df_stressorMetadata <- bioComm <-
#   bioIndex <- dir_results <- boo_plot_user <- Type <- StationID <-
#   StressSampleID <- StressSampleDate <- RespSampleID <- RespSampleDate <-
#   Quality <- Stressor <- LogTransf <- Label <- y <- y.name <- x <-
#   StressorCode <- StressorValue <- n <- SRpred_Deg <- Sc_SRlog <-
#   bioIndexName <- LoE <- Score <- BioComm"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str15 <- stringr::str_replace_all("localdir <- StationID_Master <- FinalLatitude <- FinalLongitude <- BioComm <- BioDeg <- NumRespSamples <- minIndex <- meanIndex <- maxIndex <- StressorType <- NumStressSamples <- NumStressors <- WtTot_WoE <- WtTotTS_barplot <- WtTotCO_boxplot <- WtTotSR_InCase_LogRegr <- WtTotSR_InCase_LinRegr <- WtTotSR_OutCase_LinRegr <- WtTotVP_boxplot"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str16 <- stringr::str_replace_all("bioComm <- bioIndex <- data_Stress <- bioMetricData <- data_bmiMetricsInfo <-
#   df_stressorMetadata <- StressSampleDate <- StationID <- RespSampleDate <-
#   StressSampleID <- RespSampleID <- Quality <- Value <- type <- SampleDate <-
#   bioIndexName <- StressorValue <- Stressor <- LogTransf <- Label"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str17 <- stringr::str_replace_all("df_stressorMetadata <- df_PairedStressResp <- bioComm <- bioTaxaData <-
#   bioMasterTaxa <- bioIndexGp <- dir_results <- boo.plot.user <- SSTVname <-
#   SensMin <- SensMax <- Stressor <- LogTransf <- Label <- TaxonID <-
#   IncaseYN <- Quality <- StationID <- IncaseCol <- StressSampleID <-
#   StressSampleDate <- RespSampleID <- RespSampleDate <- RefSiteFlag <-
#   BetterThan <- sstv.sensall <- sstv.sensmaxLabel <- Group <- NumInd <-
#   PctInd <- PctTaxa <- NumInds <- PctInds <- NumTaxa <- variable <- value <-
#   Min <- q25 <- q50 <- q75 <- Max <- StressorLabel <- StressorValue <-
#   Response <- ResponseValue <- Score <- Scores <- Type <- OverallMax <-
#   bioIndex <- bioIndexName <- LoE <- BioComm"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str18 <- stringr::str_replace_all("df_stressorMetadata <- df_PairedStressResp <- bioComm <- bioMetricData <-
#   bioMetricInfo <- bioIndex <- dir_results <- boo_plot_user <- Type <-
#   Stressor <- SSIndex <- Label <- LogTransf <- MetricName <- StationID <-
#   StressSampleID <- StressSampleDate <- RespSampleID <- RespSampleDate <-
#   IncaseYN <- IncaseCol <- Quality <- SSIValue <- SSIqual <- StressorValue <-
#   n <- VPpred_Deg <- Sc_VPlog <- y <- y.name <- x <- Sc_VPSSI_box <- Min <-
#   q25 <- q50 <- q75 <- Max <- TransfValue <- incaseLabel <- bioIndexName <-
#   LoE <- Score <- BioComm"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str19 <- stringr::str_replace_all("bioComm <- df_LoE <- Score <- StressorValue <- StationID <- StressSampleID <-
#   StressSampleDate <- RespSampleID <- RespSampleDate <- Stressor <- LoE <-
#   NumSupport <- NumRefute <- NumIndeterminate <- NumNotEvaluated <-
#   Evidence <- Count <- max_count <- df_stressorMetadata"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str20 <- stringr::str_replace_all("data_stressorWS <- data_stressorinfoWS <-
#   StationID <- SampleDate <- SampleType <- COMID <- WatershedValue <-
#   mySiteInfo <- StreamCatVar <- TargetValue <- Year <- WatershedValueMedian <-
#   yLoc <- Label <- TargetAboveMedian <- list.CompSites <- data_sampSummary <-
#   data_Sites <- IncaseCol"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str21 <- stringr::str_replace_all("bioIndexGp <- not_all_na <- RespSampleDate <- RespSampleID <- StationID <-
#   UseYN <- MetricName <- MetricLabel <- IndexYN <- TrendWIncStress <-
#   CutoffValue <- InclusiveIndicator <-NumInd <-
#   SampleTotAbund <- NumTaxa <- RelAbund <- SampleTotTaxa <- TaxonID <-
#   PctInd <- PctTaxa"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str22 <- stringr::str_replace_all("UseInStressorID <- StdParamName <- StationID <- StressSampleID <-
#   StressSampleDate <- TransfResult <- IQRmethod <- SDmethod <- Outlier"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str23 <- stringr::str_replace_all("UseInStressorID <- StdParamName <- data_measOutliers <- StationID <-
#   StressSampleID <- StressSampleDate <- TransfResult <- IQRmethod <-
#   SDmethod <- Outlier <- data_und_all"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# str24 <- stringr::str_replace_all("data_outliers <- data_stressInfo <- list.CompSites <-
#   StationID <- TransfResult <- Outlier <- StdParamName <- IncaseCol"," <- | <-\n  ", ",")|> stringr::str_split(",") |> unlist()
#
# all_strs <- c(str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11,str12,str13,str14,str15,str16,str17,str18,str19,str20,str21,str22,str23,str24)
#
#  all_strs |>
#        unique() |>
#        sort() |>
#       (\(x) cat(shQuote(x), sep = ",\n"))()
