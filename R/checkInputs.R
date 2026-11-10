#' Check input files
#'
#' @param dir_uploaded path to data input directory
#' @param dir_out path to results directory
#' @param fn.inputcheck path to input check file
#' @param fn.matchups path to matchups file
#'
#' @return a list containing the two check input dataframes
#' @export
checkInputs <- function(dir_uploaded,
                        dir_out,
                        fn.inputcheck = system.file("extdata",
                                                    "CASTool_InputCheck.xlsx",
                                                    package = "CASTfxn"),
                        fn.matchups = system.file("extdata",
                                                  "CASTool_MatchUps.xlsx",
                                                  package = "CASTfxn")) {

  # QC----
  if (!dir.exists(dir_uploaded)) {
    stop("ERROR: 'dir_uploaded' does not exist.")
  }

  if (!dir.exists(dir_out)) {
    stop("ERROR: 'dir_out' does not exist.")
  }

  if (!file.exists(fn.inputcheck)) {
    stop("ERROR: 'fn.inputcheck' does not exist.")
  }

  # Read and combine inputcheck and CASTool_Metadata files ----
  input_check <- readxl::read_xlsx(fn.inputcheck, sheet = "data", na = "")
  paired_check <- readxl::read_xlsx(fn.inputcheck, sheet = "paired", na = "") |>
    tibble::rowid_to_column("UID")

  metadata_file <- list.files(dir_uploaded) |> grep("_CASTool_Metadata", x = _, value = TRUE)
  fn.metadata <- file.path(dir_uploaded, metadata_file)

  if (file.exists(fn.metadata)) {
    CASTmetadata <- readCASToolData(fn = fn.metadata, NAs = c("", "NA"))
    pf.cols.CASTmeta <- compare.colnames(c("Variable", "Value"),
                                         colnames(CASTmetadata))

    if(length(setdiff(c("Variable", "Value"), colnames(CASTmetadata))) != 0){
      stop(pf.cols.CASTmeta)
    } else {
      CASTmetadata <- dplyr::select(CASTmetadata, Variable, Value)
    }

  } else {
    msg <- "_CASTool_Metadata is not uploaded or otherwise does not exist."
    stop(msg)
  }

  input_check <- merge(input_check, CASTmetadata, by.x = "FilePath",
                       by.y = "Variable", all.x = TRUE)
  input_check$Uploaded <- 0

  cluster_graphic.fn <- CASTmetadata |>
    dplyr::filter(Variable == "fn.cluster.graphic") |>
    dplyr::pull(Value)

  # Check fatal errors ----
  ## Check definition status ----
  # All required files must be defined. If not, generate a fatal error.
  req.not.defined <- input_check |>
    dplyr::filter(is.na(Value) & stringr::str_detect(Type, "required")) |>
    dplyr::pull(FilePath)

  if(length(req.not.defined) > 0){
    msg <- paste0(paste(req.not.defined, collapse = ", "), " are not defined and likely not uploaded.")
    stop(msg)
  }

  ## Check load status ----
  helperImport <- as.logical(CASTmetadata$Value[CASTmetadata$Variable == "helperImport"])

  if(helperImport == TRUE){
    input_check_cluster <- input_check |> dplyr::filter(Object == "data_cluster")
  } else{
    input_check_cluster <- NULL
  }

  input_check <- input_check |>
    dplyr::filter(!is.na(Value))|>
    dplyr::mutate(Uploaded = ifelse(file.exists(file.path(dir_uploaded, Value)), 1, 0))

  loaded <- input_check |> dplyr::filter(Uploaded == 1)
  not.loaded <- input_check |> dplyr::filter(Uploaded == 0)

  ## Check req'd files loaded ----
  req.not.loaded <- not.loaded |>
    dplyr::filter(stringr::str_detect(Type, "required")) |>
    dplyr::pull(FilePath)

  if (length(req.not.loaded) > 0) {
    msg <- paste0(paste(req.not.loaded, collapse = ", "), " are not loaded.")
    stop(msg)
  }

  opt.not.loaded <- not.loaded |>
    dplyr::filter(stringr::str_detect(Type, "optional")) |>
    dplyr::pull(FilePath)

  if (length(opt.not.loaded) > 0) {
    msg <- paste0(paste(opt.not.loaded, collapse= ", "), " are not loaded.")
    warning()
  }

  ## Check contingent files (paired) ----
 missing_paired <-  paired_check |>
    dplyr::filter(DataFile %in% input_check$FilePath) |>
    tidyr::pivot_longer(cols = c(DataFile, MetadataFile)) |>
    dplyr::left_join(input_check |> dplyr::select(FilePath, Uploaded), by = c("value" = "FilePath"))|>
    dplyr::filter(Uploaded == 0) |>
    dplyr::pull(value)

  if(length(missing_paired) > 0){
    msg <- paste0("Missing: ", paste(missing_paired, collapse = ", "))
    stop(msg)
  }

  # Read all loaded files ----
  ## Could have better syntax but leaving for now
  df.reqd.obj.cols <- data.frame("Object" = as.character(),
                                 "Result" = as.character())

  obj.list <- list()

  for (j in 1:nrow(loaded)) {

    object <- loaded$Object[j]
    fn <- loaded$Value[j]

    if(tolower(tools::file_ext(fn)) != "rda"){
      assign(object, readCASToolData(file.path(dir_uploaded, fn),
                                     NAs = c("", "NA", "na", "N/A", "n/a")))
    } else {
      readCASToolData(file.path(dir_uploaded, fn),
                      NAs = c("", "NA", "na", "N/A", "n/a"))
    }

    obj.list[[object]] <- get(object)

    ## Check required columns ----
    if(is.na(loaded$ReqCols[loaded$Object == object])){
      df.reqd.obj.cols <- rbind(df.reqd.obj.cols,
                                data.frame(object = object,
                                           reqd.v.act.cols = "Missing columns: "))
      df.coltypes.tmp <- data.frame(object = object, col = NA, errors = "pass")
      df.reqd.coltypes <- rbind(df.reqd.coltypes, df.coltypes.tmp)

      next
    } else{
      reqd.cols <- stringr::str_split_1(loaded$ReqCols[loaded$Object == object], ", ")
    }

    if ("*SSTV" %in% reqd.cols) {
      next # These will be checked below
    }
    if ("!IncaseCol" %in% reqd.cols) {
      reqd.cols <- sub("!IncaseCol",
                       CASTmetadata$Value[CASTmetadata$Variable == "incaseColName"],
                       reqd.cols)
    }
    #act.cols <- colnames(get(object))
    act.cols <- colnames(obj.list[[object]])

    if(length(intersect(act.cols, c("GEOMETRY", "Geometry", "Shape", "SHAPE", "shape")))>0){
      act.cols <- c(act.cols, "geometry")
    }

    reqd.v.act.cols <- compare.colnames(reqd.cols, act.cols)


    if (object %in% c("data_stressorWS", "data_stressorinfoWS")) {
      df.coltypes.tmp <- cbind("object" = object,
                               "col" = "all columns",
                               "errors" = "were generated by code")
    } else {
      df.obj <- obj.list[[object]] |>
        dplyr::select(dplyr::any_of(reqd.cols))

      df.coltypes.tmp <- compare.coltypes(object, df.obj)
    }
    if (j == 1) {
      df.reqd.obj.cols <- data.frame(cbind(object, reqd.v.act.cols))
      df.reqd.coltypes <- df.coltypes.tmp
    } else {
      df.reqd.obj.cols <- rbind(df.reqd.obj.cols,
                                data.frame(cbind(object, reqd.v.act.cols)))
      df.reqd.coltypes <- rbind(df.reqd.coltypes, df.coltypes.tmp)
    }

    # add missing columns
    temp_missing <- setdiff(reqd.cols, act.cols)

    if(length(temp_missing)>0){
      tempdf <- obj.list[[object]]
      tempdf[temp_missing] <- replicate(length(temp_missing), NA)
      assign(object, tempdf)
    }
  }


  # Import files from helper packages
  region <- as.character(CASTmetadata$Value[CASTmetadata$Variable == "region"])

  # Check potential file path lengths ----
  longest_siteid <- data_Sites |>
    dplyr::mutate(char_no = nchar(StationID)) |>
    dplyr::arrange(dplyr::desc(char_no)) |>
    dplyr::slice(1) |>
    dplyr::pull(StationID)

  stressor_metadata <- loaded |>
    dplyr::filter(FilePath %in% c("fn.timestamped.info", "fn.undated.info"))

  stressors_df <- data.frame(stressor = character())
  for(i in seq_len(nrow(stressor_metadata))){
    temp_stress <- get(stressor_metadata$Object[i]) |>
      dplyr::filter(UseInStressorID == 1) |>
      dplyr::rename("stressor" = "StdParamName") |>
      dplyr::select(stressor)

    stressors_df <- stressors_df |>
      dplyr::bind_rows(temp_stress)
  }

  responses_metadata <- loaded |>
    dplyr::filter(FilePath %in% c("fn.alg.metrics.info", "fn.bmi.metrics.info", "fn.fish.metrics.info"))

  responses_df <- data.frame(response_metric = character(), biocomm = character())
  for(i in seq_len(nrow(responses_metadata))){
    temp_responses <- get(responses_metadata$Object[i]) |>
      dplyr::filter(UseYN == "Y") |>
      dplyr::mutate(biocomm = responses_metadata$FilePath[i] |>
                      stringr::str_extract("(?<=fn\\.)\\w+(?=\\.metrics.info)") |>
                      toupper()) |>
      dplyr::rename("response_metric" = "MetricName") |>
      dplyr::select(response_metric, biocomm)

    responses_df <- responses_df |>
      dplyr::bind_rows(temp_responses)
  }

  prob_fp_df <- dplyr::cross_join(responses_df, stressors_df) |>
    dplyr::mutate(dir_out = dir_out,
                  region = region,
                  siteid = longest_siteid,
                  file_path = paste(dir_out,
                                    region,
                                    longest_siteid,
                                    biocomm,
                                    stressor,
                                    paste(longest_siteid,
                                          stressor,
                                          biocomm,
                                          response_metric,
                                          "GI.png",
                                          sep = "_"),
                                    sep = "/")) |>
    dplyr::mutate(num_char = nchar(file_path)) |>
    dplyr::arrange(dplyr::desc(num_char)) |>
    dplyr::filter(num_char > 245)

  ### Write region directory
  if(dir.exists(file.path(dir_out, region))==FALSE){
    dir.create(file.path(dir_out, region))
  }

  fp_problem <- file.path(dir_out, region, "Potential_Problem_FilePath.csv")
  utils::write.csv(prob_fp_df, fp_problem, row.names = FALSE)

  n_problem_fp <- prob_fp_df |>
    nrow()

  msg <- paste0("Number of potentially problematic file path lengths: ", n_problem_fp)
  message(msg)

  if(n_problem_fp > 0){
    msg <- paste0("Review ", fp_problem, " for potentially problematic file paths. Follow directions in the user guide to shorten file paths.")
    message(msg)
  }

 # Import data from helper package ----
  if(helperImport == TRUE){
    conusStates <- setdiff(datasets::state.name, c("Alaska", "Hawaii"))

    regionAvailable <- region %in% conusStates

    if(regionAvailable == TRUE){
      msg <- "Downloading files for requested region from the CASTool helper package."

      clusterNum <- as.character(CASTmetadata$Value[CASTmetadata$Variable == "clusterNumber"])
      if(is.na(clusterNum) == TRUE){
        msg <- "Number of clusters requested from the CASTool helper packages is required. Please add to metadata"
        stop(msg)
      }

      boundary <- CASToolHelperPckg::getBoundary(region)
      reaches <- CASToolHelperPckg::getReaches(region)
      data_cluster <- CASToolHelperPckg::getClusterData(region, clusterNum)
      cluster_graphic <- CASToolHelperPckg::getClusterFig(region, clusterNum)

      obj.list[["data_cluster"]] <- data_cluster

      extra.coltypes <- data.frame(
        object = c(
          "boundary",
          "reaches",
          "data_cluster",
          "data_cluster",
          "cluster_graphic"
        ),
        col = c("geometry", "geometry", "COMID", "ClusterID", NA),
        errors = rep("pass", 5)
      )

      df.reqd.coltypes <- df.reqd.coltypes |> dplyr::bind_rows(extra.coltypes)

      extra.objcols <- data.frame(
        object = c("boundary", "reaches", "data_cluster", "cluster_graphic"),
        reqd.v.act.cols = rep("Missing columns: ", 4)
      )

      df.reqd.obj.cols <- df.reqd.obj.cols |> dplyr::bind_rows(extra.objcols)

      message(msg)
    } else{
      msg <- "Requested region is not available in the CASTool helper packages. Check that region name is identical to a value in https://github.com/laura-naslund/CASToolBaseDataPckg?tab=readme-ov-file#currently-available-regions. If the region is not available, please provide a boundary, reaches, cluster, and cluster graphic file."
      stop(msg)
    }

  } else{
    if(exists("boundary") == FALSE){
      msg <- "Boundary file not found. Please provide."
      stop(msg)
    }
    if(exists("reaches") == FALSE){
      msg <- "Reaches file not found. Please provide."
      stop(msg)
    }
    if(exists("data_cluster") == FALSE){
      msg <- "Cluster file not found. Please provide."
      stop(msg)
    }
    if(exists("cluster_graphic") == FALSE){
      msg <- "Cluster graphic file not found. Please provide."
      stop(msg)
    }
  }

  MT_files <- c("data_bmiMasterTaxa", "data_algMasterTaxa", "data_fishMasterTaxa")

  ## TODO check if this breaks if you have some SSTV in modeled data and some in measured
  if (exists("data_tsInfo")) {
    for(i in seq_along(MT_files)){
      if(exists(MT_files[i])){
        object <- get(MT_files[i])
        type <- paste0("SSTVname.", stringr::str_extract("data_bmiMasterTaxa", "(?<=_)\\w+(?=M)"))
        meas.sstv.bmi <- data_tsInfo |>
          dplyr::filter(!is.na(.data[[type]])) |>
          dplyr::distinct(.data[[type]]) |>
          dplyr::pull(.data[[type]])
        reqd.v.act.cols <- compare.colnames(c("TaxonID", meas.sstv.bmi), names(object))
        df.reqd.obj.cols <- rbind(df.reqd.obj.cols,
                                  data.frame(object = MT_files[i], reqd.v.act.cols = reqd.v.act.cols))
        }
      else{
        next
        }
      }
  }
  if (exists("data_undInfo")) {
    for(i in seq_along(MT_files)){
      if(exists(MT_files[i])){
        object <- get(MT_files[i])
        type <- paste0("SSTVname.", stringr::str_extract("data_bmiMasterTaxa", "(?<=_)\\w+(?=M)"))
        mod.sstv.bmi <- data_undInfo |>
          dplyr::filter(!is.na(.data[[type]])) |>
          dplyr::distinct(.data[[type]]) |>
          dplyr::pull(.data[[type]])
        reqd.v.act.cols <- compare.colnames(c("TaxonID", mod.sstv.bmi), names(object))
        df.reqd.obj.cols <- rbind(df.reqd.obj.cols,
                                  data.frame(object = MT_files[i], reqd.v.act.cols = reqd.v.act.cols))
      }
      else{
        next
      }
    }
  }

  # Complete expected columns and datatypes ----
  df.reqd.coltypes <- df.reqd.coltypes |>
    tidyr::unite(ExpectedDatatypes, col, errors, sep = " ") |>
    dplyr::filter(grepl("pass$", ExpectedDatatypes) == FALSE)

  df.reqd.obj.cols <- merge(df.reqd.obj.cols, df.reqd.coltypes, all.x = TRUE)
  df.reqd.obj.cols <- df.reqd.obj.cols |>
    dplyr::rename(ExpectedColumns = reqd.v.act.cols) |>
    dplyr::mutate(ExpectedColumns = ifelse(ExpectedColumns == "Missing columns: ",
                                           "No missing columns", ExpectedColumns),
                  ExpectedDatatypes = ifelse(is.na(ExpectedDatatypes),
                                             "All expected datatypes confirmed",
                                             ExpectedDatatypes))
  df.reqd.obj.cols <- merge(loaded[, c("Object", "FilePath")], df.reqd.obj.cols,
                            by.x = "Object", by.y = "object", all.y = TRUE)


  # Relational integrity table ----
  match_files <- c(input_check$Object, ifelse(is.null(data_cluster), NULL, "data_cluster"))

  df_matchups <- readxl::read_xlsx(fn.matchups) |>
    dplyr::filter(FileOne %in% match_files & FileTwo %in% match_files)

  rel_int <- dplyr::bind_rows(
      purrr::pmap(df_matchups, buildRelationalIntegrity, obj.list = obj.list)
    )

  ## TableOne: Summary of file inputs ----

  df.TableOne <- df.reqd.obj.cols |>
    dplyr::select(FilePath, Object, ExpectedColumns, ExpectedDatatypes) |>
    dplyr::mutate(FilePath = dplyr::if_else(is.na(FilePath), "helper package", FilePath))

  ## TableTwo: Relational integrity ----
  df.TableTwo <- rel_int


  ## Write results directory ----
  out.folders <- c(dir_out, region, "_CheckedInputs")

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

  dir_out <- dir.path

  # Save objects ----
  saveRDS(CASTmetadata, file.path(dir_out, "CASTmetadata.rds"))
  saveRDS(loaded, file.path(dir_out, "loaded.rds"))

  objects <- setdiff(df.TableOne$Object, "cluster_graphic")

  for (o in seq_along(unique(objects))) {
    objName <- objects[o]
    fn <- paste0(objName, ".rds")
    saveRDS(get(objName), file.path(dir_out, fn))
  }

  if(is.na(cluster_graphic.fn)==FALSE){
    file.copy(file.path(dir_uploaded, cluster_graphic.fn),
              file.path(dir_out, "cluster_graphic.png"))
  } else{

    dest <- file.path(dir_out, "cluster_graphic.png")

    aws.s3::save_object(
      object = cluster_graphic,
      bucket = "dmap-data-commons-ow",
      file = dest
    )
  }

  myTables <- list(TableOne = df.TableOne, TableTwo = df.TableTwo)
  return(myTables)

}

#' Compare column names
#'
#' @param required names of required columns for a file
#' @param actual names of columns in a file
#'
#' @return a statement about whether required columns are missing
#' @export
compare.colnames <- function(required, actual) {
  result <- paste0("Missing columns: ", paste(setdiff(required, actual),
                                              collapse = ", "))
  return(result)
}

#' Compare column types
#'
#' @param object name of dataframe
#' @param df.obj dataframe
#'
#' @return dataframe describing any column type errors encountered for an imported file
#' @export
compare.coltypes <- function(object, df.obj) {

  text.cols <- c("TargetSiteID", "StationID", "StreamCatVar", "Label",
                 "StdParamName", "SourceGroup", "SSIndex", "SSTVname.bmi",
                 "SSTVname.alg", "SSTVname.fish", "DirIncStress", "TaxonID",
                 "MetricName", "MetricLabel", "IndexYN", "UseYN", "TrendWIncStress",
                 "InclusiveIndicator")

  num.cols <- c("COMID", "Latitude", "Longitude", "RefSiteFlag",
                "WatershedValue", "Year", "LogTransf", "UseInStressorID",
                "ResultValue", "NumInd", "CutoffValue")

  flex.cols <- c("ClusterID", "SensMax.bmi", "SensMin.bmi", "SensMax.alg",
                 "SensMin.alg", "SensMax.fish", "SensMin.fish", "StressSampleID",
                 "RespSampleID")

  date.cols <- c("StressSampleDate", "RespSampleDate")

  geo.cols <- c("geometry", "GEOMETRY", "Geometry", "Shape", "SHAPE", "shape")

  accept.blank.cols <- c("SSIndex", "SSTVname.bmi", "SensMax.bmi", "SensMin.bmi",
                         "SSTVname.alg", "SensMax.alg", "SensMin.alg",
                         "SSTVname.fish", "SensMax.fish", "SensMin.fish")

  check_type <- function(col, x) {
    cls <- class(x)

    if (col %in% accept.blank.cols && all(is.na(x))) {
      return("pass")
    }

    if (col %in% geo.cols) {
      return("pass")
    }

    if (col %in% flex.cols) {
      return(if (any(cls %in% c("character", "integer", "numeric"))) {
        "pass"
      } else {
        paste0(col, " is expected to be character or numeric but is ", cls[1])
      })
    }

    if (col %in% num.cols) {
      return(if (any(cls %in% c("integer", "numeric"))) {
        "pass"
      } else {
        paste0(col, " is expected to be numeric but is ", cls[1])
      })
    }

    if (col %in% text.cols) {
      return(if ("character" %in% cls) {
        "pass"
      } else {
        paste0(col, " is expected to be character but is ", cls[1])
      })
    }

    if (col %in% date.cols) {
      parsed <- tryCatch(
        lubridate::parse_date_time(x, orders = c("ymd", "mdy", "dmy")) |>
          lubridate::date(),
        error = function(e) NULL
      )

      return(if (is.null(parsed)) {
        paste0(col, " cannot be converted to date")
      } else {
        "pass"
      })
    }

    "pass"
  }

  out <- lapply(names(df.obj), function(col) {
    data.frame(
      object = as.character(object),
      col = col,
      errors = check_type(col, df.obj[[col]]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)

}

#' Build relational integrity table
#'
#' @param FileOne name of first match up file
#' @param FileTwo name of second match up file
#' @param JoinCol name of shared parameter to facilitate matchup
#' @param obj.list list containing imported files
#'
#' @return the relational integrity table
#' @export
buildRelationalIntegrity <- function(FileOne, FileTwo, JoinCol, obj.list){

  FileOne_df <- obj.list[[FileOne]]
  FileTwo_df <- obj.list[[FileTwo]]

  # if TargetSiteID is the only one
  if("TargetSiteID" %in% names(FileOne_df)){
    FileOne_df <- FileOne_df |> dplyr::rename("StationID" = "TargetSiteID")
  }

  FileOne_unique <- FileOne_df |> dplyr::distinct(.data[[JoinCol]]) |> dplyr::pull(.data[[JoinCol]])

  if(stringr::str_detect(FileTwo, "Metrics$") & JoinCol == "MetricName"){
    FileTwo_unique <- FileTwo_df |> dplyr::select(-c(StationID, RespSampleID, RespSampleDate)) |> names() |> unique()
  } else{
    FileTwo_unique <- FileTwo_df |> dplyr::distinct(.data[[JoinCol]]) |> dplyr::pull(.data[[JoinCol]])
  }

  FileOne_diff <- setdiff(FileTwo_unique, FileOne_unique) |> length()
  FileTwo_diff <- setdiff(FileOne_unique, FileTwo_unique) |> length()

  expected_one <- ifelse(FileOne == "df_targets", "(expected)", "")
  expected_two <- ifelse(FileTwo == "data_Sites", "(expected)", "")

  FileOneCondition <- ifelse(FileOne_diff == 0,
                             paste("All", JoinCol, "are in FileOne"),
                             paste(FileOne_diff, JoinCol, "are not in FileOne", expected_one))

  FileTwoCondition <- ifelse(FileTwo_diff == 0,
                             paste("All", JoinCol, "are in FileTwo"),
                             paste(FileTwo_diff, JoinCol, "are not in FileTwo", expected_two))

  df_ret <- data.frame(FileOne = FileOne,
                       FileTwo = FileTwo,
                       JoinCol = JoinCol,
                       FileOneCondition = FileOneCondition,
                       FileTwoCondition = FileTwoCondition)

  return((df_ret))
}
