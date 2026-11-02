#' Get paired stressor-response dataset
#'
#' @param df_sites sites data containing site COMID, comparator assignment, and reference status
#' @param df_stress stressor data in long format
#' @param biocomm biological community indicator
#' @param df_resp biological repsonse data
#' @param index biological data
#' @param lagdays days before and after biological response observation to which a stressor observation can be matched
#'
#' @return dataframe with paired stressor-response observations and comparator site indicator
#' @export
getCoOccurDataset <- function(df_sites,
                              df_stress,
                              biocomm,
                              df_resp,
                              index,
                              lagdays = c(0, 0)) {##FUNCTION.START

  # Debug
  boo_DEBUG <- FALSE

  if  (boo_DEBUG == TRUE) {
    df_sites = data_Sites
    df_stress = data_Stress
    biocomm = "BMI"
    df_resp = data_bmiMetrics
    index = bmiIndex
    lagdays = lagdays
  }

  # define pipe
  not_all_na <- function(x) {!all(is.na(x))}
  biocomm <- tolower(biocomm)

  df_und <- dplyr::filter(df_stress, is.na(StressSampleDate))
  df_ts <- dplyr::filter(df_stress, !is.na(StressSampleDate))

  # Read data files (stressor and response)
  df_resp <- df_resp |> dplyr::select(dplyr::all_of(c("StationID", "RespSampleDate", "RespSampleID", index, "Quality")))

  # Convert undated data to wide format  ----
  if (nrow(df_und) > 0) {
    df_und <- df_und |>
      dplyr::select(StationID, StdParamName, TransfResult) |>
      tidyr::pivot_wider(names_from = StdParamName, values_from = TransfResult)

  } else {
    msg <- "No modeled stressors"
    message(msg)
  }

  # Clean up measured data and convert to wide format ----
  if (nrow(df_ts) > 0) {
    df_ts <- df_ts |>
      tidyr::pivot_wider(names_from = StdParamName, values_from = TransfResult)
  } else {
    msg <- "No measured stressors"
    message(msg)
  }

  # Merge undated stressor data and response data (date not required)
  if(nrow(df_und) > 0){
    df_coOccur <- merge(df_resp, df_und, by.x = "StationID", by.y = "StationID", all = TRUE)
  } else{
    df_coOccur <- df_resp
  }

  # Merge timestamped stressor data and response data
  if(nrow(df_ts) > 0){
    df_coOccur <- fuzzyjoin::fuzzy_left_join(df_coOccur, df_ts,
                                             by = c("StationID" = "StationID",
                                                    "RespSampleDate" = "StressSampleDate"),
                                             match_fun = list(`==`, function(x, y)
                                               (y - x >= -lagdays[1] &
                                                  y - x <= lagdays[2]))) |>
      dplyr::filter(!is.na(StationID.y)) |> # no stressor match
      dplyr::rename(StationID = StationID.x) |>
      dplyr::mutate(DateDiff = abs(as.numeric(RespSampleDate - StressSampleDate))) |> # deal with multiple matching dates
      dplyr::group_by(RespSampleID) |>
      dplyr::arrange(DateDiff) |>
      dplyr::slice(1) |>
      dplyr::select(-DateDiff) |>
      dplyr::ungroup()
  }

  df_coOccur <- df_coOccur |>
    dplyr::mutate(BioComm = biocomm) |>
    dplyr::select(-StationID.y) |>
    dplyr::relocate(StationID, StressSampleDate, RespSampleDate, StressSampleID,
             RespSampleID, BioComm, dplyr::all_of(index), Quality) |>
    dplyr::select_if(not_all_na)


  df_sites <- dplyr::select(df_sites, StationID, dplyr::ends_with("caseCol"))
  df_coOccur <- merge(df_sites, df_coOccur)

  return(df_coOccur)

}







