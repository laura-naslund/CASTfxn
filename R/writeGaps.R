#' Write gap statement when co-occurrence eliminates all candidate causes
#'
#' @return
#' @export
#'
#' @examples
writeCoOccurElimGap <- function(TargetSiteID, bioComm){
  msg <- paste0("No candidate causes to evaluate further for ",
                TargetSiteID, " for the ", bioComm, " community.")
  message(msg)

  gapcomment <- paste0(bioComm, ": All candidate causes were eliminated by the co-occurrence line of evidence for ",
                       TargetSiteID)

  gap.statement <- data.frame(
    fxnname = "getCoOccur",
    condition = msg,
    result = as.character(0),
    comment = gapcomment)


  temp_status <- data.frame(TargetSiteID = as.character(TargetSiteID),
                            status = "Failed",
                            reason = "No stressors remaining after removing stressors with no paired stressor-response samples at the target site, stressors with insufficient number of samples at comparator sites, and co-occurrence screening. ")

  return(
    list(gap.statement = gap.statement,
         temp_status = temp_status)
  )
  }

#' Write gap statement when there are no SSTVs
#'
#' @return
#' @export
#'
#' @examples
writeSSTVGap <- function(TargetSiteID){
  msg <- "No site stressors have stressor-specific tolerance values"
  message(msg)

  gap.statement <- data.frame(
    fxnname = "getVerifiedPredictions",
    condition = TargetSiteID,
    result = "0",
    comment = msg
  )

  return(gap.statement)
}


#' Write gap statement for SSI missing from inclusion in response metadata
#'
#' @param index
#' @param bioComm
#'
#' @return
#' @export
#'
#' @examples
writeMissingSSIGap <- function(index, bioComm){
  msg <- paste0(paste(unique(index), collapse = ", "),
                " are listed as stressor-specific index in stressor metadata but are not found as a metric to be included in ",
                bioComm,
                " response metadata.")
  message(msg)

  gap.statement <- data.frame(
    fxnname = "getVPSSIscores",
    condition = "Missing stressor-specific index information",
    result = unique(index),
    comment = msg
  )

  return(gap.statement)
}

#' Write gap statement for no SSIs
#'
#' @param TargetSiteID
#'
#' @return
#' @export
#'
#' @examples
writeSSIGap <- function(TargetSiteID){
  msg <- "No site stressors have stressor-specific indices"
  message(msg)

  gap.statement <- data.frame(
    fxnname = "getVPSSIscores",
    condition = TargetSiteID,
    result = "0",
    comment = msg
  )

  return(gap.statement)
}
