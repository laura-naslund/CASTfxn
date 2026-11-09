#' Advance progress through the Shiny Run Report progress bar
#'
#' @param config list of configuration parameters
#' @param msg message to display in the progress bar
#' @param prog_cnt current step of progress bar
#' @param prog_n total number of progress bar steps
#' @param prog_sleep length of system sleep before advancing, intended to give users time to read progress
#'
#' @return updated progress bar step count
#' @export
stepProgress <- function(config, msg, prog_cnt, prog_n, prog_sleep){
  if(isTRUE(config$boo_shiny)){
    prog_cnt <- prog_cnt + 1
    shiny::incProgress(1/prog_n, message = paste0("Step ", prog_cnt), detail = msg)
    Sys.sleep(prog_sleep)
    message("Step: ", prog_cnt, "; ", msg)
  }
  return(prog_cnt)
}
