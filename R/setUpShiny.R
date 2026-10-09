#' Helper function to configure parameters for Shiny app run
#'
#' @param dir_data path to data input directory
#' @param dir_results path to results directory
#' @param dn_checked_sk name of checked data directory
#'
#' @return list of configuration parameters when run in Shiny
#' @export

setUpShiny <- function(dir_data, dir_results, dn_checked_sk){

  data_CASTmeta_prog <- readRDS(file.path(dir_data, dn_checked_sk, "CASTmetadata.rds"))
  data_CASTmeta_prog <- data_CASTmeta_prog |>
    tidyr::pivot_wider(names_from = Variable, values_from = Value)

  biocommlist_prog <- data_CASTmeta_prog |>
    dplyr::pull(biocommlist) |>
    stringr::str_split(", |,")  |>
    unlist()
  n_biocomm_prog <- length(biocommlist_prog)

  return(
    list(
      in_dir = dir_data,
      out_dir = dir_results,
      region = data_CASTmeta_prog$region,
      useBC = FALSE,
      dn_checked_sk = dn_checked_sk,
      n_biocomm_prog = n_biocomm_prog
    )
  )
}
