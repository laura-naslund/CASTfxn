#' Set up R console
#'
#' @param in_dir path to data input directory
#' @param out_dir path to results directory
#' @param region name of analysis region
#' @param useBC boolean to use biological filter
#' @param install_missing boolean to install missing CASToolHelperPckg
#' @param load_local boolean to load skeleton code (CASTfxn) locally
#'
#' @return list of configuration parameters when run in R console
#' @export

setUpRConsole <- function(
    in_dir,
    out_dir,
    region,
    useBC = FALSE,
    install_missing = TRUE,
    load_local = TRUE
) {
  # Check/install remotes
  if (!requireNamespace("remotes", quietly = TRUE)) {
    if (install_missing) {
      install.packages("remotes")
    } else {
      stop("Package 'remotes' is required but not installed.")
    }
  }

  # Check/install CASToolHelperPckg
  if (!requireNamespace("CASToolHelperPckg", quietly = TRUE)) {
    if (install_missing) {
      message("Installing CASToolHelperPckg")
      remotes::install_github("laura-naslund/CASToolHelperPckg")
    } else {
      stop("Package 'CASToolHelperPckg' is required but not installed.")
    }
  }

  # Load helper package
  library(CASToolHelperPckg)

  # Load local package/functions if requested
  if (load_local) {
    devtools::load_all()
  }

  # Return setup values
  return(list(
    in_dir = in_dir,
    out_dir = out_dir,
    region = region,
    useBC = useBC
  ))
}
