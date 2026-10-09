#' Set global parameters
#'
#' @param boo_shiny boolean to run the code in the Shiny application
#' @param boo_debug boolean to run the application in debug mode
#' @param dn_checked_sk name of the checked inputs folder
#' @param boo_plot boolean to generate plots
#' @param pt_type vector of the types of points plotted (e.g., target site, inside the case degraded)
#' @param pt_fill vector of point fill colors corresponding to pt_type
#' @param pt_shape vector of point shapes corresponding to pt_type
#' @param pt_size vector of point sizes corresponding to pt_type
#' @param pt_alpha vector of point transparencies corresponding to pt_type
#' @param refOutline_col hex code for the color of point outlines
#' @param plot_dpi resolution of plots in dots per inch
#' @param plot_units units corresponding to plot dimensions
#' @param plot_H height of plots
#' @param plot_W width of plots
#'
#' @return list of configuration parameters
#' @export
setGlobals <- function(boo_shiny = FALSE,
                       boo_debug = FALSE,
                       dn_checked_sk = "_CheckedInputs",
                       boo_plot = TRUE,
                       pt_type = c("target", "insideND", "insideD", "outsideND", "outsideD"),
                       pt_fill = c("#CC79A7", "#56B4E9", "#E69F00", "#0072B2", "#D55E00"),
                       pt_shape = c(24, 21, 25, 21, 25),
                       pt_size = c(1.75, 0.8, 1, 0.8, 1),
                       pt_alpha = c(1, 0.5, 0.7, 0.5, 0.7),
                       refOutline_col = "#000000",
                       plot_dpi = 600,
                       plot_units = "in",
                       plot_H = 6,
                       plot_W = 8) {
  data_plotvars <- data.frame(
    "Type" = pt_type,
    "Fill" = pt_fill,
    "Shape" = pt_shape,
    "Size" = pt_size,
    "Alpha" = pt_alpha
  )

  return(list(
    boo_shiny = boo_shiny,
    boo_debug = boo_debug,
    dn_checked_sk = dn_checked_sk,
    boo_plot = boo_plot,
    data_plotvars = data_plotvars,
    refOutline_col = refOutline_col,
    plot_dpi = plot_dpi,
    plot_units = plot_units,
    plot_H = plot_H,
    plot_W = plot_W
  ))

}
