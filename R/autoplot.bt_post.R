#' An autoplot (ggplot2) method for objects of class "bt_post".
#' Plots a histogram of all, or a specified number of, parameters
#' from the posterior.
#' 
#' @param posterior An object of class "bt_post", returned from
#' loadPosterior.
#' @param parameter The parameter, or vector of parameter names, to plot. 
#' Must match column names of posterior. If NULL then all parameters 
#' will be plotted.
#' @param col The colour of the bars in the histogram, or 
#' a vector of colours (for multiple plots). If absent the 
#' viridis colour scale is used.
#' @param pal If using the viridis pallettes (i.e. no user-specified
#' colours) then this determines the viridis palletter used.
#' See \link[viridis]{scale_color_viridis}. 
#' Four options are available: "magma" (or "A"), "inferno" (or "B"), 
#' "plasma" (or "C"), and "viridis" (or "D", the default option).
#' @importFrom ggplot2 autoplot
#' @export
#' @name autoplot.bt_post

autoplot.bt_post <- function(posterior, parameters = NULL, col = NULL, 
    pal = "D") {
  d <- posterior[, !colnames(posterior) %in% c("Iteration", "Harmonic.Mean", "Tree.No")]  
  if (!is.null(parameters)) {
    d <- d[, colnames(d) %in% parameters]
  } else {
    parameters <- colnames(d)
  }
  if (is.null(col)) {
    col <- viridis::viridis(length(parameters))
  }
  
  d <- reshape2::melt(d, id.vars = NULL)
  p <- ggplot2::ggplot(d, ggplot2::aes(value, fill = variable)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ggplot2::facet_wrap( ~ variable, scales = "free") +
    theme(
      legend.position = "none"
    )
  p
}

#' A plot method for the class "bt_post" that invokes the 
#' autoplot method.
#' @export

plot.bt_post <- function(x, ...) {
  autoplot(x, ...)
}
