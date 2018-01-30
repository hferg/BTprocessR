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
  
  d <- reshape2::melt(d)
  p <- ggplot2::ggplot(d, ggplot2::aes(value, fill = variable)) +
    ggplot2::geom_histogram() +
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

plot.bt_post <- function(x, ...) {
  autoplot(x, ...)
}

#' Print function for the S3 class trees_summary
#'
#' @name print.trees_summary

print.trees_summary <- function(x) {
  n <- length(x)
  nm <- names(x)
  cat(n, "phylogenetic trees\n")
  for (i in seq_along(x)) {
    cat("tree", i, ":", nm[i], length(x[[i]]$tip.label), "tips\n")
  }
}

#' A plot method for the class "trees_summary" - the default will
#' plot all trees in a 2x2 plot, otherwise a tree(s) can be specified
#' to plot alongside the time-tree.
#' @name plot.trees_summary

plot.trees_summary <- function(x, tree = NULL, tips = FALSE) {
  if (is.null(tree)) {
    nms <- names(x)
    par(mfrow = c(2, 2))
    for (i in seq_along(x)) {
      plot(x[[i]], show.tip.label = tips, main = nms[i])
    }
  } else {
    if (!tree %in% c("mean_tree", "median_tree", "mode_tree")) {
      stop("'tree' must be either 'mean_tree', 'median_tree' or 'mode_tree'.")
    }
    par(mfrow = c(1, 2))
    plot(x$original_tree, show.tip.label = tips, main = names(x)[1])
    plot(x[[tree]], show.tip.label = tips, main = tree)
  }
}
