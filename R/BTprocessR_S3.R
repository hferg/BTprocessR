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
#' @method autoplot bt_post
#' @name autoplot.bt_post
#' @export

autoplot.bt_post <- function(posterior, parameters = NULL, col = NULL,
                             pal = "D") {
  d <- posterior[, !colnames(posterior) %in% c("Iteration", "Harmonic.Mean", 
    "Tree.No")]
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
    ggplot2::geom_histogram() +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ggplot2::facet_wrap( ~ variable, scales = "free") +
    ggplot2::theme(
      legend.position = "none"
    )
  p
}

#' A print method for class bt_post
#' Displays some descriptive info about the posterior in a table
#' @method print bt_post
#' @export

print.bt_post <- function(x) {
  z <- x
  class(z) <- class(x)[-1]
  params <- colnames(x)
  params <- params[!params %in% c("Iteration", "Tree.No")]
  st <- matrix(nrow = length(params), ncol = 5)
  colnames(st) <- c("Parameter", "Median", "Mean", "Mode", "SD")
  st[ , 1] <- params
  st[ , 2] <- round(apply(x[, params], 2, mean), 3)
  st[ , 3] <- round(apply(x[, params], 2, median), 3)
  st[ , 4] <- round(apply(x[, params], 2, modeStat), 3)
  st[ , 5] <- round(apply(x[, params], 2, sd), 3)
  cat("Posterior of ", nrow(x), " samples\n\n")
  print(data.frame(st))
  cat("\n")
  print(z)
}

#' A plot method for the class "bt_post" that invokes the
#' autoplot method.
#' @method plot bt_post
#' @export

plot.bt_post <- function(x, ...) {
  autoplot(x, ...)
}

#' Print function for the S3 class trees_summary
#' @method print trees_summary
#' @name print.trees_summary
#' @export

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
#' @method plot trees_summary
#' @export
#' @name plot.trees_summary

plot.trees_summary <- function(x, tree = NULL, tips = FALSE, ...) {
  if (is.null(tree)) {
    nms <- names(x)
    par(mfrow = c(2, 2))
    for (i in seq_along(x)) {
      plot(x[[i]], show.tip.label = tips, main = nms[i], ...)
    }
  } else {
    if (!tree %in% c("mean_tree", "median_tree", "mode_tree")) {
      stop("'tree' must be either 'mean_tree', 'median_tree' or 'mode_tree'.")
    }
    par(mfrow = c(1, 2))
    plot(x$original_tree, show.tip.label = tips, main = names(x)[1], ...)
    plot(x[[tree]], show.tip.label = tips, main = tree, ...)
  }
}

#' An autoplot method for the class "bt_stones". Plots a bayes factor matrix.
#' @method autoplot bt_stones
#' @export
#' @name autoplot.bt_stones

autoplot.bt_stones <- function(x) {
  d <- matrix(ncol = nrow(x), nrow = nrow(x))
  colnames(d) <- rownames(d) <- x$logfile
  
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      d[i, j] <- round(2 * (x$marginalLh[i] - x$marginalLh[j]), 2)
    }
  }

  d <- reshape2::melt(d[nrow(d):1, ])
  d$bf_cat <- cut(d$value, 
    breaks = c(min(d$value) - 1, 0, 2, 5, 10, max(d$value)),
    labels = c("<0", "0-2", "2-5", "5-10", ">10"))
  d$value[d$value == 0] <- d$bf_cat[d$value == 0] <- NA
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Var2, y = Var1, fill = bf_cat)) +
    ggplot2::geom_tile(colour = "white", size = 0.25) +
    ggplot2::geom_text(data = d, ggplot2::aes(Var2, Var1, label = value), 
      na.rm = TRUE) +
    ggplot2::scale_y_discrete(expand=c(0,0))+
    ggplot2::scale_x_discrete(expand=c(0,0), position = "top") +
    ggplot2::scale_fill_manual(values = viridis::viridis(8)[3:7], 
      na.value = "grey90", drop = FALSE) +
    ggplot2::labs(x = "", y = "", fill = "Bayes factor") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_family = "Helvetica")
  p
}

#' A plot method for the class "bt_stones" that invokes the
#' autoplot method.
#' @method plot bt_stones
#' @export
#' @name plot.bt_stones

plot.bt_stones <- function(x) {
  autoplot(x)
}

#' A print method for the class "spkey"
#' @method print spkey
#' @name print.spkey
#' @export

print.spkey <- function(x) {
  len <- length(x)
  tips <- length(x[[1]][[2]])
  cat(paste0("\nDescendant species key for a tree of ", len, " edges and ", tips, " species.\n"))
  cat("Element names correspond to column node_id in scalars table of rjpp output.\n")
}

#' A print method for the class "rjpp".
#' @method print rjpp
#' @name print.rjpp
#' @export

print.rjpp <- function(x) {
  trans <- names(x$origins)
  if ("nodes" %in% trans | "branches" %in% trans) {
    trans <- c("variable rates", trans[!trans %in% c("nodes", "branches")])
  }
  cat("\nBayesTraits reversible-jump output testing the following transformations:\n")
  cat("\n")
  for (i in seq_along(trans)) {
    cat(paste0("\t", trans[i], "\n"))
  }
  cat("\n")
  cat(paste0("A posterior consisting of ", x$niter, " samples,\n"))
  cat(paste0("for a tree of ", length(x$species_key[[1]][[2]]), " species and ", nrow(x$tree_summary$tree_summaries$original_tree$edge), " edges."))
}

#' a plot method for class rjpp.
#' @method plot rjpp
#' @name plot.rjpp
#' @export

plot.rjpp <- function(x, plot.options = list(), ...) {
  plotShifts(x, plot.options = plot.options, ...)
}

#' a plot method for class ancstates
#' @method print ancstates
#' @name print.ancstates
#' @export

print.ancstates <- function(x) {
  nn <- nrow(x$states)
  nt <- length(x$tree$tip.label)
  ntn <- x$tree$Nnode
  print(x$states)
  cat(paste("\nAncestral states for", nn, "nodes of a phylogenetic tree with",
    nt, "tips and", ntn, "internal nodes.\n\n"))
}

#' a print method for class ancstates.
#' @method plot ancstates
#' @name plot.ancstates
#' @export

plot.ancstates <- function(x, plot.options = list(), ...) {
  opts <- list(
    state = "mean",
    palette = "plasma",
    border = "black",
    shape = "circle",
    cex = 2,
    legend.pos = "auto",
    legend = "numeric",
    legend.title = NULL
  )
  opts[names(plot.options)] <- plot.options
  if (opts$shape == "circle") {
    opts$shape <- 21
  } else if (opts$shape == "sqaure") {
    opts$shape <- 22
  } else if (opts$shape == "diamond") {
    opts$shape <- 23
  } else if (opts$shape == "uptriangle") {
    opts$shape <- 24
  } else if (opts$shape == "downtriangle") {
    opts$shape <- 25
  }
  if (is.null(opts$legend.title)) {
    opts$legend.title <- "Trait value"
  }
  pdat <- x$states[, c("node", opts$state)]
  xx <- pdat[opts$state][[1]]
  ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
  lims <- round(c(min(xx), max(xx)), 2)
  # make colours, and a scale, and add to the xx tibble.
  if (length(opts$palette) > 1) {
    pdat$cols <- plotrix::color.scale(xx, extremes = opts$palette)
    scale.cols <- plotrix::color.scale(ss, extremes = opts$palette)
  } else {
    pdat$cols <- colourvalues::colour_values(xx, palette = opts$palette)
    scale.cols <- colourvalues::colour_values(ss, palette = opts$palette)
  }
  if (opts$legend.pos == "auto") {
    pos <- c(0, 
      0, 
      round((max(ape::node.depth.edgelength(x$tree)) / 4), 2),
      round((length(x$tree$tip.label) / 60), 2))
  } else {
    pos <- opts$legend.pos
  }
  plotPhylo(x$tree)
  ape::nodelabels(node = pdat$node, bg = pdat$cols, col = opts$border,
    pch = opts$shape, cex = opts$cex)
  plotrix::color.legend(pos[1], pos[2], pos[3], pos[4],
    lims, rect.col = scale.cols, align = "rb")
  labx <- (pos[1] + pos[3]) / 2
  laby <- pos[4] + strheight("H")
  text(labx, laby, opts$legend.title)
}