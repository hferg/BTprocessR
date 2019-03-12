##############################################################################
#' rateShifts
#' Generates the edge colours to colour edges by total rate.
#' @name rateShifts
#' @keywords internal
rateShifts <- function(PP, threshold, gradientcols, colour) {
  percscaled <- apply(PP$scalars[[1]][2:nrow(PP$scalars[[1]]), ], 1, function(x) sum(x != 1)) / PP$niter

  if (threshold == 0) {
    edge.cols <- plotrix::color.scale(percscaled, extremes = gradientcols, na.color = NA)

  } else if (threshold > 0) {
    nodes <- as.numeric(names(percscaled[percscaled >= threshold]))
    edge.cols <- rep("black", nrow(PP$meantree$edge))
    edge.cols[PP$meantree$edge[ , 2] %in% nodes] <- colour
  } else if (threshold == "relative") {
    edge.cols <- plotrix::color.scale(log(PP$data$meanRate[2:nrow(PP$data)]), extremes = gradientcols, na.color = NA)
  }
  return(edge.cols)  
}

##############################################################################
#' transShifts
#' Generates the node labels, edge.colours and transparencies to plot the
#' location of transformations in a posterior.
#' @name transShifts
#' @keywords internal

transShifts <- function(PP, threshold, cl, tree, transparency, relativetrans, 
  nodescaling, colour, nodecex) {
  if (threshold == "relative") {
    stop("Relative threshold not valid for total rate scalars.")
  }

  if (threshold == 0) {
    threshold <- 1 / PP$niter
  }
  
  
  nodes <- PP$data$descNode[which((PP$data[ , cl] / PP$niter) >= threshold)]
  

  pprobs <- PP$data[which((PP$data[ , cl] / PP$niter) >= threshold) , cl] / PP$niter

  if (length(nodes) == 0) {
    stop("No scalars above threshold.")
  }

  if (transparency) {
    alphas <- pprobs
  } else {
    alphas <- rep(1, length(nodes))
  }

  if (relativetrans) {
    for (i in 1:length(alphas)) {
      alphas[i] <- (alphas[i] - min(alphas)) / (max(alphas) - min(alphas))
    }
  }

  col <- vector(mode = "character", length = length(nodes))  

  col <- sapply(1:length(alphas), 
    function(x) makeTrans(colour = colour, alpha = alphas[x]))

  if (nodescaling) {
    nodecex = nodecex * pprobs
  }

  if (cl == "nOrgnBRate") {
    nodes <- which(tree$edge[ , 2] %in% nodes)
  }

  list(nodes = nodes, colours = col, alphas = alphas, nodecex = nodecex)
}

################################################################################
#' plotShifts
#' 
#' Plots the locations of the origins of scalars from the postprocessor output of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor - delta/lambda/kappa/node/branch
#' @param threshold Threshold of probability in posterior to display deltas for, defaults to zero (i.e. shows all deltas shaded proportionally to the posterior probability)
#' @param colour The colour to use for the node circles
#' @param scaled Plot the original tree (scaled = "time", the default), or the mean/sclaed tree (scaled = "mean") or plot the tree scaled only by scalars present above the threshold (scaled = "threshold")?
#' @param nodecex The scaling factor for the size of the node circles
#' @param tips Show tip labels?
#' @param scalebar Include scale bar?
#' @param measure When plotting "siginficant" tree, what measure of the parameter? Median (default), mode or mean.
#' @param exludeones If plotting according to a threshold of significance, should 1s (i.e. no scalar) be excluded from the posterior when calculating average scalar?
#' @param relativetrans If TRUE (defaults to FALSE) the scale of transparency will go from the threshold (totally transparent) to the maximum presence (full opacity).
#' @param nodescaling Scale node symbols according to posterior probability of shift (default).
#' @param transparency Adjust node symbol transparency according to posterior probability? Defaults to FALSE.
#' @param gradientcols A vector of two colours - the min and max colours used when colouring the tree according to percentage time rate scaled (when threshold = 0) or using rate.edges.
#' @param rate.edges Takes a numeric value between 0 and 1. If NULL (default) then node shapes are plotted for a transformation. If equal to zero then node shapes are plotted along with a colour gradient on the branches for rates (if in the posterior), and if set to a threshold then the branches are coloured black/red for whether there is a scalar over the threshold (red) along with node scalars.
#' @param shp The shape of the node markers (uses the usual pch index).
#' @name plotShifts
#' @import plotrix
#' @export
#' 
plotShifts <- function(PP, scalar, threshold = 0, threshold2 = 0, nodecex = 2, scaled = "time", scalebar = TRUE,
  measure = "median", excludeones = FALSE, relativetrans = FALSE, nodescaling = TRUE, 
  transparency = FALSE, gradientcols = c("dodgerblue", "firebrick1"), rate.edges = NULL, 
  colour = "red", colour2 = "green",  shp = 21, tips = FALSE, ...) {

  if (scalar == "delta") {
    cl <- "nOrgnDelta"
    mode <- "trans"
  } else if (scalar == "kappa") {
    cl <- "nOrgnKappa"
    mode <- "trans"
  } else if (scalar == "lambda") {
    cl <- "nOrgnLambda"
    mode <- "trans"
  } else if (scalar == "branch") {
    cl <- "nOrgnBRate"
    mode <- "trans"
  } else if (scalar == "node") {
    cl <- "nOrgnNRate"
    mode <- "trans"
  } else if (scalar == "rate") {
    cl <- "nOrgnScalar"
    mode <- "rate"
  } else if (scalar == "nodebranch") {
    cl <- c("nOrgnNRate", "nOrgnBRate")
    mode <- "doubletrans"
  }
  
  if (scaled == "time") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$orgBL[2:nrow(PP$data)]
  } else if (scaled == "mean") {
      tree <- PP$meantree
  } else if (scaled == "median") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$medianBL[2:nrow(PP$data)]
  } else if (scaled == "mode") {
    tree <- PP$meantree
    tree$edge.length <- PP$data$modeBL[2:nrow(PP$data)]
  } else if (scaled == "threshold") {
    if (mode == "rate") {
      tree <- tree
    } else {
      # tree <- shiftScaledTree() #### IN DEV.
    }
  }

  if (mode == "trans") {
    edge.cols <- "black"

    if (isDefined(rate.edges)) {
      try(if(is.null(PP$scalars)) stop("No rate scalars in posterior output."))
      edge.cols <- rateShifts(PP, threshold = rate.edges, gradientcols, colour)
    }
    
    node_info <- transShifts(PP, threshold, cl, tree, transparency, relativetrans,
      nodescaling, colour, nodecex)

  } else if (mode == "doubletrans") {
    edge.cols <- "black"

    if (isDefined(rate.edges)) {
      try(if(is.null(PP$scalars)) stop("No rate scalars in posterior output."))
      edge.cols <- rateShifts(PP, threshold = rate.edges, gradientcols, colour)
    }
    
    node_info <- transShifts(PP, threshold, cl[1], tree, transparency, relativetrans,
      nodescaling, colour, nodecex)

    branch_info <- transShifts(PP, threshold2, cl[2], tree, transparency, relativetrans,
      nodescaling, colour2, nodecex)   

  } else if (mode == "rate") {
    edge.cols <- rateShifts(PP, threshold, gradientcols, colour)
  }

  plotPhylo(tree, tips = tips, edge.col = edge.cols, scale = scalebar, ...)
  
  if (scalar != "rate") {
    if (scalar == "branch") {
      edgelabels(edge = node_info$nodes, bg = node_info$col, 
        pch = shp, cex = node_info$nodecex)
    } else if (scalar == "nodebranch") {
      nodelabels(node = node_info$nodes, bg = node_info$col, 
        pch = shp, cex = node_info$nodecex)
      edgelabels(edge = branch_info$nodes, bg = branch_info$col, 
        pch = shp, cex = branch_info$nodecex)
    } else {
      nodelabels(node = node_info$nodes, bg = node_info$col, 
        pch = shp, cex = node_info$nodecex)
    }
  }
}


# I think that actually the way to go here is to use the RJPP output to generate
# the info for the plot. THEN have a class for that object, and then have a 
# plot method for it.

# this would also solve the problem of the RJPP table being absurd - there would
# be a "summarise RJPP" type of function that reduces the table down to the
# specific interests of the user, which also returns the edge labels, colours 
# and transparencies that will be needed to plot that object.