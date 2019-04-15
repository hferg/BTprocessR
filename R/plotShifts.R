##############################################################################
#' getPalette
#' Generate the palette for colour ramps.
#' @name getPalette
#' @keywords internal
getPalette <- function(opts) {
  if (opts$palette == "viridis") {
    pal <- viridis::viridis(9)
  } else if (opts$palette == "plasma") {
    pal <- viridis::plasma(9)
  } else if (opts$palette == "magma") {
    pal <- viridis::magma(9)
  } else if (opts$palette == "cividis") {
    pal <- viridis::cividis(9)
  } else {
    pal <- tryCatch({
      gplots::col2hex(opts$palette)
    },
    error = function(e) {
      message("Palette contains invalid colour names:")
      message(e)
    })
  }
  return(pal)
}

##############################################################################
#' generateTrans
#' Generate the transparencies from data and colours
#' @name getPalette
#' @keywords internal
generateTrans <- function(cols, values) {
  tt <- values / max(values)
  cols <- sapply(seq_along(cols), function(x) {
    makeTrans(cols[x], alpha = tt[x])
  })
  return(cols)
}

##############################################################################
#' branchColours
#' Generates the edge colours to colour edges by total rate.
#' @name branchColours
#' @keywords internal
branchColours <- function(PP, opts) {
  # generate data for colour scale
  if (opts$branch.colour == "mean") {
    xx <- log(PP$scalars$meanRate[2:nrow(PP$scalars)])
  } else if (opts$branch.colour == "median") {
    xx <- log(PP$scalars$medianRate[2:nrow(PP$scalars)])
  } else if (opts$branch.colour == "mode") {
    xx <- log(PP$scalars$modeRate[2:nrow(PP$scalars)])
  } else if (opts$branch.colour == "sd") {
    xx <- PP$scalars$sdRate[2:nrow(PP$scalars)]
  } else if (opts$branch.colour == "scale_pc") {
    xx <- PP$scalars$pRate[2:nrow(PP$scalars)]
  }

  # define palette
  pal <- getPalette(opts)

  # now generate edge.colours, if needed.
  if (opts$branch.colour != "none") {
    edge.cols <- plotrix::color.scale(xx, extremes = pal, na.color = NA)
    ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
    scale.cols <- plotrix::color.scale(ss, extremes = pal, na.color = NA)
  } else if (opts$branch.colour == "none") {
    edge.cols <- rep("black", length(xx))
    scale.cols <- NA
  }

  # apply transparency, if needed.
  if (opts$branch.transparency == "scale_pc") {
    edge.cols <- generateTrans(edge.cols, PP$scalars$pRate[2:nrow(PP$scalars)])
  } else if (opts$branch.transparency == "mean") {
    edge.cols <- generateTrans(edge.cols, 
      log(PP$scalars$meanRate[2:nrow(PP$scalars)]))
  } else if (opts$branch.transparency == "median") {
    edge.cols <- generateTrans(edge.cols, 
      log(PP$scalars$medianRate[2:nrow(PP$scalars)]))
  } else if (opts$branch.transparency == "mode") {
    edge.cols <- generateTrans(edge.cols, 
      log(PP$scalars$modeRate[2:nrow(PP$scalars)]))
  } else if (opts$branch.transparency == "sd") {
    edge.cols <- generateTrans(edge.cols, PP$scalars$sdRate[2:nrow(PP$scalars)])
  }

  # identify un-coloured nodes if needed.
  if (opts$coloured.branches == "threshold") {
    tree <- PP$tree_summary$tree_summaries$original_tree
    percs <- PP$scalars$pRate[2:nrow(PP$scalars)]
    names(percs) <- PP$scalars$branch[2:nrow(PP$scalars)]
    nodes <- as.numeric(names(percs[percs < opts$threshold]))
    null_edges <- tree$edge[, 2] %in% nodes
  }

  # turn off colours if coloured.branches == "threshold"
  # i.e. set the na colour.
  if (opts$branch.col != "none" && opts$coloured.branches == "threshold") {
    edge.cols[null_edges] <- opts$na.colour
  }
  return(list(edge.cols = edge.cols, scale.cols = scale.cols))
}

##############################################################################
#' nodeShapes
#' Generates the node labels, edge.colours and transparencies to plot the
#' location of transformations in a posterior.
#' @name nodeShapes
#' @keywords internal

nodeShapes <- function(PP, opts, mode) {
  
  # old arguments - threshold, cl, tree, transparency, relativetrans,
  # nodescaling, colour, nodecex

  # if threshold is zero, make the threshold the minimum threshold - this means
  # that nodes that never get a scalar don't get a little circle plotted on
  # them!
  if (opts$threshold == 0) {
    threshold <- 1 / PP$niter
  }

  # get correct information
  if (opts$transformation == "rates") {
    if (mode == "branches") {
      brates <- PP$origins$branches[ , 3:ncol(PP$origins$branches)]
      pS <- PP$scalars$nOrgnBRate / PP$niter
      pVmean <- rowMeans(brates)
      pVmedian <- apply(brates, 1, median)
      pVmode <- apply(brates, 1, modeStat)
      pVsd <- apply(brates, 1, sd)
    } else if (mode == "nodes") {
      nrates <- PP$origins$nodes[ , 3:ncol(PP$origins$nodes)]
      pS <- PP$scalars$nOrgnNRate / PP$niter
      pVmean <- rowMeans(nrates)
      pVmedian <- apply(nrates, 1, median)
      pVmode <- apply(nrates, 1, modeStat)
      pVsd <- apply(nrates, 1, sd)
    }
  } else if (opts$transformation == "delta") {
    pS <- PP$scalars$pDelta
    pVmean <- log(PP$scalars$meanDelta)
    pVmedian <- log(PP$scalars$medianDelta)
    pVmode <- log(PP$scalars$modeDelta)
    pVsd <- PP$scalars$sdDelta
  } else if (opts$transformation == "lambda") {
    pS <- PP$scalars$pLambda
    pVmean <- log(PP$scalars$meanLambda)
    pVmedian <- log(PP$scalars$medianLambda)
    pVmode <- log(PP$scalars$modeLambda)
    pVsd <- PP$scalars$sdDelta
  } else if (opts$transformation == "kappa") {
    pS <- PP$scalars$pKappa
    pVmean <- log(PP$scalars$meanKappa)
    pVmedian <- log(PP$scalars$medianKappa)
    pVmode <- log(PP$scalars$modeKappa)
    pVsd <- PP$scalars$sdKappa
  }

  # find nodes that get labels.
  node_tf <- pS > threshold
  nodes <- PP$scalars$descNode[node_tf]

  # come back to this... what is it asking?!

  if (length(nodes) == 0) {
    stop("No scalars above threshold.")
  }

  if (opts$nb.colour != "none") {
    pal <- getPalette(opts)
    makeCols <- function(xx, pal) {
      cols <- plotrix::color.scale(xx, extremes = pal, 
        na.color = NA)
      ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
      scale.cols <- plotrix::color.scale(ss, extremes = pal, na.color = NA)
      return(list(cols = cols, scale.cols = scale.cols))
    }
    # define the colours
    if (opts$nb.colour == "mean") {
      bcols <- makeCols(pVmean[node_tf], pal)
    } else if (opts$nb.colour == "median") {
      bcols <- makeCols(pVmedian[node_tf], pal)
    } else if (opts$nb.colour == "mode") {
      bcols <- makeCols(pVmode[node_tf], pal)
    } else if (opts$nb.colour == "sd") {
      bcols <- makeCols(pVsd[node_tf], pal)
    } else if (opts$nb.colour == "scale_pc") {
      bcols <- makeCols(pS[node_tf], pal)
    }
  } else {
    cols <- nb.fill
  }

  # define transparencies if needed.
  if (opts$shape.transparency == "scale_pc") {
    trans <- generateTrans(cols$cols, pS)
  } else if (opts$shape.transparency == "mean") {
    trans <- generateTrans(cols$cols, pVmean)
  } else if (opts$shape.transparency == "median") {
    trans <- generateTrans(cols$cols, pVmedian)
  } else if (opts$shape.transparency == "mode") {
    trans <- generateTrans(cols$cols, pVmode)
  } else if (opts$shape.transparency == "sd") {
    trans <- generateTrans(cols$cols, pVsd)
  }

  # define scale - then match to the nodes
  if (opts$nb.scale == "mean") {
    pcex <- opts$nb.cex * pVmean[node_tf]
  } else if (opts$nb.scale == "median") {
    pcex <- opts$nb.cex * pVmedian[node_tf]
  } else if (opts$nb.scale == "mode") {
    pcex <- opts$nb.cex * pVmode[node_tf]
  } else if (opts$nb.scale == "sd") {
    pcex <- opts$nb.cex * pVsd[node_tf]
  } else if (opts$nb.scale == "scale_pc") {
    pcex <- opts$nb.cex * pS[node_tf]
  }

# OLD

  pprobs <- PP$scalars[which(
    (PP$scalars[ , cl] / PP$niter) >= threshold
  ) , cl] / PP$niter

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
    # I think I need to rescale nodecex so that the max is 1? That way there
    # will always be something plotted...
    nodecex = nodecex * pprobs
  }

  if (cl == "nOrgnBRate") {
    nodes <- which(tree$edge[ , 2] %in% nodes)
  }

  list(nodes = as.numeric(nodes), 
    colours = cols$cols, 
    alphas = trans,
    scale.cols = cols$scale.cols, 
    nodecex = as.numeric(nodecex[, 1]))
}

################################################################################
#' plotShifts
#'
#' Plots the locations of the origins of scalars from the postprocessor output
#' of bayestraits.
#' CURRENTLY WORKS ONLY FOR DELTAS.
#' @param PP The psotprocessor (localscalrPP) output.
#' @param scalar The scalar to find and plot from the post processor -
#' delta/lambda/kappa/node/branch
#' @param threshold Threshold of probability in posterior to display deltas for,
#' defaults to zero (i.e. shows all deltas shaded proportionally to the
#' posterior probability)
#' @param colour The colour to use for the node circles
#' @param scaled Plot the original tree (scaled = "time", the default), or the
#' mean/sclaed tree (scaled = "mean") or plot the tree scaled only by scalars
#' present above the threshold (scaled = "threshold")?
#' @param nodecex The scaling factor for the size of the node circles
#' @param tips Show tip labels?
#' @param scalebar Include scale bar?
#' @param measure When plotting "siginficant" tree, what measure of the
#' parameter? Median (default), mode or mean.
#' @param exludeones If plotting according to a threshold of significance,
#' should 1s (i.e. no scalar) be excluded from the posterior when calculating
#' average scalar?
#' @param relativetrans If TRUE (defaults to FALSE) the scale of transparency
#' will go from the threshold (totally transparent) to the maximum presence
#' (full opacity).
#' @param nodescaling Scale node symbols according to posterior probability of
#' shift (default).
#' @param transparency Adjust node symbol transparency according to posterior
#' probability? Defaults to FALSE.
#' @param gradientcols A vector of two colours - the min and max colours used
#' when colouring the tree according to percentage time rate scaled (when
#' threshold = 0) or using rate.edges.
#' @param rate.edges Takes a numeric value between 0 and 1. If NULL (default)
#' then node shapes are plotted for a transformation. If equal to zero then node
#' shapes are plotted along with a colour gradient on the branches for rates (if
#' in the posterior), and if set to a threshold then the branches are coloured
#' black/red for whether there is a scalar over the threshold (red) along with
#' node scalars.
#' @param shp The shape of the node markers (uses the usual pch index).
#' @name plotShifts
#' @import plotrix
#' @export

# I think there might be too many options here - what is actually useful?!
# rework notes.
# RATE PLOTTING.
#   Threshold determines which branches have a colour applied to them. 
#   Then the coloured branches are either coloured by the percentage of time
#   they are scaled, OR by the magnitude of the scalar.
#   Third option is to plot the tree scaled according to the magnitude of the
#     scalar and then have the colour of the branch correspond to the percentage
#     of time the branch is receiving a scalar.
#   DEFAULT BEHAVIOUR.
#     Plot branches coloured by rate, and branches below the threshold remain
#     uncoloured. Scale bar shows rate. Then pass through a control list (or
#     something?!) to allow more complex mixing of parameters...

plotShifts <- function(PP, scalar, threshold = 0, threshold2 = 0, nodecex = 2,
  scaled = "time", scalebar = TRUE, measure = "median", excludeones = FALSE,
  relativetrans = FALSE, nodescaling = TRUE, transparency = FALSE,
  gradientcols = c("dodgerblue", "firebrick1"), rate.edges = NULL,
  colour = "red", colour2 = "green",  shp = 21, tips = FALSE, ...) {

  # options.
  #   threshold = numeric, 0-1
  #   transformation = rates, delta, kappa, lambda
  #   branch.colour = none, mean, median, mode, sd, scale_pc
  #   branch.transparency = none, scale_pc, mean, median, mode, sd
  #   coloured.brances = threshold, all, none
  #   nb.colour = mean, median, mode, sd, scale_pc, none
  #   nb.transparency = none, mean, median, mode, sd, scale_pc
  #   nb.fill = any single colour
  #   nb.border = any single colour
  #   nb.shape = a pch number, OR directional
  #   nb.cex = scaling factor for the node/branch shapes
  #   branch.scale = none, mean, median, mode, sd, scale_pc
  #   scaled.branches = threshold, all, none
  #   na.colour = any colour (name or hex)
  #   palette = viridis, plasma, magma, inferno, cividis, at least two colours
  #   rates.mode = all, rates, nodes, branches, nodes_branches

  transformation <- names(PP$origins)
  if ("nodes" %in% transformation) {
    transformation <- c("rates", 
      transformation[!transformation %in% c("nodes", "branches")])
  }

  opts <- list(
    threshold = 0,
    transformation = transformation,
    branch.colour = "mean",
    branch.transparency = "none",
    coloured.branches = "threshold",
    nb.colour = "mean",
    nb.transparency = "none",
    nb.fill = "red",
    nb.border = "black",
    nb.scale = "mean",
    nb.shape = 21,
    nb.cex = 2,
    branch.scale = "none",
    scaled.branches = "threshold",
    na.colour = "black",
    palette = "viridis",
    rates.mode = "all"
  )
  opts[names(plot.options)] <- plot.options

  # if there are multiple transformations, ask the user which to return.
  if (length(opts$transformation) >1 ) {
    message("Multiple transformations detected. Please select one to plot:")
    opts$transformation <- select.list(choices = opts$transformation)
  }

  if (opts$transformation == "rates") {
    if (opts$rates.mode == "all" | opts$rates.mode == "rates") {
      # get edge colours for linear rates.
      branch_cols <- branchColours(PP, opts)
    }

    if (opts$rate.mode == "nodes" | rates$rates.mode == "nodes_branches" | 
      opts$rates.mode == "all") {
      # get node shapes
      node_shapes <- nodeShapes(PP, opts, mode = "nodes")
    }

    if (opts$rates.mode == "branches" | opts$rates.mode == "nodes_branches" | 
      opts$rates.mode == "all") {
      # get branch shapes
      branch_shapes <- nodeShapes(PP, opts, mode = "branches")
    }

  } else if (opts$transformation == "delta" | opts$transformation == "lamba" |
    opts$tranformation == "kappa") {
    # get node shapes
    node_shapes <- nodeShapes(PP, opts, mode = "nodes")
  }

    node_info <- transShifts(
      PP,
      threshold,
      cl,
      tree,
      transparency,
      relativetrans,
      nodescaling,
      colour,
      nodecex
    )

  # here is the plotting - this needs refining - a multi-panel plot for 
  # rates.mode = ALL, and then the correct colour scales etc. to be added
  # as well.
  plotPhylo(tree, tips = tips, edge.col = edge.cols, scale = scalebar, ...)

  if (scalar != "rate") {
    if (scalar == "branch") {
      ape::edgelabels(edge = node_info$nodes, bg = node_info$col,
        pch = shp, cex = node_info$nodecex)
    } else if (scalar == "nodebranch") {
      ape::nodelabels(node = node_info$nodes, bg = node_info$col,
        pch = shp, cex = node_info$nodecex)
      ape::edgelabels(edge = branch_info$nodes, bg = branch_info$col,
        pch = shp, cex = branch_info$nodecex)
    } else {
      ape::nodelabels(node = node_info$nodes, bg = node_info$col,
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

# general plotShifts workflow...

# 1) identify the transformation(s) requested
# 2) generate edge colours for each of those transformations (see control list)
# 3) generate node labels for each of those transformations (see control list)
# 4) transform the tree for each of those transformations (see control list)
# 5) plot each of the trees with the right colours, node labels and branch
#   lengths according to specified control list.
# 6) add the scale bar for the colour gradient (if present).
# 7) add a time axis?!

# The plots that are generated for the presence of each transformation.
# RATES
#   1) A tree with branches coloured according to rate - not scaled.
#   2) A tree with node scalars over the threshold plotted, coloured according
#     to scalar value
#   3) A tree with branch scalars over the threshold plotted, coloured according
#     to scalar value