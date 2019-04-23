##############################################################################
#' generateTrans
#' Generate the transparencies from data and colours
#' @param cols The colours to make transparent
#' @param values The data to set transparency by
#' @name generateTrans
#' @keywords internal

generateTrans <- function(cols, values) {
  tt <- values / max(values)
  cols <- sapply(seq_along(cols), function(x) {
    makeTrans(cols[x], alpha = tt[x])
  })
  return(cols)
}

##############################################################################
#' funLegend
#' Generate some interesting legend labels.
#' @name funLegend
#' @keywords internal

funLegend <- function() {
  names <- list(
    c("Jimmy Cliff", "Jimi Hendrix"),
    c("Neil Sedaka", "Neil Young"),
    c("James Blunt", "James Chance"),
    c("David Gray", "David Byrne"),
    c("John Mayer", "John Bonham")
  )
  sample(names, 1)[[1]]
}

##############################################################################
#' branchColours
#' Generates the edge colours to colour edges by total rate.
#' @name branchColours
#' @keywords internal
branchColours <- function(PP, opts) {
  plot.dat <- list(
    scale_pc = PP$scalars$pRate[2:nrow(PP$scalars)],
    mean = log(PP$scalars$meanRate[2:nrow(PP$scalars)]),
    median = log(PP$scalars$medianRate[2:nrow(PP$scalars)]),
    mode = PP$scalars$modeRate[2:nrow(PP$scalars)],
    sd = PP$scalars$sdRate[2:nrow(PP$scalars)]
  )

  if (opts$edge.colour != "none") {
    xx <- plot.dat[[opts$edge.colour]]
    edge.cols <- colourvalues::colour_values(xx, palette = opts$edge.palette)
    ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
    scale.cols <- colourvalues::colour_values(ss, palette = opts$edge.palette)
    scale.lims <- range(ss)
  } else if (opts$edge.colour == "none") {
    edge.cols <- rep("black", length(xx))
    scale.cols <- NA
  }

  if (opts$edge.transparency != "none") {
    edge.cols <- generateTrans(edge.cols, plot.dat[[opts$edge.transparency]])
  }

  if (opts$coloured.edges == "threshold") {
    tree <- PP$tree_summary$tree_summaries$original_tree
    names(plot.dat$scale_pc) <- PP$scalars$branch[2:nrow(PP$scalars)]
    nodes <- as.numeric(
      names(plot.dat$scale_pc[plot.dat$scale_pc < opts$threshold])
    )
    null_edges <- tree$edge[, 2] %in% nodes
  }

  if (opts$edge.colour != "none" && opts$coloured.edges == "threshold") {
    edge.cols[null_edges] <- opts$na.colour
  }

  return(list(edge.cols = edge.cols, scale.cols = scale.cols, 
    scale.lims = scale.lims))
}

##############################################################################
#' shapeCols
#' Generates colours for shapes for either node or branch scalars.
#' @name shapeCols
#' @keywords internal

shapeCols <- function(opts, plot.dat, mode) {
  if (mode == "nodes") {
    if (opts$node.colour == "none") {
      cols <- scale.cols <- opts$node.fill
      scale.lims <- c(0, 0)
    } else {
      xx <- plot.dat[[opts$node.colour]]
      cols <- colourvalues::colour_values(xx, palette = opts$node.palette)
      ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
      scale.cols <- colourvalues::colour_values(ss, palette = opts$node.palette)
      scale.lims <- range(ss)
    }
  } else if (mode == "branches") {
    if (opts$branch.colour == "none") {
      cols <- scale.cols <- opts$branch.fill
      scale.lims <- c(0, 0)
    } else {
      xx <- plot.dat[[opts$branch.colour]]
      cols <- colourvalues::colour_values(xx, palette = opts$branch.palette)
      ss <- seq.int(from = min(xx), to = max(xx), length.out = length(xx))
      scale.cols <- colourvalues::colour_values(ss, 
        palette = opts$branch.palette)
      scale.lims <- range(ss)
    }
  }
  return(list(cols = cols, scale.cols = scale.cols, scale.lims = scale.lims))
}

##############################################################################
#' plotShapes
#' Generates the node labels, edge.colours and transparencies to plot the
#' location of transformations in a posterior.
#' @name plotShapes
#' @keywords internal

plotShapes <- function(PP, opts, mode) {
  
  if (opts$threshold == 0) {
    opts$threshold <- 1 / PP$niter
  }

  if (opts$transformation == "rate") {
    if (mode == "branches") {
      brates <- PP$origins$branches[-1 , 3:ncol(PP$origins$branches)]
      plot.dat <- list(
        scale_pc = PP$scalars$nOrgnBRate[-1] / PP$niter,
        mean = log(rowMeans(brates)),
        median = log(apply(brates, 1, median)),
        mode = apply(brates, 1, modeStat),
        sd = apply(brates, 1, sd)
      )
    } else if (mode == "nodes") {
      nrates <- PP$origins$nodes[-1 , 3:ncol(PP$origins$nodes)]
      plot.dat <- list(
        scale_pc = PP$scalars$nOrgnNRate[-1] / PP$niter,
        mean = log(rowMeans(nrates)),
        median = log(apply(nrates, 1, median)),
        mode = apply(nrates, 1, modeStat),
        sd = apply(nrates, 1, sd)
      )
    }
  } else if (opts$transformation == "delta") {
    plot.dat <- list(
      scale_pc = PP$scalars$pDelta,
      mean = log(PP$scalars$meanDelta),
      median = log(PP$scalars$medianDelta),
      mode = PP$scalars$modeDelta,
      sd = PP$scalars$sdDelta
    )
  } else if (opts$transformation == "lambda") {
    plot.dat <- list(
      scale_pc = PP$scalars$pLambda,
      mean = log(PP$scalars$meanLambda),
      median = log(PP$scalars$medianLambda),
      mode = PP$scalars$modeLambda,
      sd = PP$scalars$sdDelta
    )
  } else if (opts$transformation == "kappa") {
    plot.dat <- list(
      scale_pc = PP$scalars$pKappa,
      mean = log(PP$scalars$meanKappa),
      median = log(PP$scalars$medianKappa),
      mode = PP$scalars$modeKappa,
      sd = PP$scalars$sdKappa
    )
  }
  
  node_tf <- plot.dat$scale_pc > opts$threshold
  if (mode == "nodes") {
    nodes <- PP$scalars$descNode[-1][node_tf]
  } else if (mode == "branches") {
    nodes <- PP$scalars$branch[-1][node_tf]
  }
  
  # this needs to be changed to return something useful.
  if (length(nodes) == 0) {
    ret <- list(nodes = NULL,
      colours = NULL,
      scale.cols = opts$na.colour,
      scale.lims = c(0, 1))
  } else {
    plot.dat <- lapply(plot.dat, function(x) x[node_tf])
    plot.cols <- shapeCols(opts, plot.dat, mode)

    if (mode == "nodes" && opts$node.transparency != "none") {
      plot.cols$cols <- generateTrans(plot.cols$cols, 
        plot.dat[[opts$node.transparency]])
    } else if (mode == "branches" && opts$branch.transparency != "none") {
      plot.cols$cols <- generateTrans(plot.cols$cols, 
        plot.dat[[opts$branch.transparency]])
    }

    if (mode == "nodes") {
      if (opts$node.scale != "none") {
        pcex <- opts$node.cex * plot.dat[[opts$node.scale]]
      } else {
        pcex <- opts$node.cex
      }
    } else if (mode == "branches") {
      if (opts$branch.scale != "none") {
        pcex <- opts$node.cex * plot.dat[[opts$branch.scale]]
      } else {
        pcex <- opts$branch.cex
      }
    }
    ret <- list(nodes = as.numeric(nodes), 
      colours = plot.cols$cols, 
      scale.cols = plot.cols$scale.cols,
      scale.lims = plot.cols$scale.lims,
      nodecex = pcex)
  }
  return(ret)
}

##############################################################################
#' legendInfo
#' Generates legend info for colour scale legends
#' @name legendInfo
#' @keywords internal

legendInfo <- function(tree, opts, cols) {
  if (opts$legend != "numeric") {
    leg <- funLegend()
  } else {
    leg <- round(cols$scale.lims, 2)
  }
  if (opts$legend.pos == "auto") {
    pos <- c(0, 
      0, 
      round((max(ape::node.depth.edgelength(tree)) / 4), 2),
      round((length(tree$tip.label) / 30), 2))
  } else {
    pos <- opts$legend.pos
  }
  return(list(legend = leg, pos = pos))
}

##############################################################################
#' makeLegendTitle
#' Generates legend title for scale bars.
#' @name makeLegendTitle
#' @keywords internal

makeLegendTitle <- function(opts) {
  if (opts$edge.colour == "mean" | opts$edge.colour == "median") {
    edge_leg <- Hmisc::capitalize(paste("ln", opts$edge.colour, 
      opts$transformation))
  } else {
    edge_leg <- Hmisc::capitalize(paste(opts$edge.colour, 
      opts$transformation))
  }
  if (opts$node.colour == "mean" | opts$node.colour == "median") {
    node_leg <- Hmisc::capitalize(paste("ln", opts$node.colour, 
      opts$transformation))
  } else {
    node_leg <- Hmisc::capitalize(paste(opts$node.colour, 
      opts$transformation))
  }
  if (opts$branch.colour == "mean" | opts$branch.colour == "median") {
    branch_leg <- Hmisc::capitalize(paste("ln", opts$branch.colour, 
      opts$transformation))
  } else {
    branch_leg <- Hmisc::capitalize(paste(opts$branch.colour, 
      opts$transformation))
  }
  return(list(edge_leg = edge_leg, 
    node_leg = node_leg, 
    branch_leg = branch_leg))
}

##############################################################################
#' scaleTree
#' Scales the tree for plotting
#' @name scaleTree
#' @keywords internal

scaleTree <- function(PP, opts) {
  if (opts$edge.scale == "none") {
    tree <- PP$tree_summary$tree_summaries$original_tree
  } else if (opts$edge.scale == "mean") {
    tree <- PP$tree_summary$tree_summaries$mean_tree
  } else if (opts$edge.scale == "median") {
    tree <- PP$tree_summary$tree_summaries$median_tree
  } else if (opts$edge.scale == "mode") {
    tree <- PP$tree_summary$tree_summaries$mode_tree
  }

  if (opts$edge.scale != "none" && opts$scaled.edges == "threshold") {
    psc <- PP$scalars[-1, ]
    sc <- psc$pScalar < opts$threshold
    tree$edge.length[sc] <- 
      PP$tree_summary$tree_summaries$original_tree$edge.length[sc]
  }
  return(tree)
}

################################################################################
#' plotShifts
#'
#' Plots the locations of the origins of scalars from the postprocessor output
#' of bayestraits.
#' @param PP The output of the rjpp function.
#' @param plot.options A list of control options. See description.
#' @details The default behaviour of plotShifts depends on the transformations
#' present in the rjpp output. If variable rates, then 3 trees will be plotted:
#' the first has branches coloured according the log of the mean rate, the 
#' second shows all node scalars present more than once in the posterior, 
#' coloured according to the mean log rate and the third shows the same for 
#' branch scalars. If delta, kappa or lambda are present then a single tree is 
#' plotted showing all nodes that receive a scalar, coloured according to mean 
#' magnitude. If multiple transformations are present then the user will be
#' prompted to select one. 
#' The plot.options list provides a high degree of control over what
#' is plotted, allowing the default behaviour to be customised. The options, and
#' values that they can take, are as follows.
#' 
#' \itemize{
#' \item{"threshold"}{0-1 - The threshold of presence in the posterior over which 
#' a node and/or branch scalar is plotted. Also the threshold referenced by 
#' coloured.edges and scaled.edges.}
#' \item{"transformation"}{rate, delta, lambda, kappa - The transformation to 
#' plot.}
#' \item{"edge.colour"}{none, mean, median, mode, sd, scale_pc - the metric to
#' colour edges by. If none branches default to the na.colour option. Mean,
#' median, mode and sd correspond to the appropriate branch lengths from the 
#' posterior of trees and scale_pc colours edges according to the percentage of
#' time they are scaled in the posterior.}
#' \item{"edge.transparency"}{none, scale_pc, sd - the measure to make edges 
#' proportionally transparent by. None results in uniform solid branches, 
#' scale_pc gives edges that are scaled less frequently in the posterior higher
#' transparency, and sd gives branches that have higher SD of estimated branch
#' lengths more solid colours.}
#' }
#' 
#' 
#' @name plotShifts
#' @import plotrix
#' @export

plotShifts <- function(PP, plot.options = list(), ...) {

  # options.
  #   threshold = numeric, 0-1
  #   transformation = rates, delta, kappa, lambda
  #   branch.colour = none, mean, median, mode, sd, scale_pc
  #   branch.transparency = none, scale_pc, sd
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
  #   rates.mode = all_single, all_panel, rates, nodes, branches, nodes_branches
  #   legend.pos = auto, OR a vectro with coords for bottom left and top right of legend
  #   legen = numeric

  # TODO. Work out what to do when scaling factors are negative (e.g. from 
  # logged data). Just add 1? I think first finish the function and then see 
  # what comes out with different options, see if it makes sense etc.
  # For SOME reason the colour scaling for the nodes and branches has gone out
  # of whack again - fix.

  transformation <- names(PP$origins)
  if ("nodes" %in% transformation) {
    transformation <- c("rate", 
      transformation[!transformation %in% c("nodes", "branches")])
  }

  opts <- list(
    threshold = 0,
    transformation = transformation,
    edge.colour = "mean",
    edge.transparency = "none",
    coloured.edges = "threshold",
    edge.palette = "viridis",
    edge.scale = "none",
    scaled.edges = "threshold",
    node.colour = "mean",
    node.scale = "none",
    node.transparency = "none",
    node.palette = "viridis",
    node.fill = "red",
    node.border = "black",
    node.shape = 21,
    node.cex = 2,
    branch.colour = "mean",
    branch.transparency = "none",
    branch.palette = "viridis",
    branch.fill = "red",
    branch.border = "black",
    branch.scale = "none",
    branch.shape = 22,
    branch.cex = 2,
    na.colour = "black",
    rates.mode = "all_panel",
    legend.pos = "auto",
    legend = "numeric"
  )
  opts[names(plot.options)] <- plot.options

  if (length(opts$transformation) >1 ) {
    message("Multiple transformations detected. Please select one to plot:")
    opts$transformation <- select.list(choices = opts$transformation)
  }

  if (opts$transformation == "rate") {
    if (opts$rates.mode == "all_panel" | 
        opts$rates.mode == "all_single" | 
        opts$rates.mode == "rates") {
      edge.cols <- branchColours(PP, opts)
    }

    if (opts$rates.mode == "nodes" | opts$rates.mode == "nodes_branches" | 
      opts$rates.mode == "all_single" | opts$rates.mode == "all_panel") {
      node.shapes <- plotShapes(PP, opts, mode = "nodes")
    }

    if (opts$rates.mode == "branches" | opts$rates.mode == "nodes_branches" | 
      opts$rates.mode == "all_panel" | opts$rates.mode == "all_single") {
      branch.shapes <- plotShapes(PP, opts, mode = "branches")
    }

  } else if (opts$transformation == "delta" | opts$transformation == "lamba" |
    opts$tranformation == "kappa") {
    node.shapes <- plotShapes(PP, opts, mode = "nodes")
  }

  tree <- scaleTree(PP, opts)
  leg.tit <- makeLegendTitle(opts)

  if (opts$transformation == "rate" && opts$rates.mode == "all_panel") {
    par(mfrow = c(1, 3))
    plotPhylo(tree, edge.col = edge.cols$edge.cols)
      leg <- legendInfo(tree, opts, edge.cols)
      plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
        leg$legend, rect.col = edge.cols$scale.cols, align = "rb")
      labx <- (leg$pos[1] + leg$pos[3]) / 2
      laby <- leg$pos[4] + strheight("M")
      text(labx, laby, leg.tit$edge_leg)
    plotPhylo(tree)
      ape::nodelabels(node = node.shapes$nodes, 
        bg = node.shapes$colours,
        col = opts$node.border,
        cex = node.shapes$nodecex, 
        pch = opts$node.shape)
      leg <- legendInfo(tree, opts, node.shapes)
      plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
            leg$legend, rect.col = node.shapes$scale.cols, align = "rb")
      labx <- (leg$pos[1] + leg$pos[3]) / 2
      laby <- leg$pos[4] + strheight("M")
      text(labx, laby, leg.tit$node_leg)
    plotPhylo(tree)
      ape::edgelabels(edge = branch.shapes$nodes, 
        bg = branch.shapes$colours,
        col = opts$branch.border, 
        cex = branch.shapes$nodecex, 
        pch = opts$branch.shape)
      leg <- legendInfo(tree, opts, branch.shapes)
      plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
        leg$legend, rect.col = branch.shapes$scale.cols, align = "rb")
      labx <- (leg$pos[1] + leg$pos[3]) / 2
      laby <- leg$pos[4] + strheight("M")
      text(labx, laby, leg.tit$branch_leg)
  } else if (opts$transformation == "rate" && opts$rates.mode == "all_single") {
    # Plot all on single plot
  } else if (opts$transformation == "rate" && opts$rates.mode == "nodes_brancehs") {
    # plot node and branch labels together
  } else if (opts$transformation == "rate" && opts$rates.mode == "nodes") {
    plotPhylo(tree)
    ape::nodelabels(node = node.shapes$nodes, 
      bg = node.shapes$colours,
      col = opts$node.border,
      cex = node.shapes$nodecex, 
      pch = opts$node.shape)
    leg <- legendInfo(tree, opts, node.shapes)    
    plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
          leg$legend, rect.col = node.shapes$scale.cols, align = "rb")
    labx <- (leg$pos[1] + leg$pos[3]) / 2
    laby <- leg$pos[4] + strheight("M")
    text(labx, laby, leg.tit$node_leg)
  } else if (opts$transformation == "rate" && opts$rates.mode == "branches"){
    plotPhylo(tree)
    ape::edgelabels(edge = branch.shapes$nodes, 
      bg = branch.shapes$colours,
      col = opts$branch.border, 
      cex = branch.shapes$nodecex, 
      pch = opts$branch.shape)
    leg <- legendInfo(tree, opts, branch.shapes)
    plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
      leg$legend, rect.col = branch.shapes$scale.cols, align = "rb")
    labx <- (leg$pos[1] + leg$pos[3]) / 2
    laby <- leg$pos[4] + strheight("M")
    text(labx, laby, leg.tit$branch_leg)
  } else if (opts$transformation == "rate" && opts$rates.mode == "rates") {
    plotPhylo(tree, edge.col = edge.cols$edge.cols)
    leg <- legendInfo(tree, opts, edge.cols)
    plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
      leg$legend, rect.col = edge.cols$scale.cols, align = "rb")
    labx <- (leg$pos[1] + leg$pos[3]) / 2
    laby <- leg$pos[4] + strheight("M")
    text(labx, laby, leg.tit$edge_leg)
  } else if (opts$transformation == "delta" | 
             opts$transformation == "kappa" |
             opts$transformation == "lambda") {
    plotPhylo(tree)
    ape::nodelabels(node = node.shapes$nodes, 
      bg = node.shapes$colours,
      col = opts$node.border,
      cex = node.shapes$nodecex, 
      pch = opts$node.shape)
    leg <- legendInfo(tree, opts, node.shapes)
    plotrix::color.legend(leg$pos[1], leg$pos[2], leg$pos[3], leg$pos[4], 
      leg$legend, rect.col = node.shapes$scale.cols, align = "rb")
    labx <- (leg$pos[1] + leg$pos[3]) / 2
    laby <- leg$pos[4] + strheight("M")
    text(labx, laby, leg.tit$node_leg)
  }
}
