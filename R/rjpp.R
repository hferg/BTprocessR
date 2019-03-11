##############################################################################
#' loadRJ
#'
#' Returns the full mcmc object from a BayesTraits log file. This
#' is used inside plot functions and so on, but might be useful for
#' other MCMC manipulations and so on.
#' @param logfile The name of the logfile of the BayesTraits analysis.
#' @return A list containing the taxa translation table, all possible subtrees 
#' a scalar can occur on, and a data frame of the rj model configuration.

loadRJ <- function(logfile, burnin = 0, thinning = 1) {

  raw <- readLines(logfile)
  rawhead <- strsplit(raw[1:(grep("\\bIt*\\b", raw) -1)], "\t")
  rawtail <- strsplit(raw[grep("\\bIt*\\b", raw):length(raw)], "\t")
  nms1 <- rawtail[[1]][1:7]
  nms2 <- rawtail[[1]][8:length(rawtail[[1]])]
  nms2 <- gsub(" ", "", nms2)
  nms2 <- gsub("/", "", nms2)  

  for (i in 1:length(rawtail)) {
    if (length(rawtail[[i]]) == 7) {
      names(rawtail[[i]]) <- nms1
    } else {
      len <- length(rawtail[[i]][8:length(rawtail[[i]])])
      end <- vector(mode = "character", length = len)
      
      st <- 1
      ed <- 4
      for (j in 1:(len/4)) { 
        end[c(st:ed)] <- paste(nms2, j, sep = "_")
        st <- st + 4
        ed <- ed + 4
      }
      nms <- c(nms1, end)
      names(rawtail[[i]]) <- nms
    }
  }

  tipnum <- rawhead[[1]]
  taxatrans <- do.call(rbind, rawhead[c(1:tipnum+1)])
  subtreestart <- nrow(taxatrans) + 3
  subtrees <- rawhead[subtreestart:length(rawhead)]
  
  for (i in 1:length(subtrees)) {
    names(subtrees[[i]]) <- c(1:length(subtrees[[i]]))
  }
  
  output <- do.call(smartBind, rawtail)
  output <- output[seq.int(burnin, nrow(output), thinning), ]
  output <- output[2:nrow(output), ]
  subtrees <- do.call(smartBind, subtrees)
  subtrees <- data.frame(subtrees, stringsAsFactors = FALSE)
  colnames(subtrees)[c(1:2)] <- c("node", "bl")

  res <- list(taxatrans, subtrees, output)
  names(res) <- c("taxa", "subtrees", "rj_output")
  return(res)
}


##############################################################################
#' createCounts
#' Creates the counts table for the rjpp.
#' @param extree The time tree
#' @param meanbl The output of meanBranches
#' @name createCounts
#' @keywords internal

#' NOTES
#' This makes a massive and unwieldy table that is populated later on.
#' Take a look at this and see if:
#'    a) the headings are sensible and descriptive
#'    b) can this table be split into more than one table (e.g. one for 
#'      rates information, and one for branchlength information?)
#'    c) If there are possible splits, can the occur here, or is it more
#'      useful later on.
#'    d) Is the table being built in the best way?
#' The species descendent from each tip 100% needs to be a seperate list
#' and it's own element in the results - this makes for ludicrous plotting.
#' It should also have a method that prints the names properly when called,
#' and when a specific node is selected.

createCountsTable <- function(reftree, tree_summary) {
  counts <- matrix(ncol = 53, nrow = (nrow(reftree$edge) + 1))
  colnames(counts) <- c("branch", "ancNode", "descNode", "nTips", "start", 
      "end", "mid", "original_bl", "mean_bl", "median_bl", "mode_bl", 
      "range_bl", "sd_bl", "itersScaled", "itersRatescaled", "itersDelta", 
      "itersKappa", "itersLambda", "pScaled", "pRate", "pDelta", "pKappa", 
      "pLambda", "nScalar", "nRate", "nDelta", "nKappa", "nLambda", 
      "nOrgnScalar", "nOrgnNRate", "nOrgnBRate", "nOrgnDelta", "nOrgnKappa", 
      "nOrgnLambda", "rangeRate", "lqRate", "uqRate", "meanRate", "medianRate", 
      "modeRate", "rangeDelta", "meanDelta", "medianDelta", "modeDelta",
      "rangeKappa", "meanKappa", "medianKappa", "modeKappa",
      "rangeLambda", "meanLambda", "medianLambda", "modeLambda", "species")

  counts[ , "branch"] <- c(0:nrow(reftree$edge))
  counts[ , "ancNode"] <- c(0, reftree$edge[ , 1])
  counts[ , "descNode"] <- c((length(reftree$tip.label) + 1), 
    reftree$edge[ , 2])
  counts[ , "original_bl"] <- c(0, reftree$edge.length)
  
  info <- tree_summary$branchlength_info
  counts[ , "mean_bl"] <- c(0, info$mean_bl)
  counts[ , "median_bl"] <- c(0, info$median_bl)
  counts[ , "mode_bl"] <- c(0, info$mode_bl)
  counts[ , "range_bl"] <- c(0, info$range_bl)
  counts[ , "sd_bl"] <- c(0, info$sd_bl)

  hts <- phytools::nodeHeights(reftree)
  hts <- round(abs(hts - max(hts)), 4)
  counts[ , "start"] <- c(0, hts[ , 1])
  counts[ , "end"] <- c(0, hts[ , 2])
  counts <- as.data.frame(counts)

  # Deal with the root
  counts[1, 1] <- "root"
  descs <- getDescs(reftree, node = counts[1, "descNode"])
  counts[1, "nTips"] <- sum(descs <= length(reftree$tip.label))
  counts[1, "mid"] <- 0
  counts[1, "species"] <- paste0(reftree$tip.label[order(reftree$tip.label)], 
    collapse = ",")

  for (i in 2:nrow(counts)) {
    descs <- getDescs(reftree, node = counts[i, "descNode"])
    counts[i, "nTips"] <- sum(descs <= length(reftree$tip.label))
    if (counts[i, "nTips"] == 0) {
      counts[i, "nTips"] <- 1
    }
    if (counts[i, "descNode"] <= length(reftree$tip.label)) {
      counts[i, "species"] <- reftree$tip.label[counts[i, "descNode"]]
    } else {
      tips <- getDescs(reftree, counts[i, "descNode"])
      tips <- tips[tips <= length(reftree$tip.label)]
      tips <- reftree$tip.label[tips]
      counts[i, "species"] <- paste0(sort(tips), collapse = ",")
    }
    counts[i, "mid"] <- mean(c(hts[(i - 1), 1], hts[(i - 1), 2]))
  }
  
  counts[ , c(14:52)] <- 0
  return(counts)
}

##############################################################################
#' scalarSearch
#' Searches through the posterior of an RJ continuous model for scalars and 
#' returns them.
#' @param rj_output partially processed RJ output.
#' @param counts The counts table.
#' @keywords internal
#' @name scalarSearch

scalarSearch <- function(rj_output, counts, fullmrcas, verbose) {
  alltypes <- allmrcas <- vector(mode = "list", length = nrow(rj_output))

  rates <- matrix(
    rep(1, nrow(counts) * nrow(rj_output)), ncol = nrow(rj_output)
  )
  rownames(rates) <- counts[ , "descNode"]

  # make lists for the origins of deltas etc.
  .tmp <- rep(1, nrow(rj_output))
  Node <- Branch <- Delta <- Lambda <- Kappa <- Node_effects <- replicate(
    nrow(counts), as.numeric(paste(.tmp)), simplify = FALSE
  )  
  names(Node) <- names(Branch) <- names(Delta) <- names(Lambda) <- names(Kappa) <- names(Node_effects) <- counts[ , "descNode"]

  print("Searching for scalars...")
    # I want to really turn this into a function so that I can use the 
    # same progress bar style - can it be done with apply?
    if (verbose) {
      pb <- pbapply::startpb(0, nrow(rj_output))
    }
    
    for (i in 1:nrow(rj_output)) {
      # This step is extremely slow - I wonder if keeping rj_output as a long 
      # list, rather than smartBinding it might help.
      lastrates <- rj_output[i, !is.na(rj_output[i, ])]
      
      # If the number of columns is seven, there are no scalars applied this 
      # generation.
      if (length(lastrates) == 7) {
        nodes <- NA
        scales <- NA
        types <- NA
      } else {
        
        int <- lastrates[8:length(lastrates)]
   
        nodes <- unlist(c(int[grep("NodeID*", names(int))]))
        scales <- unlist(c(int[grep("Scale*", names(int))]))
        types <- unlist(c(int[grep("NodeBranch*", names(int))]))
        mrcas <- sapply(
          nodes, function(x) fullmrcas[fullmrcas$node %in% x, "mrca"]
        )
        alltypes[[i]] <- types
        allmrcas[[i]] <- mrcas

        # Is this for-loop filling the scalar objects? Do I need to make them 
        # within this function? I wonder if this could be an after-the-fact 
        # function...
        for (j in seq_along(mrcas)) {
          nm <- paste0(
            types[j], "[[\"", as.character(mrcas[j]), "\"]]", "[", i, "]"
          )
          eval(parse(text = paste0(nm, "<-", scales[j])))
        }
      }
      if (verbose) { 
        pbapply::setpb(pb, i)    
      }
    }

    res <- list(alltypes = alltypes,
                allmrcas = allmrcas,
                rates = rates,
                Node = Node,
                Branch = Branch,
                Delta = Delta,
                Lambda = Lambda,
                Kappa = Kappa,
                Node_effects = Node_effects)
  return(res)
}

##############################################################################
#' multiplyNodes
#' Works out the cumulative effect of linear scalars on branches per iteration
#' @param scales A vector of scalars for a node
#' @param name The name of the node
#' @param tree The time tree
#' @param Node_effects A list, one element per node, to fill with the cumulative 
#' scalars
#' @name multiplyNodes
#' @keywords internal
multiplyNodes <- function(scales, name, tree, Node_effects) {
  # get descendents
  descs <- c(getDescs(tree, name), as.numeric(name))
  .tmp <- lapply(Node_effects[as.character(descs)], function(x) x * scales)
  return(.tmp)
}

##############################################################################
#' rjpp
#'
#' A function that takes the output of a kappa, lambda, delta, VRates etc. RJ 
#' bayesTraits run and runs post-processing on it.
#' @param rjlog The RJ output of the run - typically suffixed with .VarRates.txt
#' @param tree The time tree the analysis was run on as an object of class 
#' "phylo", or the filename of the timetree.
#' @param burnin The burnin (if required) for the mcmc (generally worked out 
#' from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out 
#' from the raw MCMC output logfile.
#' @param meanbranches If true, calculates mean, median and mode branch lengths 
#' and returns mean tree.
#' @param ratestable 
#' @import phytools pbapply ape
#' @export
#' @name rjpp

rjpp <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1, 
   ratestable = TRUE, verbose = TRUE) {
profvis({
  pbapply::pboptions(type = "timer", style = 3, char = "*")

  if (class(tree) == "phylo") {
    reftree <- ape::ladderize(tree)
  } else {
    reftree <- ape::ladderize(ape::read.nexus(tree))
  }

  if (verbose) {
    print("Loading log file.")
  }
  rjout <- loadRJ(rjlog, burnin = burnin, thinning = thinning)

  if (verbose) {
    print("Loading posterior trees.")
  }

  if (verbose) {
    print("Calculating mean branch lengths.")
  }

  tree_summary <- summariseTrees(reftree = reftree, trees = rjtrees, 
    burnin = burnin, thinning = thinning, verbose = verbose)

  rj_output <- rjout$rj_output
  subtrees <- rjout$subtrees
  rjtaxa <- rjout$taxa
  niter <- nrow(rj_output)
  
  if (verbose) {
    print("Finding taxa.")
  }
  
  if (verbose) {
    taxa <- pbapply::pblapply(subtrees$node, function(x) getTaxa(x, 
      subtrees = subtrees))
  } else {
    taxa <- lapply(subtrees$node, function(x) getTaxa(x, subtrees = subtrees))
  }
  
  if (verbose) {
    print("Calculating MRCAs.")
    fullmrcas <- unlist(
      pbapply::pblapply(taxa, function(x) {
        getMRCAbtr(x , tree = reftree, rjtaxa = rjtaxa)
      })
    )
  } else {
    fullmrcas <- unlist(
      lapply(taxa, function(x) {
        getMRCAbtr(x , tree = reftree, rjtaxa = rjtaxa)
      })
    )
  }
  
  fullmrcas <- data.frame(node = subtrees$node, mrca = fullmrcas) 
  counts <- createCountsTable(reftree, tree_summary)

  # Find the scalars.
  all_scalars <- scalarSearch(rj_output, counts, fullmrcas, verbose = verbose)

  # Calculate cumulative node effects
  for (i in 1:length(all_scalars$Node)) {
    .tmp <- multiplyNodes(all_scalars$Node[[i]], 
      names(all_scalars$Node)[i], 
      extree, 
      all_scalars$Node_effects)
    all_scalars$Node_effects[names(.tmp)] <- .tmp
  }

  all_scalars$Node_effects <- lapply(1:length(all_scalars$Node_effects), 
    function(x) all_scalars$Node_effects[[x]] * all_scalars$Branch[[x]])
  names(all_scalars$Node_effects) <- counts[ , "descNode"]
  
  origins <- list(nodes = do.call(rbind, all_scalars$Node), 
                  branches = do.call(rbind, all_scalars$Branch), 
                  delta = do.call(rbind, all_scalars$Delta),
                  lambda = do.call(rbind, all_scalars$Lambda), 
                  kappa = do.call(rbind, all_scalars$Kappa), 
                  rates = do.call(rbind, all_scalars$Node_effects)
                  )

  alltypes <- unlist(all_scalars$alltypes)
  allmrcas <- unlist(all_scalars$allmrcas)

  bs <- table(unlist(allmrcas)[alltypes == "Branch"])
  ns <- table(unlist(allmrcas)[alltypes == "Node"])
  ds <- table(unlist(allmrcas)[alltypes == "Delta"])
  ks <- table(unlist(allmrcas)[alltypes == "Kappa"])
  ls <- table(unlist(allmrcas)[alltypes == "Lambda"])

  bstaxa <- counts$descNode %in% names(bs)
  nstaxa <- counts$descNode %in% names(ns)
  dstaxa <- counts$descNode %in% names(ds)
  kstaxa <- counts$descNode %in% names(ks)
  lstaxa <- counts$descNode %in% names(ls)

  counts$nOrgnBRate[bstaxa] <- bs[match(counts$descNode[bstaxa], names(bs))] 
  counts$nOrgnScalar[bstaxa] <- counts$nOrgnBRate[bstaxa] + 
  bs[match(counts$descNode[bstaxa], names(bs))] 
  
  counts$nOrgnNRate[nstaxa] <- ns[match(counts$descNode[nstaxa], names(ns))]
  counts$nOrgnScalar[nstaxa] <- counts$nOrgnBRate[nstaxa] + 
  ns[match(counts$descNode[nstaxa], names(ns))]
  
  counts$nOrgnDelta[dstaxa] <- ds[match(counts$descNode[dstaxa], names(ds))]
  counts$nOrgnScalar[dstaxa] <- counts$nOrgnBRate[dstaxa] + 
  ds[match(counts$descNode[dstaxa], names(ds))]
  
  counts$nOrgnKappa[kstaxa] <- ks[match(counts$descNode[kstaxa], names(ks))]
  counts$nOrgnScalar[kstaxa] <- counts$nOrgnBRate[kstaxa] + 
  ks[match(counts$descNode[kstaxa], names(ks))]
  
  counts$nOrgnLambda[lstaxa] <- ls[match(counts$descNode[lstaxa], names(ls))]
  counts$nOrgnScalar[lstaxa] <- counts$nOrgnBRate[lstaxa] + 
  ls[match(counts$descNode[lstaxa], names(ls))]
  
  # Fill in transformation detail.

  counts[ , "meanDelta"] <- rowMeans(origins$delta)
  counts[ , "medianDelta"] <- apply(origins$delta, 1, median)
  counts[ , "modeDelta"] <-  apply(origins$delta, 1, modeStat)
  counts[ , "rangeDelta"] <- suppressWarnings(apply(origins$delta, 1, max) - 
    apply(origins$delta, 1, min))

  counts[ , "meanKappa"] <- rowMeans(origins$kappa)
  counts[ , "medianKappa"] <- apply(origins$kappa, 1, median)
  counts[ , "modeKappa"] <- apply(origins$kappa, 1, modeStat)
  counts[ , "rangeKappa"] <- suppressWarnings(apply(origins$kappa, 1, max) - 
    apply(origins$kappa, 1, min))

  counts[ , "meanLambda"] <- rowMeans(origins$lambda)
  counts[ , "medianLambda"] <- apply(origins$lambda, 1, median)
  counts[ , "modeLambda"] <- apply(origins$lambda, 1, modeStat)
  counts[ , "rangeLambda"] <- suppressWarnings(apply(origins$lambda, 1, max) - 
    apply(origins$lambda, 1, min))

  counts[ , "meanRate"] <- rowMeans(origins$rates)
  counts[ , "medianRate"] <- apply(origins$rates, 1, median)
  counts[ , "modeRate"] <- apply(origins$rates, 1, modeStat)
  counts[ , "rangeRate"] <- suppressWarnings(apply(origins$rates, 1, max) - 
    apply(origins$rates, 1, min))

  counts[ , "itersScaled"] <- 
  counts[ , "itersRatescaled"] <- apply(origins$rates, 1, function(x) {
    sum(x != 1)
  })
  counts[ , "itersDelta"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersKappa"] <- counts[ , "nOrgnDelta"]
  counts[ , "itersLambda"] <- counts[ , "nOrgnDelta"]

  counts[ , "pScaled"] <- 
  counts[ , "pRate"] <- apply(origins$rates, 1, function(x) sum(x != 1)) / niter
  counts[ , "pDelta"] <- counts[ , "nOrgnDelta"] / niter
  counts[ , "pKappa"] <- counts[ , "nOrgnKappa"] / niter
  counts[ , "pLambda"] <- counts[ , "nOrgnLambda"] / niter

  counts[ , "nScalar"] <- 
  counts[ , "nRate"] <- counts[ , "nOrgnNRate"] + counts[ , "nOrgnBRate"]
  counts[ , "nDelta"] <- counts[ , "nOrgnDelta"]
  counts[ , "nKappa"] <- counts[ , "nOrgnDelta"]
  counts[ , "nLambda"] <- counts[ , "nOrgnDelta"]

  # The big table 100% needs to be a tibble, and so do each of the per-iteration
  # scalar origins. With that in mind, I think the whole thing needs to have a 
  # class
  # and a custom print method - it just makes the most sense.

  counts <- counts[ , apply(counts, 2, function(x) all(x != 1))]
  counts <- counts[ , apply(counts, 2, function(x) all(x != 0))]

  if (meanbranches) {
    meantree <- extree
    meantree$edge.length <- counts[c(2:nrow(counts)) , "meanBL"]
  }

  if (all(sapply(origins$delta, function(x) all(x == 1)))) {
    origins$delta <- NULL
  }

  if (all(sapply(origins$kappa, function(x) all(x == 1)))) {
    origins$kappa <- NULL
  }

  if (all(sapply(origins$lambda, function(x) all(x == 1)))) {
    origins$lambda <- NULL
  }

  if (all(sapply(origins$nodes, function(x) all(x == 1)))) {
    origins$nodes <- NULL
  }  

  if (all(sapply(origins$branches, function(x) all(x == 1)))) {
    origins$branches <- NULL
  }

  if (all(sapply(origins$rates, function(x) all(x == 1)))) {
    origins$rates <- NULL
  }

  if (meanbranches) { 
    res <- list(data = counts, niter = niter, meantree = meantree)
  } else {
    res <- list(data = counts, niter = niter)
  }

  if (!is.null(origins$rates)) {
    scalars <- list(rates = origins$rates)
    origins$rates <- NULL
    res <- c(res, list(scalars = scalars))
  }

  res <- c(res, list(origins = origins))
})
  return(res)
}

