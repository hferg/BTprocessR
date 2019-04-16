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
#' @param reftree The time tree
#' @param tree_summary The output of summariseTrees
#' @name createCounts
#' @keywords internal

# NOTES
# This makes a massive and unwieldy table that is populated later on.
# Take a look at this and see if:
#    a) the headings are sensible and descriptive
#    b) can this table be split into more than one table (e.g. one for 
#      rates information, and one for branchlength information?)
#    c) If there are possible splits, can the occur here, or is it more
#      useful later on.
#    d) Is the table being built in the best way?
# The species descendent from each tip 100% needs to be a seperate list
# and it's own element in the results - this makes for ludicrous plotting.
# It should also have a method that prints the names properly when called,
# and when a specific node is selected.

createCountsTable <- function(reftree) {
  counts <- matrix(ncol = 42, nrow = (nrow(reftree$edge) + 1))
  colnames(counts) <- c(
    # identifying information
    "node_id", "branch", "ancNode", "descNode", "nTips",
    # branch length info
    "start", "end", "mid",
    # total scalar info
    "nScalar", "pScalar", "nOrgnScalar",
    # linear rate information
    "nRate", "pRate", "nOrgnNRate", "nOrgnBRate", "meanRate", "medianRate", 
    "modeRate", "rangeRate", "sdRate",
    # delta info
    "nDelta", "pDelta", "meanDelta", "medianDelta", "modeDelta",
    "rangeDelta", "sdDelta",
    # kappa info
    "nKappa", "pKappa", "meanKappa", "medianKappa", "modeKappa",
    "rangeKappa", "sdKappa",
    #lambda info
    "nLambda", "pLambda", "meanLambda", "medianLambda", "modeLambda",
    "rangeLambda", "sdLambda",
    # species - removed before returned to user
    "species"
  )

  # notes.
  # itersX and nX are effectively the same, except in the case of rates, where
  # the number of rates breaks down into branch and node. That means I don't
  # need the "iters" family of names.
  # make sure to fill in SD for each scalar.
  # equally I think the "origin" columns can be dropped - they are only there
  # for the node and branch scalars. Instead I could have the measurements for
  # node and branch scalars, and then the same for total scalar (the cumulative
  # effect of all linear scalars.)

    species_key <- vector(mode = "list", length = nrow(counts))
    counts[ , "node_id"] <- names(species_key) <- 
      paste0("node_", c(0:nrow(reftree$edge)))
    counts[ , "branch"] <- c(0:nrow(reftree$edge))

    counts[ , "ancNode"] <- c(0, reftree$edge[ , 1])
    counts[ , "descNode"] <- c((length(reftree$tip.label) + 1), 
      reftree$edge[ , 2])

    hts <- phytools::nodeHeights(reftree)
    hts <- round(abs(hts - max(hts)), 4)
    counts[ , "start"] <- c(0, hts[ , 1])
    counts[ , "end"] <- c(0, hts[ , 2])
    counts <- as.data.frame(counts)

    # Deal with the root
    counts[1, 1] <- "root"
    counts[1, 2] <- "root"
    descs <- getDescs(reftree, node = counts[1, "descNode"])
    counts[1, "nTips"] <- sum(descs <= length(reftree$tip.label))
    counts[1, "mid"] <- 0
    counts[1, "species"] <- paste0(reftree$tip.label[order(reftree$tip.label)], 
      collapse = ",")

    for (i in 2:nrow(counts)) {
      descs <- phytools::getDescendants(reftree, 
        node = as.numeric(counts[i, "descNode"])
      )
      counts[i, "nTips"] <- sum(descs <= length(reftree$tip.label))
      # if (counts[i, "nTips"] == 1) {
      #   counts[i, "nTips"] <- 1
      # }
      if (counts[i, "descNode"] <= length(reftree$tip.label)) {
        counts[i, "species"] <- reftree$tip.label[descs]
      } else {
        tips <- descs[descs <= length(reftree$tip.label)]
        tips <- reftree$tip.label[tips]
        counts[i, "species"] <- paste0(sort(tips), collapse = ",")
      }
      counts[i, "mid"] <- mean(c(hts[(i - 1), 1], hts[(i - 1), 2]))
      species_key[[i]]$nodes <- counts[i, 3:4]
      species_key[[i]]$species <- strsplit(counts[i, "species"], ",")[[1]]
    }
  species_key[[1]] <- list(nodes = "root", species = reftree$tip.label)
  counts <- counts[ , names(counts) != "species"]
  class(species_key) <- append("spkey", class(species_key))
  counts[ , c(9:42)] <- 0
  return(list(counts = counts, species_key = species_key))
}

##############################################################################
#' scalarSearch
#' Searches through the posterior of an RJ continuous model for scalars and 
#' returns them.
#' @param rj_output partially processed RJ output.
#' @param counts The counts table.
#' @param fullmrcas Most recent common ancestor tables for each tip in the tree,
#' generated internally.
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
  names(Node) <- names(Branch) <- names(Delta) <- names(Lambda) <- 
    names(Kappa) <- names(Node_effects) <- counts[ , "descNode"]

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
#' @param rjtrees 
#' @param tree The time tree the analysis was run on as an object of class 
#' "phylo", or the filename of the timetree.
#' @param burnin The burnin (if required) for the mcmc (generally worked out 
#' from the other logfile)
#' @param thinning Thinning parameter for the MCMC output - again, worked out 
#' from the raw MCMC output logfile.
#' @import phytools pbapply ape
#' @export
#' @name rjpp

rjpp <- function(rjlog, rjtrees, tree, burnin = 0, thinning = 1, 
   verbose = TRUE) {
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
  x <- createCountsTable(reftree)
  counts <- x$counts
  species_key <- x$species_key

  # Find the scalars.
  all_scalars <- scalarSearch(rj_output, counts, fullmrcas, verbose = verbose)

  # Calculate cumulative node effects. This involves propogating node scalars
  # down all branches that descend from that node.
  for (i in 1:length(all_scalars$Node)) {
    .tmp <- multiplyNodes(all_scalars$Node[[i]], 
      names(all_scalars$Node)[i], 
      tree, 
      all_scalars$Node_effects)
    all_scalars$Node_effects[names(.tmp)] <- .tmp
  }

  # And now multiply in any single-branch effects to get the cumulative linear
  # rate change across the whole tree of node + branch scalars.
  all_scalars$Node_effects <- lapply(1:length(all_scalars$Node_effects), 
    function(x) all_scalars$Node_effects[[x]] * all_scalars$Branch[[x]])
  names(all_scalars$Node_effects) <- counts[ , "descNode"]
  
  # The node effects object now becomes "rates".
  origins <- list(nodes = do.call(rbind, all_scalars$Node), 
                  branches = do.call(rbind, all_scalars$Branch), 
                  delta = do.call(rbind, all_scalars$Delta),
                  lambda = do.call(rbind, all_scalars$Lambda), 
                  kappa = do.call(rbind, all_scalars$Kappa), 
                  rates = do.call(rbind, all_scalars$Node_effects)
                  )

  alltypes <- unlist(all_scalars$alltypes)
  allmrcas <- unlist(all_scalars$allmrcas)

  # This gives tables of the numbers of each types of scalar applied to each 
  # node.
  bs <- table(unlist(allmrcas)[alltypes == "Branch"])
  ns <- table(unlist(allmrcas)[alltypes == "Node"])
  ds <- table(unlist(allmrcas)[alltypes == "Delta"])
  ks <- table(unlist(allmrcas)[alltypes == "Kappa"])
  ls <- table(unlist(allmrcas)[alltypes == "Lambda"])

  # this gives vectors with T/F for whether each branch got a branch, node, 
  # delta, kappa or lambda scalar. i.e. the TRUE in bstaxa is all taxa which
  # had a branch scalar originate on them.
  bstaxa <- counts$descNode %in% names(bs)
  nstaxa <- counts$descNode %in% names(ns)
  dstaxa <- counts$descNode %in% names(ds)
  kstaxa <- counts$descNode %in% names(ks)
  lstaxa <- counts$descNode %in% names(ls)

  # Now we want to count the number of origins of each scalar for each node.
  counts$nOrgnBRate[bstaxa] <- bs[match(counts$descNode[bstaxa], names(bs))]
  counts$nOrgnNRate[nstaxa] <- ns[match(counts$descNode[nstaxa], names(ns))]
  counts$nDelta[dstaxa] <- ds[match(counts$descNode[dstaxa], names(ds))]
  counts$nKappa[kstaxa] <- ks[match(counts$descNode[kstaxa], names(ks))]
  counts$nLambda[lstaxa] <- ls[match(counts$descNode[lstaxa], names(ls))]

  # Now add the number of scalar origins of ANY kind for each node. This is just
  # the sum of all the other nOrgn* columns.
  counts$nOrgnScalar <- counts$nOrgnBRate + counts$nOrgnNRate +
    counts$nDelta + counts$nKappa + counts$nLambda

  # Fill in transformation magnitude information.
  counts[ , "meanRate"] <- rowMeans(origins$rates)
  counts[ , "medianRate"] <- apply(origins$rates, 1, median)
  counts[ , "modeRate"] <- apply(origins$rates, 1, modeStat)
  counts[ , "rangeRate"] <- suppressWarnings(apply(origins$rates, 1, max) - 
    apply(origins$rates, 1, min))
  counts[ , "sdRate"] <- apply(origins$rate, 1, sd)
  counts[ , "nRate"] <- apply(origins$rates, 1, function(x) {
    sum(x != 1)
  })

  counts[ , "meanDelta"] <- rowMeans(origins$delta)
  counts[ , "medianDelta"] <- apply(origins$delta, 1, median)
  counts[ , "modeDelta"] <-  apply(origins$delta, 1, modeStat)
  counts[ , "rangeDelta"] <- suppressWarnings(apply(origins$delta, 1, max) - 
    apply(origins$delta, 1, min))
  counts[ , "sdDelta"] <- apply(origins$delta, 1, sd)

  counts[ , "meanKappa"] <- rowMeans(origins$kappa)
  counts[ , "medianKappa"] <- apply(origins$kappa, 1, median)
  counts[ , "modeKappa"] <- apply(origins$kappa, 1, modeStat)
  counts[ , "rangeKappa"] <- suppressWarnings(apply(origins$kappa, 1, max) - 
    apply(origins$kappa, 1, min))
  counts[ , "sdKappa"] <- apply(origins$kappa, 1, sd)

  counts[ , "meanLambda"] <- rowMeans(origins$lambda)
  counts[ , "medianLambda"] <- apply(origins$lambda, 1, median)
  counts[ , "modeLambda"] <- apply(origins$lambda, 1, modeStat)
  counts[ , "rangeLambda"] <- suppressWarnings(apply(origins$lambda, 1, max) - 
    apply(origins$lambda, 1, min))
  counts[ , "sdLambda"] <- apply(origins$lambda, 1, sd)

  counts[ , "nScalar"] <- counts[ , "nRate"] + counts[ , "nDelta"] +
    counts[ , "nKappa"] + counts[ , "nLambda"]

  # Now calculate the probability of being affected by a scaler.
  counts[ , "pScalar"] <- counts[ , "nScalar"] / niter
  counts[ , "pRate"] <- counts[ , "nRate"] / niter
  counts[ , "pDelta"] <- counts[ , "nDelta"] / niter
  counts[ , "pKappa"] <- counts[ , "nKappa"] / niter
  counts[ , "pLambda"] <- counts[ , "nLambda"] / niter

  counts <- counts[ , apply(counts, 2, function(x) !all(x == 1))]
  counts <- counts[ , apply(counts, 2, function(x) !all(x == 0))]
  
  # make any origins that aren't in the model NULL.
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

  # add branch numbers and node IDs to the origins and rates tables.
  for (i in seq_along(origins)) {
    origins[[i]] <- cbind(counts[ , 1:2], origins[[i]])
  }

  res <- list(
    scalars = tibble::as_tibble(counts), 
    tree_summary = tree_summary,
    species_key = species_key
  )

  if (!is.null(origins$rates)) {
    res <- c(res, list(rates = tibble::as_tibble(origins$rates)))
    origins$rates <- NULL
  }

  origins <- lapply(origins, tibble::as_tibble)

  res <- c(res, list(origins = origins), niter = niter)

  # NOTE: Plotshifts probably has too many different versions to just be a
  # a method actually... Although it could always just return a message if there
  # isn't, for example, a scalar or transformation selected?
  class(res) <- append("rjpp", class(res))
  return(res)
}

# NEW PLAN.
# If there is a custom print method for the output of rjpp, then I can have some
# outputs in the rjpp output which can be some tables of the information 
# required to plot the shifts and whatnot - then I just need to have some custom
# plot methods for the rjpp output - I think this is definitely worth trying...


#' summariseRjpp
#' This functions takes the somewhat massive output of the rjpp function and
#' pares it down to the scalars and/or rates the user is interested in.
#' QUESTION: Should this be limited to a single scalar type? That would include
#' the rates and the origins?
#' QUESTION: Is threshold important here? Or could that be in the autoplot
#' method?
#' @param pp An object of class "rjpp" - typically the output of the rjpp 
#' function.
#' @param scalar The scalar to summarise the results of. Either:
#' node_scalar, branch_scalar, rate, lambda, delta, kappa or node_branch
#' @param threshold The probability threshold over which scalars will be show.
#' When equal to zero ALL scalars in the posterior will be returned. When equal 
#' to 0.5 only scalars present in greater than 50% of posterior samples will be 
#' returned, and so on.

summariseRjpp <- function(PP, scalar) {

}



## TODO
# for some reason rates is empty, even from a variable rates analysis. This
# shouldn't be - this should contain the per-branch rates, surely?! Check!

