#' summariseTrees
#'
#' Summarise a posterior sample of trees from a rate-varible or RJ local transformation
#' BayesTraits MCMC analysis.
#' @param reftree A tree that provides the reference topology (in most cases this is
#' time-tree the analysis was run on). Can be a filename of a tree in the working 
#' directory or an object of class "phylo"
#' @param trees The posterior sample of trees from a rate-variable or RJlocaltransformation 
#' MCMC BayesTraits analysis. Typically will have the .Output.trees extension. Either 
#' the filename of the posterior, or an object of class "multiPhylo".
#' @param verbose If TRUE a progress bar will be shown when ladderizing the posterior
#' trees (this step can be time consuming).
#' @param burnin Number of trees to discard as burnin (if, for example, the MCMC 
#' chain had not converged until later in the run).
#' @param thinning If >1 then every nth tree will be sampled - useful if the sampling
#' interval of the original MCMC analysis was too small. Note that is it preferable
#' to ensure proper chain convergence prior to analysis of the results, in which case
#' the default settings of burnin and thinning will be appropriate.
#' @export
#' @name summariseTrees

summariseTrees <- function(reftree, trees, burnin = 0, thinning = 1, verbose = TRUE) {

  if (class(reftree) != "phylo") {
    reftree <- ape::read.nexus(reftree)
  }

  reftree <- ape::ladderize(reftree)
  
  if (class(trees) == "multiPhylo") {
    trees <- trees
  } else {
    trees <- ape::read.nexus(trees)
  }
  
  trees <- trees[seq.int(burnin, length(trees), thinning)]

  # ladderize trees.
  if (verbose) {
    print("Ladderizing posterior trees:")
    trees <- pbapply::pblapply(trees, ape::ladderize)
    class(trees) <- "multiPhylo"
  } else {
    trees <- lapply(trees, ape::ladderize)
    class(trees) <- "multiPhylo"
  }
  
  # and check topology
  for (i in seq_along(trees)) {
    if (sum(reftree$tip.label == trees[[i]]$tip.label) != length(reftree$tip.label)) {
      stop(paste("Tip labels on tree", i, "do not mactch reference tree"))
    }  
    if (sum(reftree$edge == trees[[i]]$edge) != length(reftree$edge)) {
      stop(paste("Tree", i, "has a different topology to reference tree"))
    }
  }

  bls <- sapply(trees, function(x) x$edge.length)
  
  meanbl <- apply(bls, 1, mean)
  medianbl <- apply(bls, 1, median)
  modebl <- apply(bls, 1, modeStat)
  sdbl <- apply(bls, 1, sd)
  rangebl <- apply(bls, 1, function(x) max(x) - min(x))

  meantree <- mediantree <- modetree <- reftree
  meantree$edge.length <- meanbl 
  mediantree$edge.length <- medianbl
  modetree$edge.length <- modebl

  summarytrees <- list(original_tree = reftree,
                        meant_ree = meantree,
                        median_tree = mediantree,
                        mode_tree = modetree)
  class(summarytrees) <- c("trees_summary", "multiPhylo")

  bls <- tibble(original_bls = reftree$edge.length,
                    mean_bls = meanbl,
                    median_bls = medianbl,
                    mode_bls = modebl,
                    range_bls = rangebl,
                    sd_bls = sdbl)

  res <- list(tree_summaries = summarytrees,
              branchlengt_info = bls)

  return(res)
}
