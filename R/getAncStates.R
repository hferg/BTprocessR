#' getAncStates
#' Reads the ancestral states from a posterior containing ancestral state
#' estimates and returns the mean, median and mode ancestral state, as well as 
#' the SD of estimates, for each of the tagged nodes in the posterior.
#' @param logfile The filename of a posterior that contains ancestral states.
#' @param tree The tree that the analysis was run on.
#' @name getAncStates
#' @export

getAncStates <- function(logfile, tree, burnin = 0, thinning = 1) {
  raw <- readLines(logfile)
  if (class(tree) != "phylo") {
    tree <- ape::ladderize(ape::read.nexus(tree))
  } else {
    tree <- ape::ladderize(tree)
  }

  if (length(grep("Tags", raw)) > 0) {
    ts <- grep("Tags", raw)
    te <- grep("Estimating", raw)
    ee <- grep("Restrictions", raw)
    tags <- raw[(ts+1):(te-1)]
    ests <- raw[(te + 1):(ee - 1)]
    tagssp <- strsplit(tags, "\t")
    tags <- vector(mode = "list", length = length(tagssp))
    for (i in seq_along(tagssp)) {
      tips <- strsplit(tagssp[[i]][4], " ")[[1]]
      tags[[i]] <- list(
        name = tagssp[[i]][2],
        node = ape::getMRCA(tree, tips),
        tips = tips) 
    }
    ests <- strsplit(ests, "\t")
    ests <- lapply(ests, function(x) strsplit(x[[2]], " ")[[1]][1])
    if (!all(sapply(tags, function(x) x$name) == ests)) {
      cat("Warning: Names of reconstructed nodes do not match names of tags.
        Assuming the order of reconstructed nodes matches the order of tag
        definitions.")
    }
    post <- loadPosterior(logfile, burnin = burnin, thinning = thinning)
    post <- post[, grep("Est", colnames(post))]
    colnames(post) <- sapply(colnames(post), function(x) {
      strsplit(x, "\\.")[[1]][2]
    })
    ancstates <- list(states = tibble::tibble(
      tag = colnames(post),
      node = sapply(tags, function(x) x$node),
      mean = colMeans(post),
      median = apply(post, 2, median),
      mode = apply(post, 2, modeStat),
      sd = apply(post, 2, sd)
    ),
    tree = tree
  )
  class(ancstates) <- append("ancstates", class(ancstates))
  } else {
    stop("No tags and reconstructed nodes found in logfile.")
  }
  return(ancstates)
}