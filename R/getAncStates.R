#' getAncStates
#' Reads the ancestral states from a posterior containing ancestral state
#' estimates and returns the mean, median and mode ancestral state, as well as 
#' the SD of estimates, for each of the tagged nodes in the posterior.
#' @param logfile The filename of a posterior that contains ancestral states.
#' @param tree The tree that the analysis was run on.
#' @name getAncStates
#' @export

getAncStates <- function(logfile, tree, burnin = 1, thinning = 0) {
  # NOTES.
  # This SORT of doubles up on loadPosterior (and can piggyback some of the same
  # code). Is there a reason to keep them seperate from each other? Perhaps this
  # should just be built-in to loadPosterior?
  # Arguments NOT to build in - changes the S3 methods for plotting. If the 
  # output of load posterior contains both the standard posterior stuff and 
  # ancestral states stuff, then a plot() function to return ancestral states
  # cannot work well.
  # I for now keep seperate - the functions are quick in any case.

  # read the posterior in - does loadPosterior already retain all column
  # headers? It DOES!

  # The reconstructed values DO seem ludicrous... but I think just ignore that
  # for now and get on with the function and then see what the deal is later
  # down the line... It might be a problem for the work Jon is doing though?


  raw <- readLines(logfile)
  if (class(tree) != "phylo") {
    tree <- ape::read.nexus(tree)
  }

  if (length(grep("Tags", raw)) > 0) {
    # Get the tags (i.e. the names of them, and the tips that define them)
      # get tag line
    ts <- grep("Tags", raw)
    te <- grep("Estimating", raw)
    ee <- grep("Restrictions", raw)

    # in this example the tags are the same as the node reconstruction names -
    # check if this is the same if I rename them...

    # It actually DOESN'T - that means that there isn't a good way to match
    # the tags with the reconstructed nodes... 
    
    # pop in a warning/error that if the reconstructed node names don't match
    # tag names then the function is assuming they are in the same order.

    # cut out the tag section
    tags <- raw[(ts+1):(te-1)]

    # and the "estimating" section.
    ests <- raw[(te + 1):(ee - 1)]

    # then use string split to divide the tags up into the tag name, and the
    # descendent tips.
    tagssp <- strsplit(tags, "\t")
    tags <- vector(mode = "list", length = length(tagssp))
    for (i in seq_along(tagssp)) {
      tips <- strsplit(tagssp[[i]][4], " ")[[1]]
      tags[[i]] <- list(
        name = tagssp[[i]][2],
        node = ape::getMRCA(tree, tips),
        tips = tips) 
    }

    # and the same to divide ests up into the reconstructed node names.
    ests <- strsplit(ests, "\t")
    ests <- lapply(ests, function(x) strsplit(x[[2]], " ")[[1]][1])

    # check that ests match with the names in tags
    if (!all(sapply(tags, function(x) x$name) == ests)) {
      cat("Warning: Names of reconstructed nodes do not match names of tags.
        Assuming the order of reconstructed nodes matches the order of tag
        definitions.")
    }

    # Now read the reconstructions from the posterior - this can come from
    # loadPosterior.
    post <- loadPosterior(logfile, burnin = burnin, thinning = thinning)

    # now take off columns that aren't reconstructions.
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