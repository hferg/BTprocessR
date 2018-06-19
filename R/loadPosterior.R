#' loadPosterior
#'
#' Returns the full mcmc object from a BayesTraits log file. This
#' is used inside plot functions and so on, but might be useful for
#' other MCMC manipulations and so on.
#' Extracts the MCMC samples from a BayesTraits logfile (i.e. discards the
#' header information and coerces samples into a matrix.)
#' @param logfile The name of the logfile of the BayesTraits analysis.
#' @param thinning Thinning parameter for the posterior - defaults to 1
#' (all samples). 2 uses every second sample, 3 every third and so on.
#' @param burnin The number of generations to remove from the start of the
#' chain as burnin. Use if the chain has not reached convergence before sampling
#' began. Useful if the burnin parameter for the analysis itself was not long
#' enough.
#' @return A tibble (see \link[tibble]{tibble}) with the class "bt_post" containing
#' the samples from the BayesTraits MCMC chain. Headers vary on model type.
#' @export
#' @name loadPosterior

loadPosterior <- function(logfile, thinning = 1, burnin = 0) {

  # TODO Bit of an overhaul - load the logfile, and return any specified 
  # local scalar information that may have been specified. Make a class
  # for this output that then feeds into the other functions. This is important
  # to do, since currently those specified scalars are not included in the 
  # varrates output, which means that they are ignored. This DOES involve a bit
  # more work on the other functions, however (i.e. including the new object 
  # class, making a print method for the new object class etc. etc.)

  # It would make sense to provide the prefix and this functions loads the 
  # posterior, the trees, the varrates file and the stones file (if present)
  # and does an all-in-one. This object can then be variously post processed
  # the other post-processing functions. This is a much neater way to do this
  # using the S3 methods I believe.

  raw <- readLines(logfile)
  # TODO Return the model type with the output, and put this into classes.
  # Adapt other functions to deal with the new output of btmcmc, and perhaps
  # implement methods based on the class (i.e. the model) that comes out of
  # this.
  # Get model type.
  model <- gsub(" ", "", raw[2])
  # find tags.

  # Actually... this doesn't matter for tags, just for Local Rates?
  if (length(grep("Tags", raw)) > 0) {
    # Get the tags (i.e. the names of them, and the tips that define them)
      # get tag line
    
  }

  output <- do.call(rbind, strsplit(raw[grep("\\bIteration\\b", raw):length(raw)], "\t"))
  colnames(output) <- output[1, ]
  output <- output[c(2:nrow(output)), ]
  output <- data.frame(output, stringsAsFactors = FALSE)

  for (i in 1:ncol(output)) {
    if (colnames(output)[i] != "Model.string" && colnames(output)[i] != "Dep...InDep") {
      output[ ,i] <- as.numeric(output[ ,i])
    }

  }
  output <- tibble::as_tibble(output[seq.int(burnin, nrow(output), thinning), ])
  class(output) <- append("bt_post", class(output))
  return(output)
}

getTags <- function(raw) {
  tagstart <- grep("Tags", raw) + 1

  # Get end of tags block (will be EITHER "Local Rates" or "Restrictions")
  # Get all of the local rates and the names of them.

  if (length(grep("Local Rates", raw)) > 0) {
    tagend <- grep("Local Rates", raw) - 1
  } else {
    tagend <- grep("Restrictions", raw) - 1
  }

  tagblock <- raw[tagstart:tagend]

  tags <- lapply(tagblock, function(x) {
    .x <- strsplit(tagblock, "\t")[[1]]
    strsplit(.x[4], " ")[[1]]
  })

  names(tags) <- sapply(tagblock, function(x) {
    strsplit(tagblock, "\t")[[1]][2]
  })

  # Now deal with the local rates. What is needed is a) the name of the
  # local transformation, b) the type of scalar applied (node, branch, delta,
  # etc. etc.), c) whether it is fixed (and if so, the value) or estimated, 
  # and d) the name of the tag that it is applied to.

  # THEN swap in the name of the tag for the tips.

}