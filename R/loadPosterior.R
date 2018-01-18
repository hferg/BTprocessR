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

  raw <- readLines(logfile)
  # TODO Return the model type with the output, and put this into classes.
  # Adapt other functions to deal with the new output of btmcmc, and perhaps
  # implement methods based on the class (i.e. the model) that comes out of
  # this.
  model <- gsub(" ", "", raw[2])
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

