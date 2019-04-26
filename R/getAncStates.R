#' getAncStates
#' Reads the ancestral states from a posterior containing ancestral state
#' estimates and returns the mean, median and mode ancestral state, as well as 
#' the SD of estimates, for each of the tagged nodes in the posterior.
#' @param logfile The filename of a posterior that contains ancestral states.
#' @param tree The tree that the analysis was run on.
#' @name getAncStates
#' @export

getAncStates <- function() {
  # NOTES.
  # This SORT of doubles up on loadPosterior (and can piggyback some of the same
  # code). Is there a reason to keep them seperate from each other? Perhaps this
  # should just be built-in to loadPosterior?
  # Arguments NOT to build in - changes the S3 methods for plotting. If the 
  # output of load posterior contains both the standard posterior stuff and 
  # ancestral states stuff, then a plot() function to return ancestral states
  # cannot work well.
  # I for now keep seperate - the functions are quick in any case.  
}