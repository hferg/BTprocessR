#' getStones
#'
#' Get the marginal likelihoods from from one or more stepping stones log
#' files.
#' @param stones A vector of one or more log files from stepping stones 
#' analysis.
#' @param labels A vector of equal length to the number of stones files 
#' containing a name for each model. If left NULL then the filenames will be 
#' used (which may be long and unwieldy when plotting).
#' @name getStones
#' @return A data.frame with a row per stones file, and a column with the
#' marginal likelihoods of each model in.
#' @export

getStones <- function(stones, labels = NULL) {
  # TODO A function to count the parameters in the original model and return
  # them as well as the marginal likelihoods (to look at improvement in
  # marginal Lh as parameter numbers increase).

  # TODO Make a matrix of all supplied stones files, and calculate the 
  # bayes factors for all pairs, then return a plot - have an argument to 
  # either return the plot, or the table.
  res <- matrix(ncol = 2, nrow = length(stones))
  if (!is.null(labels)) {
    if (length(labels) != length(stones)) {
      stop("The number of labels must equal the number of stones files.")
    }
  }
  colnames(res) <- c("logfile", "marginalLh")
  for (i in 1:length(stones)) {
    raw <- readLines(stones[[i]])
    res[i, 1] <- stones[[i]]
    res[i, 2] <- as.numeric(strsplit(raw[length(raw)], "\t")[[1]][2])
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  if (!is.null(labels)) {
    res$logfile <- labels
  }
  res$marginalLh <- as.numeric(as.character(res$marginalLh))
  res <- tibble::as_tibble(res)
  class(res) <- append("bt_stones", class(res))
  return(res)
}
