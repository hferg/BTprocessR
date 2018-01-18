#' getPostStats
#'
#' Given a posterior of MCMC samples from a BayesTraits analysis (a logfile -
#' typically with the .log.txt extension) and the name of one or morea
#' parameters present in that logfile this will return the mean, median, mode
#' and standard deviation of the parameter estimates.
#' @param logfile The name of the logfile of the BayesTraits analysis, or an
#' object of class "bt_post" (see \link[BTprocessR]{loadPosterior})
#' @param parameter The name(s) of one or more parameters present in the logfile.
#' @param ... Additional arguments passed to \link[BTprocessR]{loadPosterior}.
#' @return An object of class data.frame with a row for each parameter, and
#' columns for the mean, median, mode and standard deviation of those
#' parameters.
#' @export
#' @examples
#' params <- getParams("cool-data.log.txt")
#' getPostMean("cool-data.log.txt", params[c(2:5)]
#' getPostMean("cool-data.log.txt", "Lh")

getPostStats <- function(logfile, parameter, ...) {

  modeStat <- function(x) {
    z <- unique(x)
    x[which.max(tabulate(match(x, z)))]
  }

  if (!"bt_post" %in% class(logfile)) {
    posterior <- loadPosterior(logfile, ...)
  }

  d <- posterior[ , colnames(posterior) %in% parameter]
  res <- matrix(ncol = 5, nrow = length(parameter))
  colnames(res) <- c("parameter", "mean", "median", "mode", "sd")
  res <- data.frame(res)
  for (i in seq_along(parameter)) {
    res[i, 1] <- parameter[i]
    res[i, 2] <- round(mean(d[parameter[i]][[1]]), 4)
    res[i, 3] <- round(median(d[parameter[i]][[1]]), 4)
    res[i, 4] <- round(modeStat(d[parameter[i]][[1]]), 4)
    res[i, 5] <- round(sd(d[parameter[i]][[1]]), 4)
  }

  return(tibble::as_tibble(res))
}


#' getStones
#'
#' Get the marginal likelihoods from from one or more stepping stones log
#' files.
#' @param stonesfile A vector of one or more log files from stepping stones analysis.
#' @param order If TRUE then the resultant data.frame is ordered according to
#' marginal likelihood.
#' @name getStones
#' @return A data.frame with a row per stones file, and a column with the
#' marginal likelihoods of each model in.
#' @export

getStones <- function(stonesfile, order = TRUE) {
  # TODO A function to count the parameters in the original model and return
  # them as well as the marginal likelihoods (to look at improvement in
  # marginal Lh as parameter numbers increase).
  res <- matrix(ncol = 2, nrow = length(stonesfile))
  colnames(res) <- c("logfile", "marginalLh")
  res <- data.frame(res)
  for (i in 1:length(stonesfile)) {
    raw <- readLines(stonesfile[[i]])
    res[i, 1] <- gsub(".txt.log.txt.Stones.txt", "", stonesfile[i])
    res[i, 2] <- as.numeric(strsplit(raw[length(raw)], "\t")[[1]][2])
  }
  res <- data.frame(res)
  res$marginalLh <- as.numeric(as.character(res$marginalLh))

  if (order) {
    res <- res[order(res$marginalLh, decreasing = TRUE), ]
  }

  return(as_tibble(res))
}
