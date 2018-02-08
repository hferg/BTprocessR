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

  if (!"bt_post" %in% class(logfile)) {
    posterior <- loadPosterior(logfile, ...)
  } else  {
    posterior <- logfile
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
