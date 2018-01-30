#' compPosts
#'
#' Compares histograms of a single parameter from two or more output files.
#' Generates a plot showing the distributions of the same parameter from two
#' or more posterior samples from BayesTraits MCMC.
#' @param logs A vector of EITHER two or more logfile names, OR a list of two 
#' or more objects of class bt_post (see \link[BTprocessR]{loadPosterior}).
#' @param parameter A vector containing the names of the parameter to compare.
#' @param alpha The transparency of the histograms.
#' @param ... Additional arguments passed to \link[BTprocessR]{loadPosterior},
#' in the case the filenames are given for logs.
#' @export
#' @examples

compPosts <- function(logs, parameter, alpha = 0.7, ...) {

  # If logs is class character, then it's filenames. Otherwise check
  # class.
  if (class(logs) == "character") {
    posteriors <- lapply(logs, function(x) loadPosterior(x, ...))
  } else {
    if (sum(sapply(logs2, class) == "bt_post") != length(logs)) {
      stop("logs must be either a vector of filenames, or a list of
        objects of class 'bt_post'.")
    }
    posteriors <- logs
  }

  # check the parameter is in all posteriors.
  if (sum(unlist(sapply(posteriors, colnames)) == parameter) != length(posteriors)) {
    stop(paste0("Parameter ", parameter, " is not present in all posteriors."))
  }

  cols <- viridis::viridis(length(posteriors))
  d <- lapply(posteriors, function(x) x[ , parameter])
  d <- reshape2::melt(d, id.vars = NULL)

  bwidth <- 3.5 * sd(d$value) * nrow(d) ^ -(1/3)
  p <- ggplot2::ggplot(d, ggplot2::aes(value, fill = factor(L1))) +
    ggplot2::geom_histogram(binwidth = bwidth, alpha = alpha) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    labs(
      x = parameter,
      fill = "Logfile"
    )
  p
}

