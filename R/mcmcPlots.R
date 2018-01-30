#' mcmcPlots
#'
#' Plots 4 MCMC diagnostics for a parameter
#' @param logfile Either an object of class "bt_post", or the name of the
#' logfile of the analysis of interest.
#' @param parameters A vector of names of parameters to show diagnostic
#' plots for. If NULL (default) then all parameters are used (depending on
#' the analysis some of these may be nonsensical).
#' @param ... Additional arguments passed to loadPosterior
#' @keywords mcmc diagnostic autocorrelation trace running mean density
#' @name mcmcPlots
#' @importFrom ggplot2 autoplot
#' @export
#' @examples
#' paramDiagnostic("cool-data.log", "Lh")
#' params <- getParams("cool-data.log")
#' paramDiagnostic("cool-data.log", params, type = "trace", cols = 3)

mcmcPlots <- function(logfile, parameters = NULL, ...) {
  if ("bt_post" %in% class(logfile)) {
    posterior <- logfile
  } else { 
    posterior <- loadPosterior(logfile)
    # posterior <- loadPosterior(logfile, ...)
  }

  posterior <- posterior[, !colnames(posterior) %in% c("Iteration", "Harmonic.Mean", "Tree.No")]  
  if (!is.null(parameters)) {
    posterior <- posterior[, colnames(posterior) %in% parameters]
  } else {
    parameters <- colnames(posterior)
  }

  # Now build a plot list that will be passed to multiplot.

  ps <- vector(mode = "list", length = (length(parameters) * 4))
  for (i in seq_along(parameters)) {
    if (i == 1) {
      titles <- c("Hist.", "Trace", "AutoCorr", "Sliding mean")
    } else {
      titles <- c("", "", "", "")
    }
    ps[[i * 4 - 3]] <- ggHist(posterior, parameters[i], title = titles[1])
    ps[[i * 4 - 2]] <- ggTrace(posterior, parameters[i], title = titles[2])
    ps[[i * 4 - 1]] <- ggAutoCor(posterior, parameters[i], title = titles[3])
    ps[[i * 4]] <- ggRunmean(posterior, parameters[i], title = titles[4])
  }

  # Then arrange these using gridExtra.
  if (length(parameters) == 1) {
    ncol <- 2
    lmatrix <- matrix(1:length(ps), ncol = 2)
  } else {
    ncol <- 4
    lmatrix <- matrix(1:length(ps), ncol = length(parameters))
  }
  gridExtra::grid.arrange(grobs = ps, layout = t(lmatrix), ncol = ncol)
}

