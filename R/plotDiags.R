#' ggHist
#' 
#' A function to plot a nice histogram - one that uses ggplot2
#' but has sensible break behaviour in R. Produces a single
#' histogram but is mostly called internally for things like
#' plotPosts.
#' @name ggHist
#' @param posterior The posterior of a bayestraits analysis (from loadPosterior)
#' @param param The parameter to plot a histogram of.
#' @return A histogram with nice ggplot aesthetics, and good bin widths.

ggHist <- function(posterior, param, title) {
  d <- posterior[ , colnames(posterior) == param]
  bwidth <- 3.5 * sd(data.frame(d)[,1], na.rm = TRUE) * nrow(d) ^ -(1/3)
  col <- viridis::viridis(4)[1]
  z <- ggplot2::ggplot(d, ggplot2::aes_string(param)) +
    ggplot2::geom_histogram(fill = col, binwidth = bwidth) +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::ggtitle(title)
  return(z)
}

#' ggTrace
#'
#' Make a trace plot for an MCMC parameter using ggplot2.
#' @param posterior BayesTraits MCMC output - output from loadPosterior function.
#' @param param The paramater you want to see an autocorrelation plot for.

ggTrace <- function(posterior, param, title){
  d <- posterior[ , colnames(posterior) == param]
  d <- data.frame(d)[ , 1]
  col <- viridis::viridis(4)[2]
  d <- data.frame(Iteration = c(1:length(d)), p = d)
  z <- ggplot2::ggplot(d, ggplot2::aes(x = Iteration, y = p)) +
    ggplot2::geom_line(colour = col) + 
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::ggtitle(title)
  return(z)
}

#' ggAutoCor
#'
#' Makes an autocorrealtion plot using ggplot2 - which looks nicer than the
#' core R plot functions.
#' @param psoterior The data for the plot, in an object of class vector.
#' @param param The parameter to do the autocorrelation plot for.
#' @param conf The confidence level you want to see the limits for.
#' @param max.lag Maximum lag
#' @param min.lag Minimum lag

ggAutoCor <- function(posterior, param, max.lag = NULL, min.lag = 0, title) {
  d <- posterior[ , colnames(posterior) == param]
  confline <- qnorm((1 - 0.95) / 2) / sqrt(length(d))
  x <- acf(d ,plot = FALSE, lag.max = NULL)
  xacf <- with(x ,data.frame(lag, acf))
  xacf <- xacf[-seq(1 ,1), ]
  sig <- (abs(xacf[ ,2]) > abs(confline)) ^ 2
  col <- viridis::viridis(4)[3]
  z <- ggplot2::ggplot(xacf, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_bar(color = "darkgray", stat = "identity", position = "identity", fill = col) +
    ggplot2::geom_hline(yintercept = -confline, color = "blue", size = 0.2) +
    ggplot2::geom_hline(yintercept = confline, color = "blue", size = 0.2) +
    ggplot2::geom_hline(yintercept = 0, color = "red", size = 0.2) +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::ggtitle(title)
  return(z)
}

#' ggRunmean
#'
#' Make a running mean plot for an MCMC parameter using ggplot2.
#' @param dat The vector of the paramater you want to see an autocorrelation plot for.

ggRunmean <- function(posterior, param, title, window.size = 10) {
  d <- posterior[ , colnames(posterior) == param]
  d <- data.frame(d)[, 1]
  windows <- seq.int(window.size, length(d), window.size)
  
  if (windows[length(windows)] != length(d)) {
    windows <- c(windows, length(d))
  }

  means <- vector(mode = "numeric", length = length(windows))
  for (i in 1:length(windows)) {
    means[i] <- mean(d[c((windows[i] - (window.size -1)):windows[i])])
  }

  col <- viridis::viridis(4)[4]
  z <- ggplot2::ggplot(data.frame(Window = c(1:length(windows)), mean = means), 
        ggplot2::aes(x = Window, y = mean)) +
    ggplot2::geom_point(color = "black", fill = col, pch = 21) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) + 
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::ggtitle(title)
  return(z)
}

