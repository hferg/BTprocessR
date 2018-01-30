##############################################################################
#' modeStat
#' Calculate the mode.
#' @name modeStat
#' @keywords internal

modeStat <- function(x, log = FALSE) 
{  
  if (length(x) < 2) {
    return("NA")
  }

  if (all(x == x[1])) {
    return(x[1])
  }

  if (log == TRUE) {
    dens <- density(log(x))
    mode.tmp <- dens$x[which(dens$y == max(dens$y))]
    mode <- exp(mode.tmp)
  } else {
    dens <- density(x)
    mode <- dens$x[which(dens$y == max(dens$y))]
  }
  return(mode)
}
