% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscFunctions.R
\name{smartBind}
\alias{smartBind}
\title{smartBind}
\usage{
smartBind(...)
}
\arguments{
\item{A}{bunch of vectors.}
}
\description{
A function that will rbind vectors of different lengths and return a matrix, 
provided each vector element is named.
}
\examples{
do.call(smartRbind, list.of.vectors)
}
