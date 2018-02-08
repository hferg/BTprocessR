---
title: "BTprocessR"
author: "H. Ferguson-Gow"
date: "7 February 2018"
output:
  github_document:
    toc: true
    toc_depth: 2
---



# Introduction

BTprocessR is an R package that provides a set of tools to help with the analysis of the output of the various MCMC models in [BayesTraits](http://www.evolution.rdg.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html). The package includes functions for visualising the posterior distribtuion of the estimated parameters, summarising and plotting posterior distributions of phylogenies (resulting from rate-variable and/or reversible jump local transformation (RJLT) models, summarising inferred rate scalars for each node and branch in a tree, and identifying and plotting rate shifts. Currently the package only deals with the output of analyses of continuous traits, with functions for the various discrete trait analyses coming soon.

# Loading, summarising and plotting posteriors

A log file containing an MCMC posterior (by default BayesTraits appends these files with .Log.txt) can be parsed with the function loadPosterior. This function returns a tibble with each sample taken during the MCMC chain, one per row. Printing this object returns some simple summary statistics for each parameter.


```r
library(BTprocessR)

post <- loadPosterior(system.file("extdata", "marsupials_brownian.txt.Log.txt", package = "BTprocessR"))
```

```
## Warning in file(con, "r"): file("") only supports open = "w+" and open = "w
## +b": using the former
```

```
## Error in grep("\\bIteration\\b", raw):length(raw): argument of length 0
```

```r
print(post)
```

```
## Error in print(post): object 'post' not found
```

It is also possible to plot histograms of each of the parameters present in the posterior, and to plot some simple plots to aid in the visual diagnosis of convergence for either a specific parameter, or for specific parameter(s).


```r
plot(post)
```

```
## Error in plot(post): object 'post' not found
```

```r
mcmcPlots(post)
```

```
## Error in "bt_post" %in% class(logfile): object 'post' not found
```

```r
mcmcPlots(post, parameters = "Lh")
```

```
## Error in "bt_post" %in% class(logfile): object 'post' not found
```

```r
mcmcPlots(post, parameters = c("Lh", "Sigma.2.1"))
```

```
## Error in "bt_post" %in% class(logfile): object 'post' not found
```


<!-- ## Contents -->
<!-- 1. Introduction -->
<!-- 2. Loading, visualising and assessing posteriors -->
<!-- 3. Loading and visualising posterior samples of trees -->
<!-- 4. Post-processing rate-variable and reverse-jump local transformation models -->
<!-- + Postprocessing *.Varrates.txt files -->
<!-- + Finding rate shifts -->
<!-- + Visualising rate shifts  -->
