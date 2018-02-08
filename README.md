---
title: "BTprocessR"
output:
  github_document:
    toc: true
    toc_depth: 2
---



# Introduction

BTprocessR is an R package that provides a set of tools to help with the analysis of the output of the various MCMC models in [BayesTraits](http://www.evolution.rdg.ac.uk/BayesTraitsV3.0.1/BayesTraitsV3.0.1.html). The package includes functions for visualising the posterior distribtuion of the estimated parameters, summarising and plotting posterior distributions of phylogenies (resulting from rate-variable and/or reversible jump local transformation (RJLT) models, summarising inferred rate scalars for each node and branch in a tree, and identifying and plotting rate shifts. Currently the package only deals with the output of analyses of continuous traits, with functions for the various discrete trait analyses coming soon.

The package includes the output of three analyses of marsupuial body size - a brownian motion analysis, a variable rates analysis and an RJLT analysis using the delta transformation. In way of a disclaimer please note that in order to keep the package size small each of these MCMC chains collected only 500 samples, which may be insufficient to ensure proper convergence. The tree and data that these analyses used are also included should one wish to repeat these analyses with a longer chain or higher sampling rate.


```r
library(BTprocessR)
data(marsupial_data)
```

# Loading, summarising and plotting posteriors

A log file containing an MCMC posterior (by default BayesTraits appends these files with .Log.txt) can be parsed with the function loadPosterior. This function returns a [tibble](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html) with each sample taken during the MCMC chain, one per row. Printing this object returns some simple summary statistics for each parameter.


```r
library(BTprocessR)

post <- loadPosterior(system.file("extdata", "marsupials_brownian.txt.Log.txt", package = "BTprocessR"))

print(post)
```

```
## Posterior of  500  samples
## 
##   Parameter   Median     Mean     Mode    SD
## 1        Lh -140.112 -139.772 -139.377 1.018
## 2   Alpha.1     2.13    2.127    2.071 0.388
## 3 Sigma.2.1    0.011    0.011    0.011 0.001
## 
## # A tibble: 500 x 5
##    Iteration    Lh Tree.No Alpha.1 Sigma.2.1
##  *     <dbl> <dbl>   <dbl>   <dbl>     <dbl>
##  1   1020000  -139    1.00    1.84   0.0111 
##  2   1040000  -139    1.00    2.38   0.0112 
##  3   1060000  -140    1.00    2.18   0.00989
##  4   1080000  -139    1.00    2.19   0.0110 
##  5   1100000  -139    1.00    2.10   0.0113 
##  6   1120000  -140    1.00    2.60   0.0104 
##  7   1140000  -141    1.00    1.60   0.0122 
##  8   1160000  -140    1.00    2.60   0.0110 
##  9   1180000  -139    1.00    1.89   0.0103 
## 10   1200000  -140    1.00    2.56   0.00985
## # ... with 490 more rows
```

```r
# This information can also be assigned to an object for one or more parameters.

post_info <- getPostStats(post, parameter = "Lh")
post_info
```

```
## # A tibble: 1 x 5
##   parameter  mean median  mode    sd
##   <chr>     <dbl>  <dbl> <dbl> <dbl>
## 1 Lh         -140   -140  -139  1.02
```

```r
post_info <- getPostStats(post, parameter = c("Lh", "Sigma.2.1"))
post_info
```

```
## # A tibble: 2 x 5
##   parameter      mean    median      mode      sd
##   <chr>         <dbl>     <dbl>     <dbl>   <dbl>
## 1 Lh        -140      -140      -139      1.02   
## 2 Sigma.2.1    0.0111    0.0110    0.0110 0.00100
```

It is also possible to plot histograms of each of the parameters present in the posterior, and to plot some simple plots to aid in the visual diagnosis of convergence for either a specific parameter, or for specific parameter(s). The resulting plots are a histogram of the parameter samples, a trace plot, an autocorrelation plot, and a plot of the sliding mean of the samples with a window size of ten and a regression line (standard linear model).


```r
plot(post)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk plotPosterior](figure/plotPosterior-1.png)

```r
mcmcPlots(post)
```

![plot of chunk plotPosterior](figure/plotPosterior-2.png)

```r
mcmcPlots(post, parameters = "Lh")
```

![plot of chunk plotPosterior](figure/plotPosterior-3.png)

```r
mcmcPlots(post, parameters = c("Lh", "Sigma.2.1"))
```

![plot of chunk plotPosterior](figure/plotPosterior-4.png)

# Model comparison

## Visually comparing posteriors
BTprocessR contains functions that allow the quick comparison of like parameters over more than one posterior for the visual comparison of models. The legend on the resulting plot numbers the logfiles in the order that they were given to the function.


```r
post_brownian <- loadPosterior(system.file("extdata", "marsupials_brownian.txt.Log.txt", package = "BTprocessR"))

post_vrates <- loadPosterior(system.file("extdata", "marsupials_vrates.txt.Log.txt", package = "BTprocessR"))

compPosts(logs = list(post_brownian, post_vrates),  parameter = "Lh")
```

![plot of chunk comparePosteriors](figure/comparePosteriors-1.png)

```r
# Can take any number of posteriors.
post_vdelta <- loadPosterior(system.file("extdata", "marsupials_variabledelta.txt.Log.txt", package = "BTprocessR"))

compPosts(logs = list(post_brownian, post_vrates, post_vdelta),  parameter = "Lh")
```

![plot of chunk comparePosteriors](figure/comparePosteriors-2.png)

## Comparing marginal likelihoods

BayesTraits can calculate marginal likelhoods for a model by use of the stepping stone sampler. This is the recommended way to calculate marginal likelihoods for RJ models, and returns a file appended with .Stones.txt. BTprocessR can retrieve the marginal likelihoods from these files and return a table of marginal likelihoods as estimated by the stepping stone sampler. Plotting the resulting object will create a bayes factor plot (row vs column). It is worth providing labels for each model, especially if plotting or using a long filepath to the stones files otherwise the plot and/or table will be unwieldy.


```r
stones <- c(system.file("extdata", "marsupials_brownian.txt.Stones.txt", package = "BTprocessR"),
  system.file("extdata", "marsupials_vrates.txt.Stones.txt", package = "BTprocessR"),
  system.file("extdata", "marsupials_variabledelta.txt.Stones.txt", package = "BTprocessR"))

marginal_likelihoods <- getStones(stones, labels = c("Brownian motion", "Var. rates", "RJ delta"))

marginal_likelihoods
```

```
## # A tibble: 3 x 2
##   logfile         marginalLh
##   <chr>                <dbl>
## 1 Brownian motion       -155
## 2 Var. rates            -120
## 3 RJ delta              -120
```

```r
plot(marginal_likelihoods)
```

![plot of chunk marginalLikelihoods](figure/marginalLikelihoods-1.png)

# Summarising and plotting posterior samples of trees

One of the outputs of the various reversible jump (RJ) models in BayesTraits is a collection of transformed trees sampled during the MCMC process. By default BayesTraits appendeds these files with ".Output.trees". BTprocessR can read these posteriors in, summarise them, and plot the original tree alongisde the mean, median and modal trees from the posterior to visualise the difference.

The original tree (generally the time tree) is provided to enable a comparison between the posterior of trees and the original branch lengths.

The function returns a list of two elements - the first is a multiPhylo object containing the original tree, and trees transformed with the mean, median and mode branch lengths and the second element is a tibble containing the per-branch information.


```r
data(marsupial_trees)

post_trees <- summariseTrees(reftree = marsupial_trees$timetree,
                              trees = marsupial_trees$vrates_trees)
```

```
## [1] "Ladderizing posterior trees:"
```

```r
names(post_trees)
```

```
## [1] "tree_summaries"    "branchlength_info"
```

```r
print(post_trees)
```

```
## $tree_summaries
## 4 phylogenetic trees
## tree 1 : original_tree 246 tips
## tree 2 : mean_tree 246 tips
## tree 3 : median_tree 246 tips
## tree 4 : mode_tree 246 tips
## 
## $branchlength_info
## # A tibble: 490 x 6
##    original_bl mean_bl median_bl mode_bl range_bl   sd_bl
##          <dbl>   <dbl>     <dbl>   <dbl>    <dbl>   <dbl>
##  1       8.90   11.3       8.90    8.95     196    13.5  
##  2       7.80   25.4       7.80    7.91    1159    63.6  
##  3       2.10    6.82      2.10    2.27     232    18.6  
##  4       0.100   0.338     0.100   0.116     13.4   0.984
##  5       9.50  106        46.8    17.0     1241   152    
##  6       8.30   19.7       8.30    8.23     295    32.7  
##  7      10.9    24.0      10.9    10.7      651    42.3  
##  8      10.1    25.4      10.1     9.74     819    56.0  
##  9       0       0         0       0          0     0    
## 10       0       0         0       0          0     0    
## # ... with 480 more rows
```

```r
plot(post_trees$tree_summaries)
```

![plot of chunk treePosteriors](figure/treePosteriors-1.png)

```r
# Compare specific tree to time tree.
plot(post_trees$tree_summaries, tree = "median_tree")
```

![plot of chunk treePosteriors](figure/treePosteriors-2.png)

# Summarising and interpreting RJ output



<!-- ## Contents -->
<!-- 1. Introduction -->
<!-- 2. Loading, visualising and assessing posteriors -->
<!-- 3. Loading and visualising posterior samples of trees -->
<!-- 4. Post-processing rate-variable and reverse-jump local transformation models -->
<!-- + Postprocessing *.Varrates.txt files -->
<!-- + Finding rate shifts -->
<!-- + Visualising rate shifts  -->
