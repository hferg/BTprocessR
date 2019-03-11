################################################################################
#' getDescs
#'
#' A function to get all descendant nodes from a given node, or vector of tip 
#' labels.
#' @param tree A tree of class phylo
#' @param node Either a single node, or a vector of tip labels
#' @name getDescs
#' @export

getDescs <- function(tree, node, nds = NULL) {

  if (length(node) > 1) {
    node <- getMRCA(tree, node)
  }
  
  if (is.null(nds)) {
    nds <- vector()
  }
    
  dtrs <- tree$edge[which(tree$edge[ , 1] == node), 2]
  nds <- c(nds, dtrs)
  now <- which(dtrs >= length(tree$tip))
  
  if (length(now) > 0) {
  
    for (i in 1:length(now)) {
      nds <- getDescs(tree, dtrs[now[i]], nds)
    }
    
  }
  return(nds)
}

################################################################################
#' getPL
#' Gets the path length from the root to the position of a particular node.
#' @export
#' @name getPL
#' @keywords internal

getPL <- function(tree, startnode = NA, node) {

  # Get all paths for this tree
  allpaths <- ape::nodepath(tree)
  
  # get all paths containing this node
  if(node > ape::Ntip(tree)) {

    paths <- allpaths[grepl(node, allpaths)]

  } else {

    paths <- allpaths[[node]]
  
  }
  
  # If the terminal node we care about is not a terminal, we want to remove the 
  # branches after our node
  if(node > ape::Ntip(tree)) {

    paths <- lapply(paths, function(x) x <- x[x > ape::Ntip(tree)])
    path <- unlist(unique(lapply(paths, function(x) x <- x[x <= node]))) 
  
  } else {

    path <- unlist(unique(paths))
  
  }
  # If the startnode is not NA, then we want to trim the paths. 
  if(!is.na(startnode)) {              
    
    while(path[1] != startnode) {
      path <- path[-1]
    }
    
    path <- path[-1] 
  } 
  # Remove the start node [ otherwise it puts the leading branch in there too ]
  # get the path length
  
  branches <- which(tree$edge[,2] %in% path)
  distance <- sum(tree$edge.length[branches])
  return(distance)
}

################################################################################
#' getTipNames
#'
#' A function to get the names of the descendant tips from a given node of a 
#' tree.
#' @param tree A tree of class phylo.
#' @param node The node number of interest.
#' @export

getTipNames <- function(tree, node) {
  descs <- getDescs(tree, node)
  descs <- descs[descs <= length(tree$tip.label)]
  tree$tip.label[descs]
}

################################################################################
#' getTaxa
#'
#' This gets the taxa names of a particular subtree from a list of all subtrees 
#' comprising a full tree.
#' @param subtrees A list of subtrees, as written by BayesTraits, and read in 
#' during post-processing
#' @param node The node number of interest.
#' @name getTaxa
#' @keywords internal

getTaxa <- function(x, subtrees) {
  taxa <- subtrees[subtrees$node == x, ]
  taxa <- taxa[ , !is.na(taxa)]
  taxa <- taxa[c(4:length(taxa))]
  return(as.numeric(unlist(taxa)))
}

################################################################################
#' getMRCAbtr
#'
#' This an extension of apes's getMRCA that enables the return of a tip, or an 
#' MRCA. Translates taxa codes (BayesTraits) to proper tip labels. Useful only 
#' in post-processing.
#' @param x A vector of taxa names
#' @param tree A phylogeny of class "phylo" (generally the time tree used as 
#' input to BayesTraits)
#' @param rjtaxa The taxa translations as output from BayesTraits
#' @name getMRCAbtr
#' @keywords internal

getMRCAbtr <- function(x, tree, rjtaxa) {
  if (length(x) == 1) {
    mrca <- which(tree$tip.label == rjtaxa[rjtaxa[ , 1] %in% x, 2])
  } else {
    mrca <- ape::getMRCA(tree, rjtaxa[rjtaxa[ , 1] %in% x, 2])
  }
  return(mrca)
}

