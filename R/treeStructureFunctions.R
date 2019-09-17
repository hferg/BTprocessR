################################################################################
#' getDescs
#'
#' A function to get all descendant nodes from a given node, or vector of tip 
#' labels.
#' @param tree A tree of class phylo
#' @param node Either a single node, or a vector of tip labels
#' @name getDescs
#' @export

# something is wrong with this function - a) it seems to call itself (weird?!)
# and b) it's returning some nodes with no tip-descendants, which should never
# happen.

getDescs <- function(tree, node, nds = NULL) {
  # If more than one node is specified, find the MRCA of those nodes to get
  # descendents from.
  if (length(node) > 1) {
    node <- getMRCA(tree, node)
  }
  
  # make a nodes vector.
  if (is.null(nds)) {
    nds <- vector()
  }
  
  # get the descendent
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

################################################################################
#' sharedBranches
#'
#' Make species key from a tree.
#' @param tree An object of class "phylo"
#' @name makeKey
#' @keywords internal

makeKey <- function(tree) {
  tb <- matrix(ncol = 2, nrow = (nrow(tree$edge) + 1))
  colnames(tb) <- c("ancNode", "descNode")
  key <- vector(mode = "list", length = nrow(tb))
  names(key) <- paste0("node_", c(0:nrow(tree$edge)))
  tb[ , "ancNode"] <- c(0, tree$edge[ , 1])
  tb[ , "descNode"] <- c((length(tree$tip.label) + 1), 
    tree$edge[ , 2])

  for (i in 2:nrow(tb)) {
    descs <- phytools::getDescendants(tree, 
      node = as.numeric(tb[i, "descNode"])
    )
    if (as.numeric(tb[i, "descNode"]) <= length(tree$tip.label)) {
      sp <- tree$tip.label[descs]
    } else {
      tips <- descs[descs <= length(tree$tip.label)]
      tips <- tree$tip.label[tips]
      sp <- paste0(sort(tips), collapse = ",")
    }
    key[[i]]$nodes <- tb[i, 1:2]
    key[[i]]$species <- strsplit(sp, ",")[[1]]
  }
  key[[1]] <- list(nodes = "root", species = tree$tip.label)
  return(key)
}

################################################################################
#' sharedBranches
#'
#' Finds shared branches between two trees where the smaller tree is a subset of
#' the larger tree.
#' @param tree1 An object of class "phylo"
#' @param tree2 An object of class "phylo"
#' @name sharedBranches
#' @export

sharedBranches <- function(tree1, tree2) {
  # find the biggest tree.
  if (tree1$Nnode > tree2$Nnode) {
    big <- tree1
    small <- tree2
  } else if (tree1$Nnode < tree2$Nnode) {
    big <- tree2
    small <- tree1
  } else if (tree1$Nnode == tree2$Nnode) {
    stop("Both trees have the same number of nodes.")
  }

  # make sure that the smaller tree is a subset of the larger.
  if (any(!small$tip.label %in% big$tip.label)) {
    stop("Smaller tree is not a subset of larger tree.")
  }
  
  bns <- big$tip.label[!big$tip.label %in% small$tip.label]
  big_key <- makeKey(big)
  small_key <- makeKey(small)
  shared_table <- matrix(nrow = length(shared), ncol = 7)
  colnames(shared_table) <- c("lg_nodeid", "lg_ancNode", "lg_descNode", 
    "shared",
    "sml_nodeid", "sml_ancNode", "sml_descNode")
  for (i in seq_along(big_key)) {
    shared[[i]]$nodes <- big_key[[i]]$nodes
    sp <- big_key[[i]]$species
    sp <- sp[!sp %in% bns]
    big_key[[i]]$small_sp <- sp
    if (i == 1) {
      shared_table[i, "lg_nodeid"] <- "root"
      shared_table[i, "lg_ancNode"] <- big_key[[i]]$nodes[1]
      shared_table[i, "lg_descNode"] <- big_key[[i]]$nodes[1]
    } else {
      shared_table[i, "lg_nodeid"] <- names(big_key)[i]
      shared_table[i, "lg_ancNode"] <- big_key[[i]]$nodes["ancNode"]
      shared_table[i, "lg_descNode"] <- big_key[[i]]$nodes["descNode"]
    }
      
    if (length(sp) != 0) {
      shared_table[i, "shared"] <- TRUE
      b_sp <- paste(sort(big_key[[i]]$small_sp), collapse = " ")
      # now put in the anc/dec nodes of the smaller tree for the shared branch.
      # find the branch from small_key that has matching species.
      sh <- sapply(small_key, function(x) {
        paste(sort(x$species), collapse = " ") == b_sp
      })
      sml_match <- small_key[[which(sh)]]
      shared_table[i, "sml_nodeid"] <- names(which(sh))
      if (i == 1) {
        shared_table[i, "sml_ancNode"] <- "root"
        shared_table[i, "sml_descNode"] <- "root"
      } else {
        shared_table[i, "sml_ancNode"] <- sml_match$nodes["ancNode"]
        shared_table[i, "sml_descNode"] <- sml_match$nodes["descNode"]
      }
        
    } else {
      shared_table[i, "shared"] <- FALSE
      shared_table[i, "sml_ancNode"] <- NA
      shared_table[i, "sml_descNode"] <- NA
    }
  }
  shared_table <- as.data.frame(shared_table, stringsAsFactors = FALSE)
  shared$shared <- as.logical(shared)
  return(tibble::as_tibble(shared_table))
}
