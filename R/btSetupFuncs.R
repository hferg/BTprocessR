#' makeInfile
#' Takes a tree and some data, ladderizes the tree, matches the tree and the
#' data and then makes a run file based on a control list. The prepared tree,
#' data, infile and command to execute the analysis will be written out to
#' the analysis_path with the prefixes analysis_name.
#' Sets up basic
#' BayesTraits analysis - for more complex analyses (e.g. alternate priors,
#' parameter restriction, fixed scalars on nodes, fossilised states etc.) users 
#' must write their own control files.
#' @param tree The tree for the analysis to be run on.
#' @param data The data to analyse
#' @param analaysis_name The name for the analysis - will be prefixed to all
#' files generated.
#' @param analysis_path The file path to put the resultant files in - defaults
#' to working directory.
#' @param control List of options for BayesTraits.
#' @name makeInfile
#' @export

makeInfile <- function(tree, data, analysis_name, analysis_path, 
  control = list()) {

  control <- list(
    model = "",
    analysis_type = "mcmc",
    varrates = FALSE,
    rjtransform = "none",
    burnin = 10000000,
    iterations = 110000000,
    sample = 20000,
    stones = c(100, 1000),
    anc_states = FALSE,
    anc_nodes = "all"
  )

  # read and ladderize the tree.
  if (class(tree) != phylo) {
    tree <- ape::read.nexus(tree)
  }
  tree <- ape::ladderize(tree)

  # match data to tree - drop missing tips/data and order data to match the
  # tree.

  # translate options into BT format.
  if (model == "cont_random") {

  } else if (model == "cont_directional") {

  } else if (model ==  "cont_regression") {

  } else if (model == "ic") {

  } else if (model == "ic_correl") {

  } else if (model == "ic_regression") {

  } else if (model == "multistate") {

  } else if (model == "discrete_ind") {

  } else if (model == "discrete_dep") {

  }

  # if model is unspecified, then have a guess at what it ought to be.
  # If the data are discrete, then discrete (i.e., if only two states).
  # If there are more than one discrete states default to independent.
  # If the data appears multistate (hard to guess?) then assume multistate.
  # If continuous and varrates is ON assume IC.
  # If continuous and varrates is OFF assume continuous random walk.
  # Print a message saying which assumption has been made.
  if (control$analysis_type == "ml") {
    control$analysis_type <- 1
  } else if (control$analysis_type == "mcmc") {
    control$analysis_type <- 2
  }
}