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

makeInfile <- function(tree, data, species_col = NULL, analysis_name, 
  analysis_path, job = FALSE, bt.control = list(), job.control = list()) {

  opts <- list(
    model = "",
    analysis_type = "mcmc",
    varrates = FALSE,
    rjtransform = "none",
    burnin = 10000000,
    iterations = 110000000,
    sample = 20000,
    stones = c(1000, 5000),
    ancstates = FALSE,
    bayestraits = "BayesTraitsV3",
    shellscript = TRUE
  )
  opts[names(bt.control)] <- bt.control

  if (job) {
    job.opts <- list(
      walltime = "4:0:0",
      mem = 16,
      tmpfs = 1,
      jobname = analysis_name
    )
    job.opts[names(job.control)] <- job.control
    if (!"destination" %in% names(job.opts)) {
      stop("Please provide a destination for analysis outputs in job.control.")
    }
    if (!"wd" %in% names(job.opts)) {
      stop("If making a job script please specify the working directory for the
        job in the job.control list.")
    }
  }

  if (class(data) == "matrix") {
    data <- as.data.frame(data)
  }

  # create the directory to put the analysis files in, if needed.
  dir.create(file.path(analysis_path), showWarnings = FALSE)

  # read and ladderize the tree.
  if (!class(tree) %in% c("phylo", "multiPhylo")) {
    tree <- ape::read.nexus(tree)
  }
  if (class(tree) == "phylo") {
    tree <- ape::ladderize(tree)
  } else if (class(tree) == "multiPhylo") {
    tree <- lapply(tree, ape::ladderize)
    class(tree) <- "multiPhylo"
  }

  # match data to tree - drop missing tips/data and order data to match the
  # tree.
  if (is.null(rownames(data)) && is.null(species_col)) {
    stop("If species names are not the rownames of data then species_col must
      be provided")
  }

  # is species_col is not provided then assume that the rownames are the species
  # names.
  if (is.null(species_col)) {
    data$species <- rownames(data)
    species_col <- "species"
  }

  tr_n_d <- tree$tip.label[!tree$tip.label %in% data[, species_col]]
  d_n_tr <- data[, species_col] %in% tree$tip.label
  tree <- ape::drop.tip(tree, tip = tr_n_d)
  data <- data[d_n_tr, ]
  data <- cbind(data[, species_col], data[, colnames(data) != species_col])

  # translate options into BT format.
  if (opts$model == "cont_random") {
    opts$ model <- 4
  } else if (opts$model == "cont_directional") {
    opts$ model <- 5
  } else if (opts$model ==  "cont_regression") {
    opts$ model <- 6
  } else if (opts$model == "ic") {
    opts$ model <- 7
  } else if (opts$model == "ic_correl") {
    opts$ model <- 8
  } else if (opts$model == "ic_regression") {
    opts$ model <- 9
  } else if (opts$model == "multistate") {
    opts$ model <- 1
  } else if (opts$model == "discrete_ind") {
    opts$ model <- 2
  } else if (opts$model == "discrete_dep") {
    opts$ model <- 3
  } else if (model == "") {
    stop("No model specified!")
  } else {
    stop(paste("Model", opts$model, "not implemented in makeInfile. Please
      specify model and setup analysis manually."))
  }

  if (opts$analysis_type == "ml") {
    opts$analysis_type <- 1
    mcmc <- FALSE
  } else if (opts$analysis_type == "mcmc") {
    opts$analysis_type <- 2
  }

  if (opts$varrates && opts$model %in% c(4, 5, 6)) {
    stop("Variable rates is not available for continuous models.")
  }
  if (opts$ancstates && opts$model %in% c(7, 8, 9)) {
    stop("Ancestral state reconstruction is not available for independent 
      contrast models.")
  }

  if (opts$ancstates) {
    tags <- list()
    commands <- list()
    descnodes <- unique(tree$edge[, 1])
    for (i in seq_along(descnodes)) {
      tips <- getTipNames(tree, descnodes[i])
      tags[[i]] <- paste0(
        "AddTag ",
        "node", descnodes[i], " ",
        paste(tips, collapse = " ")
      )
      commands[[i]] <- paste0("AddMRCA node", descnodes[i], " node", descnodes[i])
    }
  }

  # write data
  dfile <- paste0(analysis_name, ".txt")
  tfile <- paste0(analysis_name, ".trees")
  ifile <- paste0(analysis_name, ".in")
  if (opts$shellscript) {
    afile <- paste0(analysis_name, "_analysis.sh")
  } else {
    afile <- paste0(analysis_name, "_command.txt")  
  }
  jfile <- paste0(analysis_name, "_job.sh")

  write.table(data, file = file.path(analysis_path, dfile), col.names = FALSE,
    row.names = FALSE, quote = FALSE)
  
  ape::write.nexus(tree, file = file.path(analysis_path, tfile))
  
  sink(file.path(analysis_path, ifile))
    cat(paste0(opts$model, "\n"))
    cat(paste0(opts$analysis_type, "\n"))
    if (opts$rjtransform != "none") {
      cat(paste0("rjlocaltransform ", opts$rjtransform, "\n"))
    }
    if (opts$varrates) {
      cat("varrates\n")
    }
    cat(paste0("burnin ", opts$burnin, "\n"))
    cat(paste0("iterations ", opts$iterations, "\n"))
    cat(paste0("sample ", opts$sample, "\n"))
    cat(paste0("stones ", opts$stones, "\n"))
    if (opts$ancstates) {
      ilapply(tags, function(x) cat(paste0(x, "\n")))
      ilapply(commands, function(x) cat(paste0(x, "\n")))
    }
    cat("run\n")
  sink()

  # write command
  sink(file.path(analysis_path, afile))
    if (opts$shellscript) {
      cat("#!/bin/bash\n")
    }
    cat(paste0("./", opts$bayestraits, " ", tfile, " ", dfile, 
      " ", "<", " ", ifile, "\n"))
  sink()

  if (job) {
    sink(file.path(analysis_path, jfile))
      cat("#!/bin/bash -l\n")
      cat("#$ -S /bin/bash\n")
      cat(paste0("#$ -l h_rt=", job.opts$walltime, "\n"))
      cat(paste0("#$ -l mem=", job.opts$mem, "G\n"))
      cat(paste0("#$ -l tmpfs=", job.opts$tmpfs, "G\n"))
      cat(paste0("#$ -N ", job.opts$jobname, "\n"))
      cat(paste0("#$ -wd ", job.opts$wd, "\n"))
      cat("module unload compilers mpi\n")
      cat("module load r/recommended\n")
      cat(paste("cp", tfile, dfile, ifile, opts$bayestraits, "$TMPDIR\n"))
      cat(paste(opts$bayestraits, tfile, dfile, "<", ifile, "> stdout\n"))
      cat(
        paste0("tar zcvf ", job.opts$destination, "/", 
          analysis_name, "_$JOB_ID.tar.gz $TMPDIR\n")
        )
    sink()
  }
}
