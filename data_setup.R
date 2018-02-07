## setting up data for the documentation and examples.

# Posterior from a normal analysis (i.e. no varrates).

setwd("/home/hfg/pCloudDrive/BTprocessR_testdata/testdata")

artiodactyl_tree <- ape::read.nexus("Artiodactyl.trees")
devtools::use_data(artiodactyl_tree, artiodactyl_tree)

artiodactyl_traits <- read.table("Artiodactyl.txt")
devtools::use_data(artiodactyl_traits, artiodactyl_traits)
artiodactyl_post <- readLines("Artiodactyl.txt.Log.txt")
