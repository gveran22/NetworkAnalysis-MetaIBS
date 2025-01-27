# **********************************
# Purpose: Single Network Analysis
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# **********************************

myPaths <- .libPaths()
myPaths <- c(myPaths, "~/MetaIBS/MetaIBS-library")
myPaths <- c(myPaths[3], myPaths[1], myPaths[2])
.libPaths(myPaths)  # add new path

# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries

library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
library(tidyverse)
library(ggplot2)
library(SpiecEasi)
library(igraph)
library(microViz)
library(NetCoMi)


# ROOT DIRECTORY (to modify on your computer)
path.root <- "~/MetaIBS"
path.datasets    <- file.path(path.root, "data/phyloseq_without_tree")
path.phylobj    <- file.path(path.root, "build/Agglomeration/Individual")
path.outputs <- file.path(path.root, "outputs/single-network-analysis/Individual")
path.filt_phy <- file.path(path.outputs , "filtered_otus")
path.spiec_easi <- file.path(path.outputs , "spiec-easi_results")
path.assoc_mat <- file.path(path.outputs , "association_matrices")
path.properties <- file.path(path.outputs , "network_properties")
path.ranks <- file.path(path.outputs ,"ranks")
path.plots <- file.path(path.outputs , "plots")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)

source("tools/functions.R")
source("scripts/single-network-analysis/filtering.R")
source("scripts/single-network-analysis/fit_network.R")
source("scripts/single-network-analysis/assoc_mat.R")
source("scripts/single-network-analysis/net_properties.R")
source("tools/analysis_configs.R")


run_single_network_analysis <- function(agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar, datasets_names) {
  for (i in seq_along(datasets_names)) {
    filtering(datasets_names[i], agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
    fit_network(datasets_names[i], agg_level)
    get_assoc_matrix(datasets_names[i], agg_level)
    get_network_properties(datasets_names[i], agg_level)
  }
}

################ Individual Analysis #######################

for (config in analysis_configs) {
  run_single_network_analysis(config$agg_level, config$filtTax, config$filtTaxPar, config$filtSamp, config$filtSampPar, datasets_names)
}
