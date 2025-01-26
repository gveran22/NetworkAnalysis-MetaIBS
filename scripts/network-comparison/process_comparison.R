# ************************************
# Purpose: Network Comparison Analysis
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ************************************

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
source("scripts/functions.R")
source("scripts/network-comparison/preprocessing_comparison.R")
source("scripts/network-comparison/filtering_comparison.R")
source("scripts/network-comparison/fit_network_comparison.R")
source("scripts/network-comparison/assoc_mat_comparison.R")
source("scripts/network-comparison/net_properties_comparison.R")
source("scripts/analysis_configs.R")


run_analysis <- function(agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar, datasets_names) {
  for (i in seq_along(datasets_names)) {
    preprocessing_comparison(datasets_names[i], agg_level)
    filtering(datasets_names[i], agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
    fit_network_comparison(datasets_names[i], agg_level)
    get_assoc_matrix(datasets_names[i], agg_level)
    
    tryCatch({
      get_network_properties(datasets_names[i], agg_level)  # Call your netAnalyze function here
    }, error = function(e) {
      # Handle the empty network error
      warning(paste("Skipping dataset", datasets_names[i], "due to empty network."))
    })
  }
}

########################### Individual Analysis ##################################

path.datasets    <- file.path(path.root, "data/Individual")
path.phylobj    <- file.path(path.root, "data/Agglomeration/Individual")
path.outputs <- file.path(path.root, "outputs/network-comparison/Individual")
path.phylobj_sep <- file.path(path.outputs , "phyloseq_IBS")
path.filt_phy <- file.path(path.outputs , "filtered_otus")
path.spiec_easi <- file.path(path.outputs , "spiec-easi_results")
path.assoc_mat <- file.path(path.outputs , "association_matrices")
path.properties <- file.path(path.outputs , "network_properties")
path.plots <- file.path(path.outputs , "plots")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)

# Run analysis for individual datasets
for (config in analysis_configs) {
  run_analysis(config$agg_level, config$filtTax, config$filtTaxPar, config$filtSamp, config$filtSampPar, datasets_names)
}



########################### Combined Analysis ####################################

# Combined analysis paths
combined_paths <- list(
  "all" = "all",
  "variable_region" = "variable_region",
  "sample_type" = "sample_type",
  "sequencing_tech" = "sequencing_tech"
)

# Run combined analysis for each variable and configuration
for (path_suffix in names(combined_paths)) {
  path.datasets    <- file.path(path.root, "data/Combined", combined_paths[[path_suffix]])
  path.phylobj    <- file.path(path.root, "data/Agglomeration/Combined", combined_paths[[path_suffix]])
  path.outputs <- file.path(path.root, "outputs/network-comparison/Combined", combined_paths[[path_suffix]])
  path.phylobj_sep <- file.path(path.outputs , "phyloseq_IBS")
  path.filt_phy <- file.path(path.outputs , "filtered_otus")
  path.spiec_easi <- file.path(path.outputs , "spiec-easi_results")
  path.assoc_mat <- file.path(path.outputs , "association_matrices")
  path.properties <- file.path(path.outputs , "network_properties")
  path.plots <- file.path(path.outputs , "plots")
  
  datasets        <- list.files(path.datasets, pattern=".rds")
  datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)
  
  # Run the same analysis for each configuration
  for (config in analysis_configs) {
    run_analysis(config$agg_level, config$filtTax, config$filtTaxPar, config$filtSamp, config$filtSampPar, datasets_names)
  }
}

