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
path.datasets    <- file.path(path.root, "data/Individual/phyloseq_without_tree")
path.phylobj    <- file.path(path.root, "data/Agglomeration/Individual")
path.outputs <- file.path(path.root, "outputs/Individual/network-comparison")
path.phylobj_sep <- file.path(path.outputs , "phyloseq_IBS")
path.filt_phy <- file.path(path.outputs , "filtered-otus")
path.spiec_easi <- file.path(path.outputs , "spiec-easi-results")
path.assoc_mat <- file.path(path.outputs , "association_matrices")
path.properties <- file.path(path.outputs , "network_properties")
path.ranks <- file.path(path.outputs ,"ranks")
path.plots <- file.path(path.outputs , "plots")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)

source("scripts/functions.R")
source("scripts/network-comparison/preprocessing_comparison.R")
source("scripts/network-comparison/filtering_comparison.R")
source("scripts/network-comparison/fit_network_comparison.R")
source("scripts/network-comparison/assoc_mat_comparison.R")
source("scripts/network-comparison/net_properties_comparison.R")


# Order
agg_level="Order"
filtTax = "numbSamp"
filtTaxPar = list(numbSamp=5)
filtSamp = "totalReads"
filtSampPar = list(totalReads=500)

for (i in 1:length(datasets_names)) {
  preprocessing_comparison(datasets_names[i], agg_level)
  filtering(datasets_names[i], agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
  fit_network_comparison(datasets_names[i], agg_level)
  get_assoc_matrix(datasets_names[i], agg_level)
  get_network_properties(datasets_names[i], agg_level)
}


# Family
agg_level="Family"
filtTax = "numbSamp"
filtTaxPar = list(numbSamp=5)
filtSamp = "totalReads"
filtSampPar = list(totalReads=500)

for (i in 1:length(datasets_names)) {
  preprocessing_comparison(datasets_names[i], agg_level)
  filtering(datasets_names[i], agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
  fit_network_comparison(datasets_names[i], agg_level)
  get_assoc_matrix(datasets_names[i], agg_level)
  get_network_properties(datasets_names[i], agg_level)
}


# Genus
agg_level="Genus"
filtTax = "highestVar"
filtTaxPar = list(highestVar=100)
filtSamp = "totalReads"
filtSampPar = list(totalReads=500)

for (i in 1:length(datasets_names)) {
  preprocessing_comparison(datasets_names[i], agg_level)
  filtering(datasets_names[i], agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
  fit_network_comparison(datasets_names[i], agg_level)
  get_assoc_matrix(datasets_names[i], agg_level)
  get_network_properties(datasets_names[i], agg_level)
}

