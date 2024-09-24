# ********************************
# Purpose: Meta Analysis
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ********************************


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
path.outputs <- file.path(path.root, "outputs/Individual/single-network-analysis")
path.assoc_mat <- file.path(path.outputs , "association_matrices")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)

source("scripts/functions.R")
source("scripts/single-network-analysis/filtering.R")
source("scripts/single-network-analysis/fit_network.R")
source("scripts/single-network-analysis/assoc_mat.R")
source("scripts/single-network-analysis/net_properties.R")

matrices.gl <- lapply(datasets_names, function(x) {
  load(file=file.path(path.assoc_mat, agg_level, paste0("AssocMat_",x,".RData")))
  objects <- ls()
  return(get(objects[1]))
  })

matrices.mb <- lapply(datasets_names, function(x) {
  load(file=file.path(path.assoc_mat, agg_level, paste0("AssocMat_",x,".RData")))
  objects <- ls()
  return(get(objects[2]))
})

matrices.slr <- lapply(datasets_names, function(x) {
  load(file=file.path(path.assoc_mat, agg_level, paste0("AssocMat_",x,".RData")))
  objects <- ls()
  return(get(objects[3]))
})



meta.gl <- meta_analysis(matrices.gl, "mean")
props_asso.gl <- network_construct(assoMat.gl)






