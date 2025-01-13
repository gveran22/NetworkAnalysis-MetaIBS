# *************************************
# Purpose: Venn-Diagramms Comparison
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# *************************************

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

source("scripts/functions.R")

########################### Individual Analysis ##################################
# ROOT DIRECTORY (to modify on your computer)
path.root <- "~/MetaIBS"
path.datasets    <- file.path(path.root, "data/Individual")
path.outputs <- file.path(path.root, "outputs/network-comparison/Individual")
path.assoc_mat <- file.path(path.outputs , "association_matrices")
path.venn_diag <- file.path(path.outputs , "venn-diagramm")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)


agg_level <- c("Order", "Family","Genus")

for (i in 1:length(datasets_names)) {
  
  load(file.path(path.assoc_mat, agg_level, 
                 paste0("AssocMat_",physeq_name,".RData")))
  
  if (!dir.exists(file.path(path.venn_diag, agg_level))) {
    dir.create(file.path(path.venn_diag, agg_level), recursive = TRUE)
  }
  get_venn_diagram(assoMat_IBS.gl, assoMat_IBS.mb, assoMat_IBS.slr, paste0(datasets_names[i], "_IBS"))
  get_venn_diagram(assoMat_H.gl, assoMat_H.mb, assoMat_H.slr, paste0(datasets_names[i],  "_Healthy"))
  
}

########################### Combined Analysis ##################################

# Combined analysis paths
combined_paths <- list(
  "all" = "Combined/all",
  "variable_region" = "Combined/variable_region",
  "sample_type" = "Combined/sample_type",
  "sequencing_tech" = "Combined/sequencing_tech"
)

agg_level <- "Order"
# Run combined analysis for each variable and configuration
for (path_suffix in names(combined_paths)) {
  path.datasets    <- file.path(path.root, "data/Combined/", combined_paths[[path_suffix]])
  path.outputs <- file.path(path.root, "outputs/network-comparison", combined_paths[[path_suffix]])
  path.assoc_mat <- file.path(path.outputs , "association_matrices")
  path.venn_diag <- file.path(path.outputs , "venn-diagramm")
  
  datasets        <- list.files(path.datasets, pattern=".rds")
  datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)
  
  # Run the same analysis for each configuration
  for (i in 1:length(datasets_names)) {
    
    load(file.path(path.assoc_mat, agg_level, 
                   paste0("AssocMat_",physeq_name,".RData")))
    
    if (!dir.exists(file.path(path.venn_diag, agg_level))) {
      dir.create(file.path(path.venn_diag, agg_level), recursive = TRUE)
    }
    get_venn_diagram(assoMat_IBS.gl, assoMat_IBS.mb, assoMat_IBS.slr, paste0(datasets_names[i], "_IBS"))
    get_venn_diagram(assoMat_H.gl, assoMat_H.mb, assoMat_H.slr, paste0(datasets_names[i],  "_Healthy"))
    
  }
}

