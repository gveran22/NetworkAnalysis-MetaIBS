# *************************************
# Purpose: Venn-Diagramms
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

# ROOT DIRECTORY (to modify on your computer)
path.root <- "~/MetaIBS"
path.datasets    <- file.path(path.root, "data/Individual")
path.outputs <- file.path(path.root, "outputs/single-network-analysis/Individual")
path.assoc_mat <- file.path(path.outputs , "association_matrices")
path.venn_diag <- file.path(path.outputs , "venn-diagramm")

source("scripts/functions.R")

datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)

agg_level <- c("Order", "Family", "Genus")


for (agg in agg_level){
  for (i in 1:length(datasets_names)) {
    
    load(file.path(path.assoc_mat, agg_level, 
                   paste0("AssocMat_",datasets_names[i],".RData")))
    
    if (!dir.exists(file.path(path.venn_diag, agg_level))) {
      dir.create(file.path(path.venn_diag, agg_level), recursive = TRUE)
    }
    get_venn_diagram(assoMat.gl, assoMat.mb, assoMat.slr, datasets_names[i])
    
  }
}
