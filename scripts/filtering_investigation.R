# ********************************
# Purpose: Filtering Investigation
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

library(phyloseq) 
library(tidyverse)
library(ggplot2)
library(SpiecEasi)
library(igraph)
library(pulsar)
library(microViz)
library(NetCoMi)

path.root <- "~/MetaIBS"
path.datasets    <- file.path(path.root, "data/phyloseq_without_tree")
path.phylobj    <- file.path(path.root, "build/Agglomeration/Individual")
datasets        <- list.files(path.datasets, pattern=".rds")
datasets_names  <- sub(".*_(.*)\\..*", "\\1", datasets)
path.outputs <- file.path(path.root, "outputs/investigation")

# ****************
# 2. FILTERING
# ****************

#################################################################
################### Single Analysis #############################
#################################################################

filtering_netcomi <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.datasets, paste0("physeq_",physeq_name,".rds")))
  
  dimensions <- vector("integer")
  dimensions <- dim(physeq@otu_table)
  
  physeq_agglo <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  dimensions <- c(dimensions, dim(physeq_agglo@otu_table))
  
  net_asso <- netConstruct(physeq_agglo,
                           filtTax =  filtTax,
                           filtTaxPar = filtTaxPar,
                           filtSamp =   filtSamp,
                           filtSampPar =  filtSampPar,
                           measure = NULL,
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3)
  
  dimensions <- c(dimensions, dim(net_asso$countMat1))
  names(dimensions) <- c("samples", "ASVs", "samples_agglo", "ASVs_agglo",
                         "samples_filtered", "ASVs_filtered")
  
  return(dimensions)
}


all_filtering_netcomi <- function(physeq, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar){
  
  prs <- expand_grid(physeq, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
  
  dimensions <- vector("integer")
  dimensions_final <- vector("integer")
  
  for (i in 1:nrow(prs)) {
    dimensions <- filtering_netcomi(physeq_name = prs$physeq[i],
                            agg_level = prs$agg_level[i],
                            filtTax =  prs$filtTax[i],
                            filtTaxPar = prs$filtTaxPar[i],
                            filtSamp =   prs$filtSamp[i],
                            filtSampPar = prs$filtSampPar[i])
    
    dimensions_final <- rbind(dimensions_final,dimensions, deparse.level = 0)
    
  }

  prs$filtTaxPar <- sapply(prs$filtTaxPar, function(x) x[[1]])
  prs$filtSampPar <- sapply(prs$filtSampPar, function(x) x[[1]])
    
  
  prs_final <- cbind(prs,dimensions_final)
  #dimensions_final <- data.frame(dimensions_final)
  return(prs_final)
}

##### Filtering: numbSamp (5,10) - totalReads (0,500,1000)
simulation <- all_filtering_netcomi(physeq=datasets_names,
                            agg_level=c("Class","Phylum", "Order","Family", "Genus"),
                            filtTax =  "numbSamp",
                            filtTaxPar = c(list(numbSamp=5),list(numbSamp=10)),
                            filtSamp =   "totalReads",
                            filtSampPar = c(list(totalReads=0),list(totalReads=500),list(totalReads=1000)))

write.csv(simulation, file.path(path.outputs, "single_numbSamp_totalReads.csv"), 
          row.names = FALSE)


##### Filtering: highestVar (40,100,150) - totalReads (0,500,1000)
simulation <- all_filtering_netcomi(physeq=datasets_names,
                                    agg_level=c("Class","Phylum", "Order","Family", "Genus"),
                                    filtTax =  "highestVar",
                                    filtTaxPar = c(list(highestVar=40), list(highestVar=100), list(highestVar=150)),
                                    filtSamp =   "totalReads",
                                    filtSampPar = c(list(totalReads=0),list(totalReads=500),list(totalReads=1000)))

write.csv(simulation, file.path(path.outputs, "single_highestVar_totalReads.csv"), 
          row.names = FALSE)

#################################################################
################### Comparison Networks #########################
#################################################################

filtering_comparison <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.datasets, paste0("physeq_",physeq_name,".rds")))
  
  dimensions <- vector("integer")
  dimensions <- dim(physeq@otu_table)
  
  # Read the agglomerated phyloseq object
  phyloseq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  dimensions <- vector("integer")
  dimensions <- dim(phyloseq@otu_table)
  
  # Renaming ASVs for agg_level name taxonomy
  otutab <- as.data.frame(otu_table(phyloseq))
  
  if(identical(colnames(otutab), taxa_names(phyloseq))){
    colnames(otutab) <- as.data.frame(tax_table(phyloseq))[[agg_level]]
  }else{
    print("ASVs in OTU Table does not match the ASVs in the Taxonomy Table")
  }
  
  otutab <- otu_table(as.matrix(otutab), taxa_are_rows=FALSE)
  phyloseq@otu_table <- otutab
  taxa <- as.data.frame(phyloseq@tax_table)
  
  rownames(phyloseq@tax_table) <- taxa[,agg_level]
  
  variable_names <- c("IBS", "Healthy")
  
  phyloseq_IBS <- phyloseq::subset_samples(phyloseq, host_disease == "IBS")
  phyloseq_IBS <- prune_taxa(taxa_sums(phyloseq_IBS)>0, phyloseq_IBS)
  
  phyloseq_H <- phyloseq::subset_samples(phyloseq, host_disease == "Healthy")
  phyloseq_H <- prune_taxa(taxa_sums(phyloseq_H)>0, phyloseq_H)
  
  dimensions <- c(dimensions, dim(phyloseq_IBS@otu_table), dim(phyloseq_H@otu_table))
  
  net_IBS <- netConstruct(data= phyloseq_IBS,
                          data2= phyloseq_H,
                          filtTax = filtTax,
                          filtTaxPar = filtTaxPar, 
                          filtSamp = filtSamp,
                          filtSampPar = filtSampPar, 
                          measure = NULL,
                          normMethod = "none", 
                          zeroMethod = "none",
                          sparsMethod = "none", 
                          dissFunc = "signed",
                          jointPrepro = FALSE,
                          verbose = 3)
  
  dimensions <- c(dimensions, dim(net_IBS$countMat1), dim(net_IBS$countMat2))
  names(dimensions) <- c("samples_agglo", "ASVs_agglo", "samples_IBS", "ASVs_IBS", "samples_H", "ASVs_H",
                        "samples_filtered_IBS", "ASVs_filtered_IBS", "samples_filtered_H", "ASVs_filtered_H")
  return(dimensions)
}

all_filtering_comparison <- function(physeq, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar){
  
  prs <- expand_grid(physeq, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar)
  
  dimensions <- vector("integer")
  dimensions_final <- vector("integer")
  
  for (i in 1:nrow(prs)) {
    dimensions <- filtering_comparison(physeq_name = prs$physeq[i],
                                    agg_level = prs$agg_level[i],
                                    filtTax =  prs$filtTax[i],
                                    filtTaxPar = prs$filtTaxPar[i],
                                    filtSamp =   prs$filtSamp[i],
                                    filtSampPar = prs$filtSampPar[i])
    
    dimensions_final <- rbind(dimensions_final,dimensions, deparse.level = 0)
    
  }
  
  #prs <- as.data.frame(prs)
  
  prs$filtTaxPar <- sapply(prs$filtTaxPar, function(x) x[[1]])
  prs$filtSampPar <- sapply(prs$filtSampPar, function(x) x[[1]])
  
  
  prs_final <- cbind(prs,dimensions_final)
  #dimensions_final <- data.frame(dimensions_final)
  return(prs_final)
}

##### Filtering: numbSamp (5,10) - totalReads (0,500)
simulation <- all_filtering_comparison(physeq=datasets_names,
                                agg_level=c("Class","Phylum", "Order","Family", "Genus"),
                                filtTax =  "numbSamp",
                                filtTaxPar = c(list(numbSamp=5),list(numbSamp=10)),
                                filtSamp =   "totalReads",
                                filtSampPar = c(list(totalReads=0),list(totalReads=500)))

write.csv(simulation, file.path(path.outputs, "comparison_numbSamp-totalReads.csv"), 
          row.names = FALSE)


##### Filtering: highestVar (40,100,150) - totalReads (0,500)
simulation <- all_filtering_comparison(physeq=datasets_names,
                                    agg_level=c("Class","Phylum", "Order","Family", "Genus"),
                                    filtTax =  "highestVar",
                                    filtTaxPar = c(list(highestVar=40), list(highestVar=100), list(highestVar=150)),
                                    filtSamp =   "totalReads",
                                    filtSampPar = c(list(totalReads=0),list(totalReads=500)))

write.csv(simulation, file.path(path.outputs, "comparison_highestVar-totalReads.csv"), 
          row.names = FALSE)
