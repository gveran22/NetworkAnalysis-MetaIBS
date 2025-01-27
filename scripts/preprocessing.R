# **********************************
# Purpose: Merging & Separation
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
path.phylobj    <- file.path(path.root, "build/Combined")
path.output <- file.path(path.root, "build/Agglomeration")

# ****************
# 2. PROCESSING ##
# ****************

############# 2.1. Read merge phyloseq object ########################################
physeq.all <- readRDS(file.path(path.phylobj, "physeq_all.rds")) # 2,584 samples 76,218 taxa
cat("Total Nb of samples:", nsamples(physeq.all))

# Defining taxonomy level to merge
agg_level <- "Phylum"

phyloseq <- tax_glom(physeq.all, taxrank = agg_level, NArm = FALSE)

# Renaming missing 'agg_level'
phyloseq <- renameTaxa(phyloseq, unclass = NULL,
                       unknown = c(NA, " ", "unclassified"), 
                       pat = "<name>", 
                       substPat = "<name>_<subst_name>(<subst_R>)",
                       numDupli = agg_level)

# Create directories if they don't exist
if (!dir.exists(file.path(path.output, "Combined/all", agg_level))) {
  dir.create(file.path(path.output, "Combined/all", agg_level), recursive = TRUE)
}

saveRDS(phyloseq, file.path(path.output, "Combined/all", agg_level,
                            paste0("agglo_all.rds")))


############# 2.1.1. Dividing datasets #################################

if (!dir.exists(file.path(path.output, "Individual", agg_level))) {
  dir.create(file.path(path.output, "Individual", agg_level), recursive = TRUE)
}

physeq <- subset_samples(phyloseq, author == "AGP")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","agp",".rds")))

physeq <- subset_samples(phyloseq, author == "Fukui")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","fukui",".rds")))

physeq <- subset_samples(phyloseq, author == "Hugerth")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","hugerth",".rds")))

physeq <- subset_samples(phyloseq, author == "Labus")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","labus",".rds")))

physeq <- subset_samples(phyloseq, author == "Liu")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","liu",".rds")))

physeq <- subset_samples(phyloseq, author == "LoPresti")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","lopresti",".rds")))

physeq <- subset_samples(phyloseq, author == "Mars")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","mars",".rds")))

physeq <- subset_samples(phyloseq, author == "Nagel")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","nagel",".rds")))

physeq <- subset_samples(phyloseq, author == "Pozuelo")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","pozuelo",".rds")))

physeq <- subset_samples(phyloseq, author == "Zeber-Lubecka")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","zeber",".rds")))

physeq <- subset_samples(phyloseq, author == "Zhu")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","zhu",".rds")))

physeq <- subset_samples(phyloseq, author == "Zhuang")
physeq <- prune_taxa(taxa_sums(physeq)>0, physeq)
cat("Total Nb of samples:", nsamples(physeq))
saveRDS(physeq, file.path(path.output,  "Individual", agg_level, 
                          paste0("agglo_","zhuang",".rds")))


############# 2.2. Agglomeration depending on variable region #################################

###### Separate variable region V4 (agp - mars - nagel -pozuelo - zhu) #####

if (!dir.exists(file.path(path.output, "Combined/variable_region", agg_level))) {
  dir.create(file.path(path.output, "Combined/variable_region", agg_level), recursive = TRUE)
}

physeq.v4 <- subset_samples(phyloseq, variable_region == 'V4') # 1,592 samples
physeq.v4 <- prune_taxa(taxa_sums(physeq.v4)>0, physeq.v4) # remove ASVs that are not present anymore
cat("Nb of V4 samples:", nsamples(physeq.v4))

saveRDS(physeq.v4, file.path(path.output, "Combined/variable_region", agg_level, "agglo_v4.rds"))


############# 2.3. Agglomeration depending on sample_type ##########################

if (!dir.exists(file.path(path.output, "Combined/sample_type", agg_level))) {
  dir.create(file.path(path.output, "Combined/sample_type", agg_level), recursive = TRUE)
}

physeq.fecal <- subset_samples(phyloseq, sample_type == 'stool') # 2,153 samples
physeq.fecal <- prune_taxa(taxa_sums(physeq.fecal)>0, physeq.fecal) # remove ASVs that are not present anymore
cat("Nb of fecal samples:", nsamples(physeq.fecal))

physeq.sigmoid <- subset_samples(phyloseq, sample_type == 'sigmoid') # 431 samples
physeq.sigmoid <- prune_taxa(taxa_sums(physeq.sigmoid)>0, physeq.sigmoid) # remove ASVs that are not present anymore
cat("Nb of sigmoid samples:", nsamples(physeq.sigmoid))


saveRDS(physeq.fecal, file.path(path.output, "Combined/sample_type", agg_level, "agglo_fecal.rds"))
saveRDS(physeq.sigmoid, file.path(path.output, "Combined/sample_type", agg_level, "agglo_sigmoid.rds"))


############# 2.4. Agglomeration depending on sequencing_tech #####################

if (!dir.exists(file.path(path.output, "Combined/sequencing_tech", agg_level))) {
  dir.create(file.path(path.output, "Combined/sequencing_tech", agg_level), recursive = TRUE)
}

physeq.pyroseq <- subset_samples(phyloseq, sequencing_tech == '454 pyrosequencing') # 109 samples
physeq.pyroseq <- prune_taxa(taxa_sums(physeq.pyroseq)>0, physeq.pyroseq) # remove ASVs that are not present anymore
cat("Nb of 454 pyrosequencing samples:", nsamples(physeq.pyroseq))

physeq.ion <- subset_samples(phyloseq, sequencing_tech == 'Ion Torrent') # 120 samples
physeq.ion <- prune_taxa(taxa_sums(physeq.ion)>0, physeq.ion) # remove ASVs that are not present anymore
cat("Nb of Ion Torrent samples:", nsamples(physeq.ion))

physeq.illumina <- subset_samples(phyloseq, grepl("^Illumina", sequencing_tech)) #  2,355 samples
physeq.illumina <- prune_taxa(taxa_sums(physeq.illumina)>0, physeq.illumina) # remove ASVs that are not present anymore
cat("Nb of Illumina samples:", nsamples(physeq.illumina))

saveRDS(physeq.pyroseq, file.path(path.output, "Combined/sequencing_tech", agg_level, "agglo_pyroseq.rds"))
saveRDS(physeq.ion, file.path(path.output, "Combined/sequencing_tech", agg_level, "agglo_ion.rds"))
saveRDS(physeq.illumina, file.path(path.output, "Combined/sequencing_tech", agg_level, "agglo_illumina.rds"))

