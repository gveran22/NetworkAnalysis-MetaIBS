# *********************************
# Purpose: Preprocessing Comparison
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# *********************************

preprocessing_comparison <- function(physeq_name, agg_level){
  
  # Read the phyloseq object
  phyloseq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
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
  
  physeq_IBS <- phyloseq::subset_samples(phyloseq, host_disease == "IBS")
  physeq_IBS <- prune_taxa(taxa_sums(physeq_IBS)>0, physeq_IBS)
  
  physeq_H <- phyloseq::subset_samples(phyloseq, host_disease == "Healthy")
  physeq_H <- prune_taxa(taxa_sums(physeq_H)>0, physeq_H)
  
  if (!dir.exists(file.path(path.phylobj_sep, agg_level))) {
    dir.create(file.path(path.phylobj_sep, agg_level), recursive = TRUE)
  }

  #Saving phyloseq objects
  save(physeq_IBS, physeq_H, file=file.path(path.phylobj_sep, agg_level, 
                                            paste0("HostDisease_",physeq_name,".RData")))
  
}

