# *********************************
# Purpose: Preprocessing Comparison
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# *********************************

preprocessing_comparison <- function(physeq_name, agg_level){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  # Rename taxa using taxonomy name at agg_level
  tax_names <- as.character(
    as.data.frame(phyloseq::tax_table(physeq))[[agg_level]]
  )
  
  phyloseq::taxa_names(physeq) <- tax_names
  
  # Split groups
  physeq_IBS <- phyloseq::subset_samples(physeq, host_disease == "IBS")
  physeq_H <- phyloseq::subset_samples(physeq, host_disease == "Healthy")
  
  # Remove taxa absent within each group
  physeq_IBS <- prune_taxa(taxa_sums(physeq_IBS)>0, physeq_IBS)
  physeq_H <- prune_taxa(taxa_sums(physeq_H)>0, physeq_H)
  
  if (!dir.exists(file.path(path.phylobj_sep, agg_level))) {
    dir.create(file.path(path.phylobj_sep, agg_level), recursive = TRUE)
  }

  #Saving phyloseq objects
  save(physeq_IBS, physeq_H, file=file.path(path.phylobj_sep, agg_level, 
                                            paste0("HostDisease_",physeq_name,".RData")))
}
