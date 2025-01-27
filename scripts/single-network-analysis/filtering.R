# ***********************************
# Purpose: Filtering Phyloseq Objects
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ***********************************

filtering <- function(physeq_name, agg_level,
                      filtTax, filtTaxPar, 
                      filtSamp, filtSampPar){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  # Filtering
  net_asso <- netConstruct(physeq,
                           taxRank = agg_level,
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
  
  
  if (!dir.exists(file.path(path.filt_phy, agg_level))) {
    dir.create(file.path(path.filt_phy, agg_level), recursive = TRUE)
  }
  
  physeq_filt <- net_asso$countMat1
  #Saving filtered OTU
  save(physeq_filt, file=file.path(path.filt_phy, agg_level, paste0("Filt_",physeq_name,".RData")))

}
