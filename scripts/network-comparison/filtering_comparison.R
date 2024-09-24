# ***********************************
# Purpose: Filtering Phyloseq Objects
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ***********************************

filtering <- function(physeq_name, agg_level,
                      filtTax, filtTaxPar, 
                      filtSamp, filtSampPar){
  
  # Read the phyloseq object
  load(file.path(path.phylobj_sep, agg_level, paste0("HostDisease_",physeq_name,".RData")))
  
  # Filtering
  net_IBS <- netConstruct(data= physeq_IBS,
                          data2= physeq_H,
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
  
  physeq_filt_IBS <- net_IBS$countMat1
  physeq_filt_H <- net_IBS$countMat2
  
  #Saving filtered OTU
  save(physeq_filt_IBS, physeq_filt_H, file=file.path(path.filt_phy, agg_level, paste0("Filt_",physeq_name,".RData")))
}








