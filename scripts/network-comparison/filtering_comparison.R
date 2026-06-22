# ***********************************
# Purpose: Filtering Phyloseq Objects
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ***********************************

filtering <- function(physeq_name, agg_level,
                      filtTax = "none", filtTaxPar = NULL,
                      filtSamp = "none", filtSampPar = NULL) {
  
  load(file.path(path.phylobj_sep, agg_level, paste0("HostDisease_", physeq_name, ".RData")))
  
  # 1. Filtering samples
  if (!identical(filtSamp, "none") && !is.null(filtSampPar$totalReads_quant)) {
    physeq_IBS <- filter_samples_totalReads(physeq_IBS, filtSampPar$totalReads_quant)$physeq
    physeq_H   <- filter_samples_totalReads(physeq_H,   filtSampPar$totalReads_quant)$physeq
  }
  
  # 2. Filtering taxa
  if (!identical(filtTax, "none") && !is.null(filtTaxPar$numbSamp_prop)) {
    physeq_IBS <- filter_taxa_numbSamp(physeq_IBS, filtTaxPar$numbSamp_prop)$physeq
    physeq_H   <- filter_taxa_numbSamp(physeq_H,   filtTaxPar$numbSamp_prop)$physeq
  }
  
  # 3. Keeping common taxa
  common_taxa <- intersect(
    phyloseq::taxa_names(physeq_IBS),
    phyloseq::taxa_names(physeq_H)
  )
  
  physeq_IBS <- phyloseq::prune_taxa(common_taxa, physeq_IBS)
  physeq_H   <- phyloseq::prune_taxa(common_taxa, physeq_H)
  
  # 4. Using NetCoMi to align matrices
  net_IBS <- NetCoMi::netConstruct(
    data = physeq_IBS,
    data2 = physeq_H,
    filtTax = "none",
    filtSamp = "none",
    measure = NULL,
    normMethod = "none",
    zeroMethod = "none",
    sparsMethod = "none",
    dissFunc = "signed",
    verbose = 3
  )
  
  if (!dir.exists(file.path(path.filt_phy, agg_level))) {
    dir.create(file.path(path.filt_phy, agg_level), recursive = TRUE)
  }
  
  physeq_filt_IBS <- net_IBS$countMat1
  physeq_filt_H   <- net_IBS$countMat2
  
  save(
    physeq_filt_IBS,
    physeq_filt_H,
    common_taxa,
    file = file.path(path.filt_phy, agg_level, paste0("Filt_", physeq_name, ".RData"))
  )
}
