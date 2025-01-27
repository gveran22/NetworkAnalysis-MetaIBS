# **********************************
# Purpose: Fitting Networks
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# **********************************

source("tools/functions.R")

fit_network <- function(physeq_name, agg_level){
  
  # Read the phyloseq object
  load(file.path(path.filt_phy, agg_level, 
                 paste0("Filt_",physeq_name,".RData")))
  
  pargs <- list(thresh=0.1, rep.num=60, seed=10010, ncores=20)
  
  se.gl <- spiec.easi(physeq_filt, method='glasso', nlambda=100,
                      lambda.min.ratio=1e-3, #lambda.log=FALSE,
                      pulsar.select=TRUE, pulsar.params=pargs)
  
  se.mb <- spiec.easi(physeq_filt, method='mb', nlambda=100,
                      lambda.min.ratio=1e-3, #lambda.log=FALSE,
                      pulsar.select=TRUE, pulsar.params=pargs)
  
  ranks <- rank.slr(physeq_filt)
  
  se.slr <- spiec.easi(physeq_filt , method='slr', nlambda=100,
                       lambda.min.ratio=1e-3, r=ranks, 
                       pulsar.select=TRUE, pulsar.params=pargs)
  
  se.slr$ebic <- sapply(se.slr, function(x)
    ebic(x$refit$stars, x$est$data,
         x$est$loglik[x$select$stars$opt.index]))

  if (!dir.exists(file.path(path.spiec_easi, agg_level))) {
    dir.create(file.path(path.spiec_easi, agg_level), recursive = TRUE)
  }
  save(se.gl, se.mb, se.slr, file=file.path(path.spiec_easi, agg_level, 
                                            paste0("NetFits_",physeq_name,".RData")))

}
