# **********************************
# Purpose: Fitting Networks
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# **********************************

source("tools/functions.R")

fit_network_comparison <- function(physeq_name, agg_level){
  
  # Read the phyloseq object
  load(file.path(path.filt_phy, agg_level, 
                 paste0("Filt_",physeq_name,".RData")))
  
  pargs <- list(thresh=0.1, rep.num=60, seed=10010, ncores=20)
  
  se_IBS.gl <- spiec.easi(physeq_filt_IBS, method='glasso', nlambda=100,
                      lambda.min.ratio=1e-3, #lambda.log=FALSE,
                      pulsar.select=TRUE, pulsar.params=pargs)
  se_H.gl <- spiec.easi(physeq_filt_H, method='glasso', nlambda=100,
                      lambda.min.ratio=1e-3, #lambda.log=FALSE,
                      pulsar.select=TRUE, pulsar.params=pargs)
  
  se_IBS.mb <- spiec.easi(physeq_filt_IBS, method='mb', nlambda=100,
                      lambda.min.ratio=1e-3, #lambda.log=FALSE,
                      pulsar.select=TRUE, pulsar.params=pargs)
  se_H.mb <- spiec.easi(physeq_filt_H, method='mb', nlambda=100,
                          lambda.min.ratio=1e-3, #lambda.log=FALSE,
                          pulsar.select=TRUE, pulsar.params=pargs)
  
  
  ranks_IBS <- rank.slr(physeq_filt_IBS)
  ranks_H <- rank.slr(physeq_filt_H)
  
  se_IBS.slr <- spiec.easi(physeq_filt_IBS, method='slr', nlambda=100,
                       lambda.min.ratio=1e-3, r=ranks_IBS, #lambda.log=FALSE,
                       pulsar.select=TRUE, pulsar.params=pargs)
  se_H.slr <- spiec.easi(physeq_filt_H, method='slr', nlambda=100,
                       lambda.min.ratio=1e-3, r=ranks_H, #lambda.log=FALSE,
                       pulsar.select=TRUE, pulsar.params=pargs)
  
  se_IBS.slr$ebic <- sapply(se_IBS.slr, function(x)
    ebic(x$refit$stars, x$est$data,
         x$est$loglik[x$select$stars$opt.index]))
  se_H.slr$ebic <- sapply(se_H.slr, function(x)
    ebic(x$refit$stars, x$est$data,
         x$est$loglik[x$select$stars$opt.index]))
  
  if (!dir.exists(file.path(path.spiec_easi, agg_level))) {
    dir.create(file.path(path.spiec_easi, agg_level), recursive = TRUE)
  }
  
  save(se_IBS.gl, se_IBS.mb, se_IBS.slr, se_H.gl, se_H.mb, se_H.slr,
       file=file.path(path.spiec_easi, agg_level, paste0("NetFits_",physeq_name,".RData")))
}
