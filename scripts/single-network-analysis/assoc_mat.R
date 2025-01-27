# *************************************
# Purpose: Getting Association Matrices
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# *************************************

source("tools/functions.R")

get_assoc_matrix <- function(physeq_name, agg_level){
  
  load(file.path(path.spiec_easi, agg_level, 
                 paste0("NetFits_",physeq_name,".RData")))
  
  # Association Matrix glasso
  secor <- stats::cov2cor(as.matrix(getOptCov(se.gl)))
  assoMat.gl <- as.matrix(secor * SpiecEasi::getRefit(se.gl))
  
  # Association matrix mb
  assoMat.mb <- as.matrix(SpiecEasi::symBeta(SpiecEasi::getOptBeta(se.mb)))
  
  # Association matrix slr
  opt <- getOptSLR(se.slr) # getting optimal se.slr
  icov <- Matrix::drop0(opt$est$icov[[getOptInd(opt)]])
  secor <- cov2cor(prec2cov(icov))
  assoMat.slr <- as.matrix(secor*getRefit(opt))
  
  rownames(assoMat.gl) <- colnames(assoMat.gl) <- colnames(se.gl$est$data)
  rownames(assoMat.mb) <- colnames(assoMat.mb) <- colnames(se.mb$est$data)
  rownames(assoMat.slr) <- colnames(assoMat.slr) <- colnames(opt$est$data)
  
  diag(assoMat.gl) <- 1
  diag(assoMat.mb) <- 1
  diag(assoMat.slr) <- 1
  
  if (!dir.exists(file.path(path.assoc_mat, agg_level))) {
    dir.create(file.path(path.assoc_mat, agg_level), recursive = TRUE)
  }
  save(assoMat.gl, assoMat.mb, assoMat.slr, file=file.path(path.assoc_mat, agg_level, 
                                            paste0("AssocMat_",physeq_name,".RData")))
}
