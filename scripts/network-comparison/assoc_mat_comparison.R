# *************************************
# Purpose: Getting Association Matrices
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# *************************************

source("tools/functions.R")

get_assoc_matrix <- function(physeq_name, agg_level){
  
  load(file.path(path.spiec_easi, agg_level, 
                 paste0("NetFits_",physeq_name,".RData")))
  
  # Association Matrices glasso
  secor_IBS <- stats::cov2cor(as.matrix(getOptCov(se_IBS.gl)))
  assoMat_IBS.gl <- as.matrix(secor_IBS * SpiecEasi::getRefit(se_IBS.gl))
  
  secor_H <- stats::cov2cor(as.matrix(getOptCov(se_H.gl)))
  assoMat_H.gl <- as.matrix(secor_H * SpiecEasi::getRefit(se_H.gl))
  
  # Association matrix mb
  assoMat_IBS.mb <- as.matrix(SpiecEasi::symBeta(SpiecEasi::getOptBeta(se_IBS.mb)))
  assoMat_H.mb <- as.matrix(SpiecEasi::symBeta(SpiecEasi::getOptBeta(se_H.mb)))
  
  # Association matrix slr
  opt_IBS <- getOptSLR(se_IBS.slr) # getting optimal se.slr
  icov <- Matrix::drop0(opt_IBS$est$icov[[getOptInd(opt_IBS)]])
  secor <- cov2cor(prec2cov(icov))
  assoMat_IBS.slr <- as.matrix(secor*getRefit(opt_IBS))
  
  opt_H <- getOptSLR(se_H.slr) # getting optimal se.slr
  icov <- Matrix::drop0(opt_H$est$icov[[getOptInd(opt_H)]])
  secor <- cov2cor(prec2cov(icov))
  assoMat_H.slr <- as.matrix(secor*getRefit(opt_H))
  
  rownames(assoMat_IBS.gl) <- colnames(assoMat_IBS.gl) <- colnames(se_IBS.gl$est$data)
  rownames(assoMat_IBS.mb) <- colnames(assoMat_IBS.mb) <- colnames(se_IBS.mb$est$data)
  rownames(assoMat_IBS.slr) <- colnames(assoMat_IBS.slr) <- colnames(opt_IBS$est$data)
  
  rownames(assoMat_H.gl) <- colnames(assoMat_H.gl) <- colnames(se_H.gl$est$data)
  rownames(assoMat_H.mb) <- colnames(assoMat_H.mb) <- colnames(se_H.mb$est$data)
  rownames(assoMat_H.slr) <- colnames(assoMat_H.slr) <- colnames(opt_H$est$data)
  
  diag(assoMat_IBS.gl) <- diag(assoMat_IBS.mb) <- diag(assoMat_IBS.slr) <- 1
  diag(assoMat_H.gl) <- diag(assoMat_H.mb) <- diag(assoMat_H.slr) <- 1
  
  if (!dir.exists(file.path(path.assoc_mat, agg_level))) {
    dir.create(file.path(path.assoc_mat, agg_level), recursive = TRUE)
  }
  save(assoMat_IBS.gl, assoMat_IBS.mb, assoMat_IBS.slr, 
       assoMat_H.gl, assoMat_H.mb, assoMat_H.slr,
       file=file.path(path.assoc_mat, agg_level, 
                      paste0("AssocMat_",physeq_name,".RData")))
}
