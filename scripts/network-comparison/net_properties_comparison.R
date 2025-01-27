# **************************************
# Purpose: Network Properties Comparison
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# **************************************

source("tools/functions.R")

get_network_properties <- function(physeq_name, agg_level){
  
  load(file.path(path.assoc_mat, agg_level, 
                 paste0("AssocMat_",physeq_name,".RData")))
  
  props_asso.gl <- summary(network_construct_comparison(assoMat_H.gl, assoMat_IBS.gl))
  props_asso.mb <- summary(network_construct_comparison(assoMat_H.mb, assoMat_IBS.mb))
  props_asso.slr <- summary(network_construct_comparison(assoMat_H.slr, assoMat_IBS.slr))
  
  props_asso.gl$group1 <- "Healthy"
  props_asso.gl$group2 <- "IBS"
  props_asso.mb$group1 <- "Healthy"
  props_asso.mb$group2 <- "IBS"
  props_asso.slr$group1 <- "Healthy"
  props_asso.slr$group2 <- "IBS"
  
  if (!dir.exists(file.path(path.properties, agg_level))) {
    dir.create(file.path(path.properties, agg_level), recursive = TRUE)
  }
  
  save(props_asso.gl, props_asso.mb, props_asso.slr,
       file=file.path(path.properties, agg_level, 
                      paste0("NetProp_",physeq_name,".RData")))
}
