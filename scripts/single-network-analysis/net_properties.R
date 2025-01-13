# ********************************
# Purpose: Network Properties
# Date: Mai 2024
# Author: Gilary Evans, Vera Nunez
# ********************************

source("scripts/functions.R")

get_network_properties <- function(physeq_name, agg_level){
  
  load(file.path(path.assoc_mat, agg_level, 
                 paste0("AssocMat_",physeq_name,".RData")))
  
  props_asso.gl <- summary(network_construct(assoMat.gl))
  props_asso.mb <- summary(network_construct(assoMat.mb))
  props_asso.slr <- summary(network_construct(assoMat.slr))
  
  if (!dir.exists(file.path(path.properties, agg_level))) {
    dir.create(file.path(path.properties, agg_level), recursive = TRUE)
  }
  save(props_asso.gl, props_asso.mb, props_asso.slr, file=file.path(path.properties, agg_level, 
                                                           paste0("NetProp_",physeq_name,".RData")))
}

