
# ********************************
# Purpose: Help Functions -NetCoMi
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ********************************


network_construct <- function(assoMat1,
                              assoMat2 = NULL,
                              dataType = "condDependence",
                              thresh = 0) {
  
  if (is.null(assoMat2)) {
    
    net_asso <- netConstruct(
      data = assoMat1,
      dataType = dataType,
      sparsMethod = "threshold",
      thresh = thresh,
      normMethod = "none",
      zeroMethod = "none",
      dissFunc = "signed",
      verbose = 0
    )
    
  } else {
    
    net_asso <- netConstruct(
      data = assoMat1,
      data2 = assoMat2,
      dataType = dataType,
      sparsMethod = "threshold",
      thresh = thresh,
      normMethod = "none",
      zeroMethod = "none",
      dissFunc = "signed",
      verbose = 0
    )
  }
  
  netAnalyze(
    net_asso,
    centrLCC = FALSE,
    avDissIgnoreInf = TRUE,
    sPathNorm = FALSE,
    clustMethod = "cluster_fast_greedy",
    hubPar = "eigenvector",
    hubQuant = 0.9,
    gcmHeat = FALSE
  )
}

extract_netcomi_summary <- function(props, method, dataset, agg_level) {
  
  gp <- props$glob_probs
  gp_lcc <- props$glob_probs_lcc
  
  data.frame(
    dataset = dataset,
    agg_level = agg_level,
    method = method,
    
    n_components = as.numeric(gp["Number of components", 1]),
    edge_density_global = as.numeric(gp["Edge density", 1]),
    clustering_global = as.numeric(gp["Clustering coefficient", 1]),
    modularity_global = as.numeric(gp["Modularity", 1]),
    natural_connectivity_global = as.numeric(gp["Natural connectivity", 1]),
    
    rel_lcc_size = as.numeric(gp_lcc["Relative LCC size", 1]),
    edge_density_lcc = as.numeric(gp_lcc["Edge density", 1]),
    avg_path_lcc = as.numeric(gp_lcc["Average path length**", 1])
  )
}