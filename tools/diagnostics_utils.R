
# *****************************************
# Purpose: Help Functions - SLR Diagnostics
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# *****************************************

get_slr_rank_diagnostics <- function(physeq_name, agg_level) {
  
  load(file.path(path.spiec_easi, agg_level,
                 paste0("NetFits_", physeq_name, ".RData")))
  
  out <- lapply(seq_along(se.slr), function(i) {
    
    fit <- se.slr[[i]]
    adj <- SpiecEasi::getRefit(fit)
    
    p <- ncol(adj)
    edges <- sum(adj) / 2
    density <- edges / (p * (p - 1) / 2)
    
    data.frame(
      dataset = physeq_name,
      agg_level = agg_level,
      rank = ranks[i],
      ebic = se.slr$ebic[i],
      nodes = p,
      edges = edges,
      density = density
    )
  })
  
  dplyr::bind_rows(out)
}

collect_slr_rank_diagnostics <- function(datasets_names, agg_levels) {
  
  out <- list()
  k <- 1
  
  for (tax in agg_levels) {
    for (ds in datasets_names) {
      out[[k]] <- get_slr_rank_diagnostics(ds, tax)
      k <- k + 1
    }
  }
  
  dplyr::bind_rows(out)
}