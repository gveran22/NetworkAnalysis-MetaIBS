
# ***********************************
# Purpose: Help Functions - Spieceasi
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ***********************************

# Getting best slr output
getOptSLR <- function(x) {
  
  ebic <- x$ebic
  
  valid <- is.finite(ebic) & ebic > 1e-1
  
  if (!any(valid)) {
    warning("No finite SLR models with ebic > 1e-1. Using minimum finite EBIC over all ranks.")
    
    finite <- is.finite(ebic)
    
    if (!any(finite)) {
      stop("No finite EBIC values available for SLR model selection.")
    }
    
    opt_index <- which(finite)[which.min(ebic[finite])]
    
  } else {
    opt_index <- which(valid)[which.min(ebic[valid])]
  }
  
  list(
    fit = x[[opt_index]],
    opt_index = opt_index,
    opt_ebic = ebic[opt_index],
    valid_ebic = valid
  )
}

#getOptSLR <- function(x) {
#  ind <- x$ebic>1e-1
#  x[ind][[which.min(x$ebic[ind])]]
#}

# Defining ranks path for slr spiec-easi method
rank.slr <- function(otu) {
  
  p <- ncol(otu)
  max_r <- p - 1
  
  if (max_r < 2) {
    stop("Too few taxa for SLR.")
  }
  
  if (max_r <= 20) {
    
    ranks <- 1:max_r
    
  } else if (max_r <= 80) {
    
    ranks <- c(
      1:20,
      seq(25, max_r, by = 5)
    )
    
  } else {
    
    ranks <- c(
      1:20,
      seq(25, 80, by = 5),
      seq(90, max_r, by = 10)
    )
  }
  
  unique(ranks[ranks < p])
}
