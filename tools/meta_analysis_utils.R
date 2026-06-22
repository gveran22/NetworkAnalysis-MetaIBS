
# ***************************************
# Purpose: Help Functions - Meta Analysis
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ***************************************


# Function to load matrices and filter NULL entries
load_rdata_object <- function(file, object_name) {
  e <- new.env()
  load(file, envir = e)
  e[[object_name]]
}

load_assoc_matrices <- function(datasets_names, path.assoc_mat, agg_level, object_name) {
  
  matrices <- lapply(datasets_names, function(ds) {
    tryCatch({
      load_rdata_object(
        file.path(path.assoc_mat, agg_level, paste0("AssocMat_", ds, ".RData")),
        object_name
      )
    }, error = function(e) {
      message("Skipping ", ds, ": ", e$message)
      NULL
    })
  })
  
  names(matrices) <- datasets_names
  Filter(Negate(is.null), matrices)
}

############# Functions for meta-network ####################
read_assoc.matrices <- function(agg_level, method){
  # Specify the path to the folder containing CSV files
  order_path <- file.path(path.assoc_mat, agg_level)
  
  # Get a list of all CSV files in the folder
  assoc.mat_files <- list.files(path = order_path, pattern = paste0("*",method,".csv"), 
                                full.names = TRUE)
  
  # Read all CSV files into a list of matrices
  matrices <- lapply(assoc.mat_files, function(file) {
    data <- as.matrix(read.csv(file, header = T, row.names=1))
    colnames(data) <- rownames(data)
    data
  })
  
  return(matrices)
}

summary_assoc_matrix <- function(matrices,
                                 summary_type = c("proportion", "sum", "mean", "mean_present",
                                                  "min", "max", "median", "median_present"),
                                 symmetrize = TRUE,
                                 diag_zero = FALSE) {
  
  summary_type <- match.arg(summary_type)
  
  if (length(matrices) == 0) {
    stop("No matrices provided.")
  }
  
  # Basic checks
  matrices <- lapply(matrices, function(mat) {
    mat <- as.matrix(mat)
    
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("All matrices must have rownames and colnames.")
    }
    
    if (symmetrize) {
      common <- intersect(rownames(mat), colnames(mat))
      mat <- mat[common, common, drop = FALSE]
      mat <- (mat + t(mat)) / 2
    }
    
    if (diag_zero) {
      diag(mat) <- 0
    }
    
    mat
  })
  
  # Union of all taxa
  all_taxa <- sort(unique(unlist(lapply(matrices, function(x) {
    union(rownames(x), colnames(x))
  }))))
  
  # Expand each matrix to the same taxa universe
  full_matrices <- lapply(matrices, function(mat) {
    
    full <- matrix(
      0,
      nrow = length(all_taxa),
      ncol = length(all_taxa),
      dimnames = list(all_taxa, all_taxa)
    )
    
    full[rownames(mat), colnames(mat)] <- mat
    
    if (diag_zero) diag(full) <- 0
    
    full
  })
  
  arr <- simplify2array(full_matrices)
  
  presence_arr <- arr != 0
  presence_count <- apply(presence_arr, c(1, 2), sum)
  
  if (summary_type == "proportion") {
    
    final_matrix <- presence_count / length(matrices)
    
  } else if (summary_type == "sum") {
    
    final_matrix <- apply(arr, c(1, 2), sum)
    
  } else if (summary_type == "mean") {
    
    # Mean across all datasets, absences counted as zero
    final_matrix <- apply(arr, c(1, 2), mean)
    
  } else if (summary_type == "mean_present") {
    
    # Mean only among datasets where edge exists
    final_matrix <- apply(arr, c(1, 2), function(x) {
      x_present <- x[x != 0]
      if (length(x_present) == 0) return(0)
      mean(x_present)
    })
    
  } else if (summary_type == "median") {
    
    # Median across all datasets, absences counted as zero
    final_matrix <- apply(arr, c(1, 2), median)
    
  } else if (summary_type == "median_present") {
    
    # Median only among datasets where edge exists
    final_matrix <- apply(arr, c(1, 2), function(x) {
      x_present <- x[x != 0]
      if (length(x_present) == 0) return(0)
      median(x_present)
    })
    
  } else if (summary_type == "min") {
    
    final_matrix <- apply(arr, c(1, 2), min)
    
  } else if (summary_type == "max") {
    
    final_matrix <- apply(arr, c(1, 2), max)
  }
  
  if (symmetrize) {
    final_matrix <- (final_matrix + t(final_matrix)) / 2
  }
  
  if (diag_zero) {
    diag(final_matrix) <- 0
  }
  
  return(final_matrix)
}

build_meta_matrix <- function(matrices, summary_type = "mean") {
  summary_assoc_matrix(matrices, summary_type)
}

get_quantile_threshold <- function(mat, prob = 0.90) {
  
  vals <- abs(mat[upper.tri(mat)])
  vals <- vals[vals > 0]
  
  if (length(vals) == 0) {
    return(0)
  }
  
  unname(stats::quantile(vals, probs = prob, na.rm = TRUE))
}

process_global_properties <- function(datasets_names,
                                      path.properties,
                                      agg_level,
                                      object_name,
                                      comparison = FALSE,
                                      group_names = c("Healthy", "IBS")) {
  
  properties <- lapply(datasets_names, function(ds) {
    tryCatch({
      load_rdata_object(
        file = file.path(path.properties, agg_level, paste0("NetProp_", ds, ".RData")),
        object_name = object_name
      )
    }, error = function(e) {
      message("Skipping ", ds, ": ", e$message)
      NULL
    })
  })
  
  names(properties) <- datasets_names
  properties <- Filter(Negate(is.null), properties)
  
  glob_probs <- lapply(properties, function(x) {
    
    out <- x$glob_probs_lcc
    
    if (comparison) {
      colnames(out) <- group_names
    }
    
    out
  })
  
  merged_df <- do.call(cbind, glob_probs)
  
  if (!comparison) {
    colnames(merged_df) <- names(properties)
  }
  
  merged_df
}

get_common_associations <- function(prop.gl,
                                    prop.mb,
                                    prop.slr,
                                    thresh = 0.4,
                                    min_methods = 2) {
  
  gl <- prop.gl >= thresh
  mb <- prop.mb >= thresh
  slr <- prop.slr >= thresh
  
  votes <- gl + mb + slr
  common <- votes >= min_methods
  
  common <- matrix(
    as.numeric(common),
    nrow = nrow(common),
    dimnames = dimnames(common)
  )
  
  diag(common) <- 0
  
  common
}

build_meta_matrix <- function(matrices,
                              agg_level,
                              thresh = 0.1,
                              summary_type = "proportion") {
  
  meta_layout_matrix <- summary_assoc_matrix(
    matrices,
    summary_type = summary_type
  )
  
  layout_matrix <- meta_layout_matrix
  layout_matrix[layout_matrix < thresh] <- 0
  
  layout_matrix
}
