
# ********************************
# Purpose: Help Functions
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ********************************


######## Functions for Single and Comparison Network Analysis  ##################

# Getting best slr output
getOptSLR <- function(x) {
  ind <- x$ebic>1e-1
  x[ind][[which.min(x$ebic[ind])]]
}

# Defining ranks path for slr spiec-easi method
rank.slr <- function(otu){
  dim <- ncol(otu) - 4
  
  #Definig ranks
  if(dim>64){
    ranks <- round(seq(2, 64, len=10))
  }else{
    ranks <- round(seq(2, dim, len=8))
  }
  
  return(ranks)
}

network_construct <- function(assoMat, dataType = "condDependence", thresh=0){
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoMat,
                             dataType = dataType,
                             sparsMethod = "threshold",
                             thresh = thresh,
                             normMethod = "none", 
                             zeroMethod = "none",
                             dissFunc = "signed",
                             verbose = 0)
  
  props_asso <- netAnalyze(net_asso_p, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           hubQuant = 0.9,
                           gcmHeat = FALSE)
  return(props_asso)
}


network_construct_comparison <- function(assoMat1,assoMat2, dataType = "condDependence", thresh = 0){
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoMat1,
                             data2 = assoMat2,
                             dataType = dataType,
                             sparsMethod = "threshold",
                             thresh = thresh,
                             normMethod = "none", 
                             zeroMethod = "none",
                             dissFunc = "signed",
                             verbose = 0)

  props_asso <- netAnalyze(net_asso_p, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           hubQuant = 0.9,
                           gcmHeat = FALSE)
  return(props_asso)
}


get_network <- function(assoc.matrix, network_name, method, 
                        dataType = "condDependence", sparsMethod = "none",
                        thresh = 0, doPlot=T){
  
  #### dataType can also be: "proportionality", "counts", "correlation", "partialCorr" 
  #### sparsMethod can also be: "threshold", "softThreshold"
  #### thresh: if sparsMethod is changed to "threshold", then you need to define a value.
  #### doPlot:If FALSE, the network plot is suppressed. Useful for saving the output (e.g., the layout) without plotting.
  
  # Network construction and analysis
  
  net_asso <- netConstruct(assoc.matrix,
                           dataType = dataType,
                           sparsMethod = sparsMethod, 
                           thresh = thresh, 
                           normMethod = "none", 
                           zeroMethod = "none",
                           dissFunc = "signed",
                           verbose = 3)
  
  
  props_asso <- netAnalyze(net_asso, 
                           #centrLCC = FALSE,
                           #avDissIgnoreInf = TRUE,
                           #sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           #hubPar = "eigenvector", #c("degree", 
                           #hubQuant = 0.9,
                           #lnormFit = TRUE,
                           #normDeg = FALSE,
                           #normBetw = FALSE,
                           #normClose = FALSE,
                           #normEigen = FALSE,
                           gcmHeat = FALSE)
  
  plot(props_asso, 
       sameLayout = TRUE, 
       repulsion = 0.7,
       nodeColor = "cluster",
       labelScale = FALSE,
       rmSingles = TRUE,
       nodeSize = "eigenvector",
       cexNodes = 1, 
       cexLabels = 1,
       cexHubLabels = 1,
       title1 = paste0(str_to_title(network_name), " Network - Method: ", method), 
       showTitle = TRUE,
       cexTitle = 2,
       hubBorderCol  = "gray40")
  
}


get_venn_diagram <- function(matrix1, matrix2, matrix3, physeq_name){
  
  get_named_associations <- function(matrix) {
    # Find non-zero elements in the upper triangle (excluding diagonal)
    non_zero_indices <- which(matrix != 0 & upper.tri(matrix), arr.ind = TRUE)
    
    # Convert these indices to unique named identifiers
    unique_pairs <- apply(non_zero_indices, 1, function(x) paste(rownames(matrix)[x[1]], colnames(matrix)[x[2]], sep = "-"))
    
    return(unique(unique_pairs))
  }
  
  # Get unique pairs of non-zero elements for each matrix
  set1_ids <- get_named_associations(matrix1)
  set2_ids <- get_named_associations(matrix2)
  set3_ids <- get_named_associations(matrix3)
  
  #png(file.path(path.plots, paste0("Venn_diagram_",physeq_name,".png")), width = 800, height = 600) 
  # Create the Venn diagram
  venn.plot <- venn.diagram(
    x = list(Matrix1 = set1_ids, Matrix2 = set2_ids, Matrix3 = set3_ids),
    category.names = c("Glasso", "Mb", "Slr"),
    filename = NULL,
    output = TRUE,
    
    # Custom colors for circles
    lwd = 2,                        # Line width
    col = c("red", "green", "blue"),# Colors of the circles
    fill = c("red", "green", "blue"),# Fill colors of the circles
    alpha = 0.3,                    # Transparency of fills
    
    # Custom text
    cex = 2,                        # Size of category labels
    fontface = "bold",              # Bold category labels
    fontfamily = "serif",           # Font family of category labels
    
    # Custom labels inside the diagram
    cat.col = c("darkred", "darkgreen", "darkblue"), # Colors for the category labels
    cat.cex = 2,                     # Size of category labels
    cat.fontface = "bold",           # Bold category labels
    cat.fontfamily = "serif",         # Font family of category labels
    
    main = paste0("Venn Diagram ",physeq_name),
    main.cex = 2,  # Adjust the size of the title text
    main.fontface = "bold",  # Make the title bold
    main.fontfamily = "serif"  # Specify the font family
  )
  #dev.off()
  grid.newpage()
  png(file.path(path.venn_diag, agg_level, paste0("Venn_diagram_",physeq_name,".png")), width = 800, height = 600) 
  grid.draw(venn.plot)
  dev.off()
  
}

########################### Functions for meta-analysis ########################################

# Function to load matrices and filter NULL entries
load_matrices <- function(datasets_names, path.assoc_mat, agg_level, object_index) {
  matrices <- lapply(datasets_names, function(x) {
    tryCatch({
      load(file = file.path(path.assoc_mat, agg_level, paste0("AssocMat_", x, ".RData")))
      objects <- ls()  # Get loaded objects
      return(get(objects[object_index]))  # Return the specified object
    }, error = function(e) {
      message(paste("Skipping file:", x, "due to error:", e$message))
      return(NULL)
    })
  })
  Filter(Negate(is.null), matrices)  # Remove NULL entries
}

# Function to plot meta-analysis
plot_network <- function(meta, layout = NULL, title = "Meta-Analysis", dataType = "condDependence", 
                         thresh = 0, repulsion=0.7) {
  props_asso_meta <- network_construct(meta, dataType = dataType , thresh = thresh)
  plot(props_asso_meta,
       layout = layout,
       repulsion = repulsion,
       sameLayout = TRUE,
       nodeColor = "cluster",
       labelScale = FALSE,
       rmSingles = TRUE,
       nodeSize = "eigenvector",
       cexNodes = 0.58,
       cexLabels = 0.7,
       cexHubLabels = 1,
       title1 = title,
       showTitle = TRUE,
       cexTitle = 1.5,
       hubBorderCol = "gray40")
}

# Function to plot meta-analysis comparison
plot_network_comparison <- function(meta_H, meta_IBS, layout = NULL, groupNames = c("Healthy", "IBS"), 
                                    dataType = "condDependence", thresh = 0, repulsion =0.7) {
  props_asso_meta <- network_construct_comparison(meta_H, meta_IBS, dataType = dataType , thresh = thresh)
  plot(props_asso_meta,
       layout = layout,
       repulsion = repulsion,
       sameLayout = TRUE,
       nodeColor = "cluster",
       labelScale = FALSE,
       rmSingles = TRUE,
       nodeSize = "eigenvector",
       cexNodes = 0.58,
       cexLabels = 0.7,
       cexHubLabels = 1,
       #title1 = title,
       groupNames = groupNames,
       showTitle = TRUE,
       cexTitle = 1.5,
       hubBorderCol = "gray40")
}


# Function to process individual plots
plot_individual_network <- function(matrix, meta, layout,  datasets_names, method) {
  
  sub_adj <- matrix
  
  # Skip if no connected nodes
  if (sum(sub_adj)== 0) {
    message(paste0("Skipping network ", datasets_names[i], " due to no connected nodes."))
    next
  }
  
  common_rows <- intersect(rownames(meta), rownames(sub_adj))
  common_cols <- intersect(colnames(meta), colnames(sub_adj))
  
  filtered_meta_adj <- meta
  filtered_meta_adj[,] <- 0
  
  filtered_meta_adj[common_rows, common_cols] <- meta[common_rows, common_cols] * (sub_adj[common_rows, common_cols] != 0)
  
  props_asso <- network_construct(filtered_meta_adj)
  plot(props_asso,
       layout = layout,
       nodeColor = "cluster",
       labelScale = FALSE,
       rmSingles = "none",#TRUE,
       nodeSize = "eigenvector",
       cexNodes = 0.58,
       cexLabels = 0.7,
       cexHubLabels = 1,
       title1 = paste0("Network Analysis ",  str_to_title(datasets_names)," ",method),
       showTitle = TRUE,
       cexTitle = 1.5,
       hubBorderCol  = "gray40")
  
}

# Function to process individual plots as loop
plot_individual_networks <- function(matrices, meta, layout, datasets_names, method) {
  for (i in seq_along(matrices)) {
    plot_individual_network(matrices[[i]], meta, layout,  datasets_names[[i]], method)
  }
}

# Function to process individual plots comparison
plot_individual_network_comparison <- function(matrix_H, matrix_IBS, meta_H, meta_IBS, 
                                                layout, datasets_names, method) {
  sub_adj_H <- matrix_H
  sub_adj_IBS <- matrix_IBS
  
  # Skip if no connected nodes
  if (sum(sub_adj_H) == 0 && sum(sub_adj_IBS) == 0) {
    message(paste0("Skipping network ", datasets_names, " due to no connected nodes."))
    next
  }
  
  common_rows_H <- intersect(rownames(meta_H), rownames(sub_adj_H))
  common_cols_H <- intersect(colnames(meta_H), colnames(sub_adj_H))
  common_rows_IBS <- intersect(rownames(meta_IBS), rownames(sub_adj_IBS))
  common_cols_IBS <- intersect(colnames(meta_IBS), colnames(sub_adj_IBS))
  
  filtered_meta_adj_H <- meta_H
  filtered_meta_adj_IBS <- meta_IBS
  filtered_meta_adj_H[,] <- 0
  filtered_meta_adj_IBS[,] <- 0
  
  filtered_meta_adj_H[common_rows_H, common_cols_H] <- meta_H[common_rows_H, common_cols_H] * (sub_adj_H[common_rows_H, common_cols_H] != 0)
  filtered_meta_adj_IBS[common_rows_IBS, common_cols_IBS] <- meta_IBS[common_rows_IBS, common_cols_IBS] * (sub_adj_IBS[common_rows_IBS, common_cols_IBS] != 0)
  
  props_asso <- network_construct_comparison(filtered_meta_adj_H, filtered_meta_adj_IBS)
  plot(props_asso,
       layout = layout,
       #sameLayout = TRUE,
       nodeColor = "cluster",
       labelScale = FALSE,
       rmSingles = TRUE,
       nodeSize = "eigenvector",
       cexNodes = 0.58,
       cexLabels = 1.5,
       cexHubLabels = 1.7,
       groupNames = c(paste0("Healthy - ", str_to_title(datasets_names)," ",method), 
                      paste0("IBS - ", str_to_title(datasets_names)," ",method)),
       cexTitle = 2,
       hubBorderCol = "gray40")
  
}

# Function to process individual plots comparison as loop
plot_individual_networks_comparison <- function(matrices_H, matrices_IBS, meta_H, meta_IBS, 
                                                layout, datasets_names, method) {
  for (i in seq_along(matrices_H)) {
    plot_individual_network_comparison(matrices_H[[i]], matrices_IBS[[i]], meta_H, meta_IBS, 
                                        layout, datasets_names[i], method)
  }
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


summary_assoc_matrix <- function(matrices, summary){
  
  #### summary can be: "proportion", "sum", "mean", "min", "max", o "median"
  
  # Find the union of all row and column names
  all_rows <- unique(unlist(lapply(matrices, rownames)))
  all_cols <- unique(unlist(lapply(matrices, colnames)))
  
  # Initialize an empty matrix with the unified structure
  full_matrix <- matrix(0, nrow = length(all_rows), ncol = length(all_cols))
  rownames(full_matrix) <- all_rows
  colnames(full_matrix) <- all_cols
  
  if (summary=="proportion") {
    # Initialize the count matrix to count the presence of non-zero elements
    count_matrix <- full_matrix
    
    # Fill the full matrices and update the count matrix
    for (mat in matrices) {
      temp <- full_matrix
      temp[rownames(mat), colnames(mat)] <- mat
      count_matrix <- count_matrix + (temp != 0)
    }
    
    # Calculate the proportion of presence
    final_matrix <- count_matrix / length(matrices)
    
  }else{
    # Initialize a list of full matrices with zeros
    full_matrices <- lapply(matrices, function(x) {
      temp <- full_matrix
      temp[rownames(x), colnames(x)] <- x
      return(temp)
    })
    
    # Sum the matrices element-wise
    sum_matrix <- Reduce("+", full_matrices)
    
    if (summary=="sum") {
      final_matrix <- sum_matrix
      
    }else if (summary=="mean"){
      # Calculate the mean matrix
      final_matrix <- sum_matrix / length(matrices)
    }else if (summary %in% c("min", "max", "median")) {
      # Convierte la lista de matrices en un array tridimensional
      matrix_array <- simplify2array(full_matrices)
      
      # Calcula min, max o median usando apply
      if (summary == "min") {
        final_matrix <- apply(matrix_array, c(1, 2), min)
      } else if (summary == "max") {
        final_matrix <- apply(matrix_array, c(1, 2), max)
      } else if (summary == "median") {
        final_matrix <- apply(matrix_array, c(1, 2), median)
      }
    } else {
      stop("El método de resumen no es válido. Usa 'proportion', 'sum', 
           'mean', 'min', 'max', o 'median'.")
    }
  }
  
  return(final_matrix)
}

process_global_properties <- function(datasets_names, path.properties, 
                                      agg_level, object_index) {
  # Load properties
  properties <- lapply(datasets_names, function(x) {
    tryCatch({
      load(file = file.path(path.properties, agg_level, paste0("NetProp_", x, ".RData")))
      objects <- ls()
      return(get(objects[object_index]))
    }, error = function(e) {
      # Handle the empty network error
      NULL
    })
    
  })
  
  # Assign dataset names to the properties
  names(properties) <- datasets_names
  properties <- Filter(Negate(is.null), properties)
  
  # Extract and process global probabilities
  glob_probs <- lapply(properties, function(data) {
    data <- data$glob_probs_lcc
    data
  })
  
  # Combine into a single data frame
  merged_df <- do.call(cbind, glob_probs)
  colnames(merged_df) <- datasets_names
  return(merged_df)
}

process_global_properties_comparison <- function(datasets_names, path.properties, 
                                      agg_level, object_index) {
  # Load properties
  properties <- lapply(datasets_names, function(x) {
    tryCatch({
      load(file = file.path(path.properties, agg_level, paste0("NetProp_", x, ".RData")))
      objects <- ls()
      return(get(objects[object_index]))
    }, error = function(e) {
      # Handle the empty network error
      NULL
    })
    
  })
  
  # Assign dataset names to the properties
  names(properties) <- datasets_names
  properties <- Filter(Negate(is.null), properties)
  
  # Extract and process global probabilities
  glob_probs <- lapply(properties, function(data) {
    data <- data$glob_probs_lcc
    colnames(data) <- c("Healthy", "IBS")
    data
  })
  
  # Combine into a single data frame
  merged_df <- do.call(cbind, glob_probs)
  return(merged_df)
}

