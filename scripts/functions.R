
# ********************************
# Purpose: Help Functions
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ********************************


### Functions for Single Network Analysis


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

network_construct <- function(assoMat, thresh=0){
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoMat,
                             dataType = "condDependence",
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


network_construct_comparison <- function(assoMat1,assoMat2){
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoMat1,
                             data2 = assoMat2,
                             dataType = "condDependence",
                             sparsMethod = "none",
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

meta_analysis <- function(association_matrices, aggregation_method){
  
  # Union of all nodes (row/column names from matrices)
  all_nodes <- unique(unlist(lapply(association_matrices, rownames)))
  
  # Create an empty adjacency matrix for the meta-network
  meta_adj <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
  rownames(meta_adj) <- colnames(meta_adj) <- all_nodes
  
  # Loop over all pairs of nodes to aggregate the edge weights
  for (i in 1:length(all_nodes)) {
    for (j in i:length(all_nodes)) {
      # Collect edge weights for the node pair across all matrices
      edge_weights <- unlist(lapply(association_matrices, function(matrix) {
        if (all(c(all_nodes[i], all_nodes[j]) %in% rownames(matrix))) {
          return(matrix[all_nodes[i], all_nodes[j]])
        }
        return(NA)
      }))
      
      # Remove NA values
      edge_weights <- edge_weights[!is.na(edge_weights)]
      
      # Aggregate edge weights based on the chosen method
      if (length(edge_weights) > 0) {
        if (aggregation_method == "min") {
          meta_adj[i, j] <- min(edge_weights)
          meta_adj[j, i] <- min(edge_weights)
        } else if (aggregation_method == "max") {
          meta_adj[i, j] <- max(edge_weights)
          meta_adj[j, i] <- max(edge_weights)
        } else if (aggregation_method == "mean") {
          meta_adj[i, j] <- mean(edge_weights)
          meta_adj[j, i] <- mean(edge_weights)
        } else if (aggregation_method == "median") {
          meta_adj[i, j] <- median(edge_weights)
          meta_adj[j, i] <- median(edge_weights)
        }
      }
    }
  }
  return(meta_adj)
}

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
  
  #### summary can be: "proportion", "sum" or "mean"
  
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
    }
  }
  
  return(final_matrix)
}
