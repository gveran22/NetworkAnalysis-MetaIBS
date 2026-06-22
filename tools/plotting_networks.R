# ********************************
# Purpose: Help Functions - Plotting
# Date: Juli 2024
# Author: Gilary Evans, Vera Nunez
# ********************************

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


get_venn_diagram <- function(matrix1, matrix2, matrix3, physeq_name, agg_level){
  
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

plot_network <- function(meta,
                         meta2 = NULL,
                         layout = NULL,
                         title = "Meta-Analysis",
                         groupNames = c("Healthy", "IBS"),
                         dataType = "condDependence",
                         thresh = 0,
                         repulsion = 0.7,
                         sameLayout = TRUE,
                         nodeColor = "cluster",
                         labelScale = FALSE,
                         rmSingles = TRUE,
                         nodeSize = "eigenvector",
                         cexNodes = 0.58,
                         cexLabels = 0.5,
                         cexHubLabels = 0.7,
                         showTitle = TRUE,
                         cexTitle = 1.3,
                         hubBorderCol = "gray40") {
  
  props <- network_construct(
    assoMat1 = meta,
    assoMat2 = meta2,
    dataType = dataType,
    thresh = thresh
  )
  
  plot_args <- list(
    x = props,
    layout = layout,
    repulsion = repulsion,
    sameLayout = sameLayout,
    nodeColor = nodeColor,
    labelScale = labelScale,
    rmSingles = rmSingles,
    nodeSize = nodeSize,
    cexNodes = cexNodes,
    cexLabels = cexLabels,
    cexHubLabels = cexHubLabels,
    showTitle = showTitle,
    cexTitle = cexTitle,
    hubBorderCol = hubBorderCol
  )
  
  if (is.null(meta2)) {
    plot_args$title1 <- title
  } else {
    plot_args$groupNames <- groupNames
  }
  
  do.call(plot, plot_args)
  
}

mask_meta_by_matrix <- function(meta, matrix) {
  
  common_rows <- intersect(rownames(meta), rownames(matrix))
  common_cols <- intersect(colnames(meta), colnames(matrix))
  
  filtered_meta <- meta
  filtered_meta[,] <- 0
  
  filtered_meta[common_rows, common_cols] <-
    meta[common_rows, common_cols] *
    (matrix[common_rows, common_cols] != 0)
  
  filtered_meta
}

plot_individual_network <- function(matrix,
                                    meta,
                                    layout,
                                    dataset_name,
                                    method,
                                    cexNodes = 0.58,
                                    cexLabels = 0.7,
                                    cexHubLabels = 1,
                                    cexTitle = 1.5,
                                    repulsion = 0.7,
                                    rmSingles = "none") {
  
  if (sum(matrix != 0, na.rm = TRUE) == 0) {
    message("Skipping network ", dataset_name, " due to no connected nodes.")
    return(NULL)
  }
  
  filtered_meta <- mask_meta_by_matrix(meta, matrix)
  props <- network_construct(filtered_meta)
  
  plot(
    props,
    layout = layout,
    repulsion = repulsion,
    nodeColor = "cluster",
    labelScale = FALSE,
    rmSingles = rmSingles,
    nodeSize = "eigenvector",
    cexNodes = cexNodes,
    cexLabels = cexLabels,
    cexHubLabels = cexHubLabels,
    title1 = paste0(
      "Network Analysis ",
      stringr::str_to_title(dataset_name),
      " ",
      method
    ),
    showTitle = TRUE,
    cexTitle = cexTitle,
    hubBorderCol = "gray40"
  )
  
  invisible(props)
}

plot_individual_networks <- function(matrices,
                                     meta,
                                     layout,
                                     datasets_names,
                                     method,
                                     cexNodes = 0.58,
                                     cexLabels = 0.7,
                                     cexHubLabels = 1,
                                     cexTitle = 1.5,
                                     repulsion = 0.7,
                                     rmSingles = "none") {
  
  for (i in seq_along(matrices)) {
    plot_individual_network(
      matrix = matrices[[i]],
      meta = meta,
      layout = layout,
      dataset_name = datasets_names[[i]],
      method = method,
      cexNodes = cexNodes,
      cexLabels = cexLabels,
      cexHubLabels = cexHubLabels,
      cexTitle = cexTitle,
      repulsion = repulsion,
      rmSingles = rmSingles
    )
  }
}

plot_individual_network_comparison <- function(matrix_H,
                                               matrix_IBS,
                                               meta_H,
                                               meta_IBS,
                                               layout,
                                               dataset_name,
                                               method,
                                               cexNodes = 0.58,
                                               cexLabels = 1.5,
                                               cexHubLabels = 1.7,
                                               cexTitle = 2,
                                               repulsion = 0.7,
                                               rmSingles = TRUE,
                                               sameLayout = FALSE) {
  
  if (sum(matrix_H != 0, na.rm = TRUE) == 0 &&
      sum(matrix_IBS != 0, na.rm = TRUE) == 0) {
    message("Skipping network ", dataset_name, " due to no connected nodes.")
    return(NULL)
  }
  
  filtered_meta_H <- mask_meta_by_matrix(meta_H, matrix_H)
  filtered_meta_IBS <- mask_meta_by_matrix(meta_IBS, matrix_IBS)
  
  props <- network_construct(
    assoMat1 = filtered_meta_H,
    assoMat2 = filtered_meta_IBS
  )
  
  plot(
    props,
    layout = layout,
    repulsion = repulsion,
    sameLayout = sameLayout,
    nodeColor = "cluster",
    labelScale = FALSE,
    rmSingles = rmSingles,
    nodeSize = "eigenvector",
    cexNodes = cexNodes,
    cexLabels = cexLabels,
    cexHubLabels = cexHubLabels,
    groupNames = c(
      paste0("Healthy - ", stringr::str_to_title(dataset_name), " ", method),
      paste0("IBS - ", stringr::str_to_title(dataset_name), " ", method)
    ),
    showTitle = TRUE,
    cexTitle = cexTitle,
    hubBorderCol = "gray40"
  )
  
  invisible(props)
}

plot_individual_networks_comparison <- function(matrices_H,
                                                matrices_IBS,
                                                meta_H,
                                                meta_IBS,
                                                layout,
                                                datasets_names,
                                                method,
                                                cexNodes = 0.58,
                                                cexLabels = 1.5,
                                                cexHubLabels = 1.7,
                                                cexTitle = 2,
                                                repulsion = 0.7,
                                                rmSingles = TRUE,
                                                sameLayout = FALSE) {
  
  for (i in seq_along(matrices_H)) {
    plot_individual_network_comparison(
      matrix_H = matrices_H[[i]],
      matrix_IBS = matrices_IBS[[i]],
      meta_H = meta_H,
      meta_IBS = meta_IBS,
      layout = layout,
      dataset_name = datasets_names[[i]],
      method = method,
      cexNodes = cexNodes,
      cexLabels = cexLabels,
      cexHubLabels = cexHubLabels,
      cexTitle = cexTitle,
      repulsion = repulsion,
      rmSingles = rmSingles,
      sameLayout = sameLayout
    )
  }
}