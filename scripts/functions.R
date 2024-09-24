
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

network_construct <- function(assoMat){
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoMat,
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


network_analysis_netcomi <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar, 
                                     method, sel.criterion ="stars"){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  pargs <- list(thresh=0.05, rep.num=60, seed=10010, ncores=20)
  
  
  for (m in 1:length(method)) {

    spieceasi_par <- list(nlambda=100,
                             lambda.min.ratio=1e-3,
                             sel.criterion = sel.criterion,
                             method = method[m], 
                             pulsar.select=TRUE,
                             pulsar.params=pargs)
    
    net_asso <- netConstruct(physeq, 
                             taxRank = agg_level,
                             filtTax =  filtTax,
                             filtTaxPar = filtTaxPar,
                             filtSamp =   filtSamp,
                             filtSampPar =  filtSampPar,
                             measure = "spieceasi",
                             measurePar = spieceasi_par,
                             normMethod = "none", 
                             zeroMethod = "none",
                             sparsMethod = "none", 
                             dissFunc = "signed",
                             verbose = 3)
    
    
    #Saving association matrix
    write.csv(net_asso$assoMat1, file = file.path(path.assoc_mat, agg_level,
                                                paste0("AssocMat_",str_to_title(physeq_name),"_",
                                                       method[m],".csv")), row.names = T)
    
    # props_asso <- netAnalyze(net_asso, clustMethod = "hierarchical", gcmHeat = FALSE)
    
    props_asso <- netAnalyze(net_asso, 
                             centrLCC = FALSE,
                             avDissIgnoreInf = TRUE,
                             sPathNorm = FALSE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = c("degree", "eigenvector"),
                             hubQuant = 0.9,
                             lnormFit = TRUE,
                             normDeg = FALSE,
                             normBetw = FALSE,
                             normClose = FALSE,
                             normEigen = FALSE,
                             gcmHeat = FALSE)
    
    # Saving the network properties
    saveRDS(summary(props_asso), file.path(path.properties,  agg_level, 
                                           paste0("Summary_",str_to_title(physeq_name),"_",
                                                  method[m],".rds")))
    
    png(file.path(path.plots, agg_level, paste0("Network_",physeq_name,"_",method[m],".png")),
        width=900, height=700)
    
    plot(props_asso,  
         rmSingles = TRUE,
         repulsion = 0.7,
         nodeColor = "cluster",
         nodeSize = "eigenvector",
         labelScale = FALSE,
         cexNodes = 1, 
         cexLabels = 1,
         cexHubLabels = 1,
         cexTitle = 2,
         title1 = paste0(str_to_title(physeq_name), " Network - Method: ", method[m]), 
         showTitle = TRUE,
         hubBorderCol  = "gray40")

    dev.off()
    
  }
  
}

network_analysis_slr <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar, 
                                 method, sel.criterion ="stars"){
  
  # Read the phyloseq object
  physeq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  # Filtering
  net_asso <- netConstruct(physeq,
                           taxRank = agg_level,
                           filtTax =  filtTax,
                           filtTaxPar = filtTaxPar,
                           filtSamp =   filtSamp,
                           filtSampPar =  filtSampPar,
                           measure = NULL,
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3)
  
  dim <- ncol(net_asso$countMat1) - 4
  
  #Definig ranks
  if(dim>64){
    ranks <- round(seq(2, 64, len=10))
  }else{
    ranks <- round(seq(2, dim, len=8))
  }
  
  se.slr <- spiec.easi(net_asso$countMat1 , method='slr', nlambda=50,
                       lambda.min.ratio=1e-2, r=ranks, lambda.log=TRUE,
                       pulsar.params=list(thresh=0.1, ncores=20, rep.num=30))
  
  se.slr$ebic <- sapply(se.slr, function(x)
    ebic(x$refit$stars, x$est$data,
         x$est$loglik[x$select$stars$opt.index]))
  
  # Getting optimal se
  opt <- getOptSLR(se.slr)
  
  r <- opt$select$envir$ocall$...$r
  
  # Write to a text file
  write.csv(r, file = file.path(path.ranks,agg_level, paste0("rank_",physeq_name,"_slr.csv")), row.names = F)
  
  icov <- Matrix::drop0(opt$est$icov[[getOptInd(opt)]])
  secor <- cov2cor(prec2cov(icov))
  assoc.matrix <- as.matrix(secor*getRefit(opt))
  
  rownames(assoc.matrix) <- colnames(assoc.matrix) <- colnames(net_asso$countMat1)
  
  #Saving association matrix
  write.csv(assoc.matrix, file = file.path(path.assoc_mat, agg_level,
                                                paste0("AssocMat_",str_to_title(physeq_name),"_slr.csv")), row.names = T)
  
  # Network construction and analysis
  net_asso_p <- netConstruct(data = assoc.matrix,
                             dataType = "condDependence",
                             sparsMethod = "none",
                             verbose = 0)
  
  props_asso <- netAnalyze(net_asso_p, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE,
                           gcmHeat = FALSE)
  
  # Saving the network properties
  saveRDS(summary(props_asso), file.path(path.properties,  agg_level, 
                                         paste0("Summary_",str_to_title(physeq_name),"_slr.rds")))
  
  png(file.path(path.plots, agg_level, paste0("Network_",physeq_name,"_slr.png")),
      width=900, height=700)
  
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
       title1 = paste0(str_to_title(physeq_name), " Network - Method: slr"),
       showTitle = TRUE,
       cexTitle = 2,
       hubBorderCol  = "gray40")
  
  dev.off()
}



## Network comparison Analysis

network_comparison_netcomi <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar,
                               method, self.criterion ="stars",  assoMat = T){
  
  # Read the phyloseq object
  phyloseq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  # Renaming ASVs for agg_level name taxonomy
  otutab <- as.data.frame(otu_table(phyloseq))
  
  if(identical(colnames(otutab), taxa_names(phyloseq))){
    colnames(otutab) <- as.data.frame(tax_table(phyloseq))[[agg_level]]
  }else{
    print("ASVs in OTU Table does not match the ASVs in the Taxonomy Table")
  }
  
  otutab <- otu_table(as.matrix(otutab), taxa_are_rows=FALSE)
  phyloseq@otu_table <- otutab
  taxa <- as.data.frame(phyloseq@tax_table)
  
  rownames(phyloseq@tax_table) <- taxa[,agg_level]
  
  variable_names <- c("IBS", "Healthy")
  
  phyloseq_IBS <- phyloseq::subset_samples(phyloseq, host_disease == "IBS")
  phyloseq_IBS <- prune_taxa(taxa_sums(phyloseq_IBS)>0, phyloseq_IBS)
  
  phyloseq_H <- phyloseq::subset_samples(phyloseq, host_disease == "Healthy")
  phyloseq_H <- prune_taxa(taxa_sums(phyloseq_H)>0, phyloseq_H)
  
  # Defining pulsar parameters
  pargs <- list(thresh=0.1, rep.num=60, seed=10010, ncores=20)
  
  for (m in 1:length(method)) {
    
    net_IBS <- netConstruct(data= phyloseq_IBS,
                            data2= phyloseq_H,
                            #group = group_vec, 
                            filtTax = filtTax, # "numbSamp"
                            filtTaxPar = filtTaxPar, # list(numbSamp = 5)
                            filtSamp = filtSamp, # "totalReads",
                            filtSampPar = filtSampPar, # list(totalReads = 500),
                            measure = "spieceasi",
                            measurePar = list(nlambda=100,
                                              lambda.min.ratio=1e-3,
                                              sel.criterion = sel.criterion,
                                              method = method[m], 
                                              pulsar.select=TRUE,
                                              pulsar.params=pargs),
                            normMethod = "none", 
                            zeroMethod = "none",
                            sparsMethod = "none", 
                            dissFunc = "signed",
                            verbose = 3)
    
    
    #Saving association matrix
    write.csv(net_IBS$assoMat1, file = file.path(path.assoc_mat, agg_level,
                                                  paste0("AssocMat_",str_to_title(physeq_name), "_" ,variable_names[1],"_",
                                                        method[m],".csv")), row.names = T)
    write.csv(net_IBS$assoMat2, file = file.path(path.assoc_mat, agg_level,
                                                 paste0("AssocMat_",str_to_title(physeq_name), "_" ,variable_names[2],"_",
                                                        method[m],".csv")), row.names = T)

    
    props_IBS <- netAnalyze(net_IBS, 
                            centrLCC = FALSE,
                            avDissIgnoreInf = TRUE,
                            sPathNorm = FALSE,
                            hubPar = c("degree", "eigenvector"),
                            hubQuant = 0.9,
                            lnormFit = TRUE,
                            normDeg = FALSE,
                            normBetw = FALSE,
                            normClose = FALSE,
                            normEigen = FALSE,
                            gcmHeat = FALSE)
    
    # Saving the network properties
    saveRDS(summary(props_IBS), file.path(path.properties,  agg_level, 
                                           paste0("Summary_",str_to_title(physeq_name),"_",
                                                  method[m],".rds")))
    
    
    plot_names <- paste(str_to_title(physeq_name), variable_names, method[m], sep=" ")
      
    png(file.path(path.plots, agg_level, paste0("Network_comparison_",physeq_name,"_",method[m],".png")),
               width=1600, height=700)
    
    plot(props_IBS, 
         sameLayout = TRUE, 
         layoutGroup = "union",
         rmSingles = "inboth",
         repulsion = 0.7,
         nodeColor = "cluster",
         nodeSize = "eigenvector",
         labelScale = FALSE,
         cexNodes = 1, 
         cexLabels = 1,
         cexHubLabels = 1,
         cexTitle = 2,
         groupNames = plot_names,
         #title1 = paste0(str_to_title(physeq_name), " Network - Method: ", method[m]), 
         #showTitle = TRUE,
         hubBorderCol  = "gray40")
    
    dev.off()
    
  }
  
}

get_assoc.matrix_slr <- function(otutab){
  
  dim <- ncol(otutab) - 4
  
  #Definig ranks
  if(dim>64){
    ranks <- round(seq(2, 64, len=10))
  }else{
    ranks <- round(seq(2, dim, len=8))
  }
  
  se.slr <- spiec.easi(otutab , method='slr', nlambda=50,
                       lambda.min.ratio=1e-2, r=ranks, lambda.log=TRUE,
                       pulsar.params=list(thresh=0.1, ncores=20, rep.num=30))
  
  se.slr$ebic <- sapply(se.slr, function(x)
    ebic(x$refit$stars, x$est$data,
         x$est$loglik[x$select$stars$opt.index]))
  
  opt <- getOptSLR(se.slr)
  
  icov <- Matrix::drop0(opt$est$icov[[getOptInd(opt)]])
  secor <- cov2cor(prec2cov(icov))
  assoc.matrix <- as.matrix(secor*getRefit(opt))
  
  return(assoc.matrix)
}

network_comparison_slr <- function(physeq_name, agg_level, filtTax, filtTaxPar, filtSamp, filtSampPar,
                                       method, self.criterion ="stars",  assoMat = T){
  
  # Read the phyloseq object
  phyloseq <- readRDS(file.path(path.phylobj, agg_level, paste0("agglo_",physeq_name,".rds")))
  
  # Renaming ASVs for agg_level name taxonomy
  otutab <- as.data.frame(otu_table(phyloseq))
  
  if(identical(colnames(otutab), taxa_names(phyloseq))){
    colnames(otutab) <- as.data.frame(tax_table(phyloseq))[[agg_level]]
  }else{
    print("ASVs in OTU Table does not match the ASVs in the Taxonomy Table")
  }
  
  otutab <- otu_table(as.matrix(otutab), taxa_are_rows=FALSE)
  phyloseq@otu_table <- otutab
  taxa <- as.data.frame(phyloseq@tax_table)
  
  rownames(phyloseq@tax_table) <- taxa[,agg_level]
  
  variable_names <- c("IBS", "Healthy")
  
  phyloseq_IBS <- phyloseq::subset_samples(phyloseq, host_disease == "IBS")
  phyloseq_IBS <- prune_taxa(taxa_sums(phyloseq_IBS)>0, phyloseq_IBS)
  
  phyloseq_H <- phyloseq::subset_samples(phyloseq, host_disease == "Healthy")
  phyloseq_H <- prune_taxa(taxa_sums(phyloseq_H)>0, phyloseq_H)
  
  # Filtering
  net_IBS <- netConstruct(data= phyloseq_IBS,
                          data2= phyloseq_H,
                          filtTax =  filtTax,
                          filtTaxPar = filtTaxPar,
                          filtSamp =   filtSamp,
                          filtSampPar =  filtSampPar,
                           measure = NULL,
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3)
  
  assoc.matrix_IBS <- get_assoc.matrix_slr(net_IBS$countMat1)
  assoc.matrix_H <- get_assoc.matrix_slr(net_IBS$countMat2)
  
  rownames(assoc.matrix_IBS) <- colnames(assoc.matrix_IBS) <- colnames(net_IBS$countMat1)
  rownames(assoc.matrix_H) <- colnames(assoc.matrix_H) <- colnames(net_IBS$countMat2)
  
  #Saving association matrix
  write.csv(assoc.matrix_IBS, file = file.path(path.assoc_mat, agg_level,
                                               paste0("AssocMat_",str_to_title(physeq_name), "_" ,variable_names[1],"_slr.csv")), row.names = T)
  write.csv(assoc.matrix_H, file = file.path(path.assoc_mat, agg_level,
                                               paste0("AssocMat_",str_to_title(physeq_name), "_" ,variable_names[2],"_slr.csv")), row.names = T)
  
  net_IBS_r <- netConstruct(data= assoc.matrix_IBS,
                          data2= assoc.matrix_H,
                          dataType = "condDependence",
                          sparsMethod = "none",
                          verbose = 0)
  
  props_IBS <- netAnalyze(net_IBS_r, 
                          centrLCC = FALSE,
                          avDissIgnoreInf = TRUE,
                          sPathNorm = FALSE,
                          clustMethod = "cluster_fast_greedy",
                          hubPar = c("degree", "eigenvector"),
                          hubQuant = 0.9,
                          lnormFit = TRUE,
                          normDeg = FALSE,
                          normBetw = FALSE,
                          normClose = FALSE,
                          normEigen = FALSE,
                          gcmHeat = FALSE)
  
  # Saving the network properties
  saveRDS(summary(props_IBS), file.path(path.properties,  agg_level, 
                                         paste0("Summary_",str_to_title(physeq_name),"_slr.rds")))
  
  plot_names <- paste(str_to_title(physeq_name), variable_names, "slr", sep=" ")
  
  png(file.path(path.plots, agg_level, paste0("Network_comparison_",physeq_name,"_slr.png")),
      width=1600, height=700)
  
  plot(props_IBS, 
       sameLayout = TRUE, 
       layoutGroup = "union",
       rmSingles = "inboth",
       repulsion = 0.7,
       nodeColor = "cluster",
       nodeSize = "eigenvector",
       labelScale = FALSE,
       cexNodes = 1, 
       cexLabels = 1,
       cexHubLabels = 1,
       cexTitle = 2,
       groupNames = plot_names,
       hubBorderCol  = "gray40")
  
  dev.off()
  
  
}

## Network Analysis Summary

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
  png(file.path(path.plots, paste0("Venn_diagram_",physeq_name,".png")), width = 800, height = 600) 
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