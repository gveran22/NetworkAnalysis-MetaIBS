library(shiny)
library(NetCoMi)
library(stringr)
library(ggplot2)
library(gridExtra)

# Helper: Load matrices from RData files
load_matrices <- function(datasets_names, path.assoc_mat, agg_level, method) {
  object_name <- switch(method,
                        "gl" = "assoMat.gl",
                        "mb" = "assoMat.mb",
                        "slr" = "assoMat.slr",
                        stop("Invalid method selected."))
  
  matrices <- lapply(datasets_names, function(dataset) {
    file_path <- file.path(path.assoc_mat, agg_level, paste0("AssocMat_", dataset, ".RData"))
    
    tryCatch({
      env <- new.env()
      load(file = file_path, envir = env)
      
      if (!exists(object_name, envir = env)) {
        warning(paste("Object", object_name, "not found in", file_path))
        return(NULL)
      }
      
      mat <- get(object_name, envir = env)
      if (!is.matrix(mat)) {
        warning(paste("Object", object_name, "in", file_path, "is not a matrix"))
        return(NULL)
      }
      
      return(mat)
      
    }, error = function(e) {
      message(paste("Skipping", dataset, "due to error:", e$message))
      return(NULL)
    })
  })
  
  names(matrices) <- datasets_names
  Filter(Negate(is.null), matrices)
}

# Helper: Summarize matrices
summary_assoc_matrix <- function(matrices, summary){
  all_rows <- unique(unlist(lapply(matrices, rownames)))
  all_cols <- unique(unlist(lapply(matrices, colnames)))
  full_matrix <- matrix(0, nrow = length(all_rows), ncol = length(all_cols))
  rownames(full_matrix) <- all_rows
  colnames(full_matrix) <- all_cols
  
  if (summary == "proportion") {
    count_matrix <- full_matrix
    for (mat in matrices) {
      temp <- full_matrix
      temp[rownames(mat), colnames(mat)] <- mat
      count_matrix <- count_matrix + (temp != 0)
    }
    final_matrix <- count_matrix / length(matrices)
  } else {
    full_matrices <- lapply(matrices, function(x) {
      temp <- full_matrix
      temp[rownames(x), colnames(x)] <- x
      return(temp)
    })
    sum_matrix <- Reduce("+", full_matrices)
    if (summary == "sum") {
      final_matrix <- sum_matrix
    } else if (summary == "mean") {
      final_matrix <- sum_matrix / length(matrices)
    } else if (summary %in% c("min", "max", "median")) {
      matrix_array <- simplify2array(full_matrices)
      final_matrix <- apply(matrix_array, c(1, 2), match.fun(summary))
    } else {
      stop("Invalid summary method.")
    }
  }
  
  return(final_matrix)
}

# Combine networks from different methods
combine_networks <- function(matrices_by_method, summary = "intersection", threshold = 0) {
  binarized <- lapply(matrices_by_method, function(mat) {
    mat[abs(mat) < threshold] <- 0
    (mat != 0) * 1
  })
  
  dims <- dim(binarized[[1]])
  result <- matrix(0, nrow = dims[1], ncol = dims[2])
  names <- rownames(binarized[[1]])
  rownames(result) <- colnames(result) <- names
  
  if (summary == "intersection") {
    result <- Reduce("*", binarized)
  } else if (summary == "union") {
    result <- Reduce("+", binarized)
    result[result > 0] <- 1
  } else if (summary == "majority") {
    result <- Reduce("+", binarized)
    result[result < 2] <- 0
    result[result >= 2] <- 1
  } else if (summary %in% c("mean", "median")) {
    weight_stack <- simplify2array(matrices_by_method)
    result <- apply(weight_stack, c(1, 2), match.fun(summary))
    result[abs(result) < threshold] <- 0
  }
  
  return(result)
}

# UI
ui <- fluidPage(
  titlePanel("Meta-Network Visualizer (NetCoMi Edition)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("agg_level", "Aggregation Level", choices = c("Order", "Family", "Genus")),
      checkboxGroupInput("methods", "SpiecEasi Method(s)", choices = c("gl", "mb", "slr"), selected = "gl"),
      selectInput("summary", "Summary Method", choices = c("mean", "median", "min", "max", "proportion", "sum")),
      selectInput("combine_summary", "Combine Networks Summary", choices = c("intersection", "union", "majority", "mean", "median")),
      uiOutput("dataset_selector"),
      actionButton("update", "Update Network")
    ),
    mainPanel(
      uiOutput("network_ui"),
      plotOutput("combined_network_plot", height = "600px")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  observeEvent(input$agg_level, {
    path <- "~/MetaIBS/outputs/single-network-analysis/association_matrices"
    files <- list.files(path = file.path(path, input$agg_level), pattern = "^AssocMat_.*\\.RData$")
    datasets <- sort(sub("^AssocMat_(.*)\\.RData$", "\\1", files))
    
    output$dataset_selector <- renderUI({
      checkboxGroupInput("datasets", "Select Dataset(s)", choices = datasets, selected = datasets)
    })
  })
  
  observeEvent(input$update, {
    datasets <- input$datasets
    
    output$network_ui <- renderUI({
      plot_output_list <- lapply(input$methods, function(method) {
        matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, method)
        summarized <- summary_assoc_matrix(matrices, input$summary)
        edge_values <- summarized[upper.tri(summarized)]
        edge_values <- edge_values[edge_values != 0]
        max_val <- ifelse(length(edge_values) > 0, round(max(abs(edge_values)), 2), 0.2)
        
        tagList(
          h4(paste("Method:", method)),
          sliderInput(paste0("threshold_", method), paste("Threshold for", method),
                      min = 0, max = max_val, value = 0.001, step = 0.001),
          downloadButton(paste0("download_plot_", method), "Download Network Plot"),
          plotOutput(paste0("network_plot_", method), height = "600px"),
          plotOutput(paste0("hist_plot_", method), height = "400px")
        )
      })
      do.call(tagList, plot_output_list)
    })
    
    lapply(input$methods, function(method) {
      local({
        m <- method
        
        output[[paste0("network_plot_", m)]] <- renderPlot({
          matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
          if (length(matrices) > 0) {
            summarized <- summary_assoc_matrix(matrices, input$summary)
            threshold <- input[[paste0("threshold_", m)]]
            dataType <- if (input$summary == "proportion") "proportionality" else "condDependence"
            filtered_mat <- summarized
            filtered_mat[abs(filtered_mat) < threshold] <- 0
            props_asso <- network_construct(filtered_mat, dataType = dataType, thresh = 0)
            plot(props_asso,
                 repulsion = 0.7,
                 sameLayout = TRUE,
                 nodeColor = "cluster",
                 labelScale = FALSE,
                 rmSingles = TRUE,
                 nodeSize = "eigenvector",
                 cexNodes = 0.58,
                 cexLabels = 0.7,
                 cexHubLabels = 1,
                 title1 = paste("Meta-Network (", toupper(m), ")", sep = ""),
                 showTitle = TRUE,
                 cexTitle = 2,
                 hubBorderCol = "gray40")
          }
        })
        
        output[[paste0("download_plot_", m)]] <- downloadHandler(
          filename = function() paste0("network_plot_", m, ".png"),
          content = function(file) {
            png(file, width = 800, height = 800)
            matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
            summarized <- summary_assoc_matrix(matrices, input$summary)
            threshold <- input[[paste0("threshold_", m)]]
            dataType <- if (input$summary == "proportion") "proportionality" else "condDependence"
            filtered_mat <- summarized
            filtered_mat[abs(filtered_mat) < threshold] <- 0
            props_asso <- network_construct(filtered_mat, dataType = dataType, thresh = 0)
            plot(props_asso,
                 repulsion = 0.7,
                 sameLayout = TRUE,
                 nodeColor = "cluster",
                 labelScale = FALSE,
                 rmSingles = TRUE,
                 nodeSize = "eigenvector",
                 cexNodes = 0.58,
                 cexLabels = 0.7,
                 cexHubLabels = 1,
                 title1 = paste("Meta-Network (", toupper(m), ")", sep = ""),
                 showTitle = TRUE,
                 cexTitle = 2,
                 hubBorderCol = "gray40")
            dev.off()
          }
        )
        
        output[[paste0("hist_plot_", m)]] <- renderPlot({
          matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
          if (length(matrices) > 0) {
            summarized <- summary_assoc_matrix(matrices, input$summary)
            edge_values <- summarized[upper.tri(summarized)]
            edge_values <- edge_values[edge_values != 0]
            threshold <- input[[paste0("threshold_", m)]]
            
            df <- data.frame(value = edge_values)
            df$included <- abs(df$value) >= threshold
            
            total <- nrow(df)
            kept <- sum(df$included)
            mean_all <- round(mean(df$value), 4)
            sd_all <- round(sd(df$value), 4)
            mean_kept <- round(mean(df$value[df$included]), 4)
            sd_kept <- round(sd(df$value[df$included]), 4)
            percent_kept <- round((kept / total) * 100, 2)
            
            ggplot(df, aes(x = value, fill = included)) +
              geom_histogram(bins = 50, color = "white") +
              scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray"), labels = c("FALSE" = "NO", "TRUE" = "YES")) +
              labs(title = paste("Edge Weights (", toupper(m), ")", sep = ""), x = "Edge Weight", y = "Frequency", fill = "Included") +
              theme_minimal() +
              annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2, size = 4, fontface = "bold",
                       label = paste0("Total edges: ", total,
                                      "\nTotal mean: ", mean_all,
                                      "\nTotal SD: ", sd_all,
                                      "\nKept edges: ", kept,
                                      "\nKept mean: ", mean_kept,
                                      "\nKept SD: ", sd_kept,
                                      "\nKept edges (%): ", percent_kept, "%"))
          }
        })
      })
    })
    
    output$combined_network_plot <- renderPlot({
      method_matrices <- lapply(input$methods, function(m) {
        mats <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
        summarized <- summary_assoc_matrix(mats, input$summary)
        threshold <- input[[paste0("threshold_", m)]]
        summarized[abs(summarized) < threshold] <- 0
        summarized
      })
      
      names(method_matrices) <- input$methods
      common_mat <- combine_networks(method_matrices, summary = input$combine_summary, threshold = 0)
      net <- network_construct(common_mat, dataType = if (input$summary == "proportion") "proportionality" else "condDependence", thresh = 0)
      
      plot(net,
           title1 = paste("Combined Network (", input$combine_summary, ")"),
           showTitle = TRUE,
           cexTitle = 2,
           nodeColor = "cluster",
           sameLayout = TRUE,
           cexNodes = 0.6,
           cexLabels = 3,
           repulsion = 0.8,
           rmSingles = TRUE)
    })
  })
}

shinyApp(ui, server)
