# app.R

library(shiny)
library(NetCoMi)
library(ggplot2)
library(stringr)
source("modules/helpers.R")

# --- UI ---

ui <- fluidPage(
  titlePanel("MetaIBS Network Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "analysis_type",
        tabPanel("Single Analysis",
                 selectInput("agg_level", "Aggregation Level", choices = c("Order", "Family", "Genus")),
                 uiOutput("dataset_selector"),
                 checkboxGroupInput("methods", "Select Method(s)", choices = c("gl", "mb", "slr"), selected = c("gl")),
                 selectInput("summary", "Summary Type", choices = c("mean", "median", "proportion", "sum", "min", "max")),
                 selectInput("combine_summary", "Combine Networks Summary", choices = c("(none)", "intersection", "union", "majority", "mean", "median")),
                 actionButton("update", "Update")
        ),
        tabPanel("Comparison",
                 selectInput("comparison_mode", "Comparison Mode", choices = c("Meta-Network", "Group Plots")),
                 conditionalPanel(
                   condition = "input.comparison_mode == 'Meta-Network'",
                   selectInput("agg_level_comp", "Aggregation Level", choices = c("Order", "Family", "Genus")),
                   uiOutput("dataset_selector_comp"),
                   checkboxGroupInput("methods_comp", "Select Method(s)", choices = c("gl", "mb", "slr"), selected = c("gl")),
                   selectInput("summary_comp", "Summary Type", choices = c("mean", "median", "proportion", "sum", "min", "max")),
                   radioButtons("diff_net_comp", "Differential Network", choices = c("No", "Yes"), selected = "No"),
                   actionButton("update_comp", "Update Meta-Network")
                 ),
                 conditionalPanel(
                   condition = "input.comparison_mode == 'Group Plots'",
                   selectInput("agg_level_group", "Aggregation Level", choices = c("Order", "Family", "Genus")),
                   checkboxGroupInput("methods_group", "Select Method(s)", choices = c("gl", "mb", "slr"), selected = c("gl")),
                   selectInput("comparison_variable", "Comparison Variable", choices = c("all", "sample_type", "sequencing_tech", "variable_region")),
                   radioButtons("diff_net_group", "Differential Network", choices = c("No", "Yes"), selected = "No"),
                   actionButton("update_group", "Update Group Plots")
                 )
        )
      )
    ),
    mainPanel(
      uiOutput("network_ui"),
      
      conditionalPanel(
        condition = "input.analysis_type == 'Single Analysis'",
        plotOutput("combined_network_plot", height = "600px")
      ),
      
      uiOutput("comparison_ui"),
      uiOutput("group_ui")
    )
  )
)

# --- Server ---

server <- function(input, output, session) {
  
  # === Populate dataset selector for Single Analysis ===
  observeEvent(input$agg_level, {
    path <- "~/MetaIBS/outputs/single-network-analysis/association_matrices"
    files <- list.files(path = file.path(path, input$agg_level), pattern = "^AssocMat_.*\\.RData$")
    datasets <- sort(sub("^AssocMat_(.*)\\.RData$", "\\1", files))
    
    output$dataset_selector <- renderUI({
      checkboxGroupInput("datasets", "Select Dataset(s)", choices = datasets, selected = datasets)
    })
  })
  
  # === Single Analysis Update ===
  observeEvent(input$update, {
    
    output$comparison_ui <- renderUI(NULL)
    output$group_ui <- renderUI(NULL)
    
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
    
    # Render plots per method
    lapply(input$methods, function(m) {
      output[[paste0("network_plot_", m)]] <- renderPlot({
        matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
        summarized <- summary_assoc_matrix(matrices, input$summary)
        threshold <- input[[paste0("threshold_", m)]]
        dataType <- if (input$summary == "proportion") "proportionality" else "condDependence"
        network <- network_construct(summarized, dataType=dataType, thresh = threshold)
        plot(network, 
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
      })
      
      output[[paste0("hist_plot_", m)]] <- renderPlot({
        matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
        summarized <- summary_assoc_matrix(matrices, input$summary)
        edge_values <- summarized[upper.tri(summarized)]
        edge_values <- edge_values[edge_values != 0]
        threshold <- input[[paste0("threshold_", m)]]
        dataType <- if (input$summary == "proportion") "proportionality" else "condDependence"
        
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
      })
      
      output[[paste0("download_plot_", m)]] <- downloadHandler(
        filename = function() paste0("network_plot_", m, ".png"),
        content = function(file) {
          png(file, width = 800, height = 800)
          matrices <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
          summarized <- summary_assoc_matrix(matrices, input$summary)
          threshold <- input[[paste0("threshold_", m)]]
          dataType <- if (input$summary == "proportion") "proportionality" else "condDependence"
          network <- network_construct(summarized, dataType = dataType ,thresh = threshold)
          plot(network, 
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
    })
    
    output$combined_network_plot <- renderPlot({
      if (input$combine_summary == "(none)") return(NULL)
      
      method_matrices <- lapply(input$methods, function(m) {
      mats <- load_matrices(datasets, "~/MetaIBS/outputs/single-network-analysis/association_matrices", input$agg_level, m)
      summarized <- summary_assoc_matrix(mats, input$summary)
      threshold <- input[[paste0("threshold_", m)]]
      summarized[abs(summarized) < threshold] <- 0
      summarized
      })
      
      names(method_matrices) <- input$methods
      common_mat <- combine_networks(method_matrices, summary = input$combine_summary)
      net <- network_construct(common_mat, thresh = 0)
      
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
  
  # === Populate dataset selector for Meta-Network Comparison ===
  observeEvent(input$agg_level_comp, {
    path <- "~/MetaIBS/outputs/network-comparison/Individual/association_matrices"
    files <- list.files(path = file.path(path, input$agg_level_comp), pattern = "^AssocMat_.*\\.RData$")
    datasets <- sort(sub("^AssocMat_(.*)\\.RData$", "\\1", files))
    
    output$dataset_selector_comp <- renderUI({
      checkboxGroupInput("datasets_comp", "Select Dataset(s)", choices = datasets, selected = datasets)
    })
  })
  
  # === Meta-Network Comparison Update ===
  observeEvent(input$update_comp, {
    output$network_ui <- renderUI(NULL)
    output$group_ui <- renderUI(NULL)
    
    output$comparison_ui <- renderUI({
      plot_output_list <- lapply(input$methods_comp, function(method) {
        group1_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, method, group = "IBS")
        group2_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, method, group = "H")
        if (length(group1_matrices) > 0 && length(group2_matrices) > 0) {
          summarized1 <- summary_assoc_matrix(group1_matrices, input$summary_comp)
          summarized2 <- summary_assoc_matrix(group2_matrices, input$summary_comp)
          edge_values <- c(summarized1[upper.tri(summarized1)], summarized2[upper.tri(summarized2)])
          edge_values <- edge_values[edge_values != 0]
          max_val <- ifelse(length(edge_values) > 0, round(max(abs(edge_values)), 2), 0.2)
          
        output_list <- list(
          h4(paste("Meta-Network for Method:", method)),
          sliderInput(paste0("threshold_comp_", method), paste("Threshold for", method),
                      min = 0, max = max_val, value = 0.001, step = 0.001),
          downloadButton(paste0("download_plot_comp_", method), "Download Comparison Plot"),
          plotOutput(paste0("comp_plot_", method), height = "600px"),
          plotOutput(paste0("hist_plot_ibs_", method), height = "300px"),
          plotOutput(paste0("hist_plot_h_", method), height = "300px")
        )
        
        # Conditional differential 
        if (input$diff_net_comp == "Yes") {
          output_list <- append(output_list, list(
            plotOutput(paste0("diff_plot_comp_", method), height = "600px")
          ))
        }
        return(tagList(output_list))
      }
      })
      do.call(tagList, plot_output_list)
    })
    
    
    lapply(input$methods_comp, function(m) {
      output[[paste0("comp_plot_", m)]] <- renderPlot({
        group1_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "IBS")
        group2_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "H")
        if (length(group1_matrices) > 0 && length(group2_matrices) > 0) {
          summarized1 <- summary_assoc_matrix(group1_matrices, input$summary_comp)
          summarized2 <- summary_assoc_matrix(group2_matrices, input$summary_comp)
          dataType <- if (input$summary_comp == "proportion") "proportionality" else "condDependence"
          threshold <- input[[paste0("threshold_comp_", m)]]
          net <- network_construct_comparison(summarized1, summarized2,dataType = dataType, thresh = threshold)
          plot(net,
               repulsion = 0.7,
               sameLayout = TRUE,
               nodeColor = "cluster",
               labelScale = FALSE,
               rmSingles = "inboth",
               nodeSize = "eigenvector",
               cexNodes = 0.58,
               cexLabels = 0.7,
               cexHubLabels = 1,
               groupNames = c("IBS", "H"),
               showTitle = TRUE,
               cexTitle = 2,
               hubBorderCol = "gray40")
          }
      })
      
      output[[paste0("hist_plot_ibs_", m)]] <- renderPlot({
        group1_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "IBS")
        if (length(group1_matrices) > 0) {
          summarized1 <- summary_assoc_matrix(group1_matrices, input$summary_comp)
          edge_values <- summarized1[upper.tri(summarized1)]
          edge_values <- edge_values[edge_values != 0]
          threshold <- input[[paste0("threshold_comp_", m)]]
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
            theme_minimal() +
            labs(title = paste("IBS Edges (", m, ")", sep = ""), x = "Edge Weight", y = "Frequency") +
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

      output[[paste0("hist_plot_h_", m)]] <- renderPlot({
        group2_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "H")
        if (length(group2_matrices) > 0) {
          summarized2 <- summary_assoc_matrix(group2_matrices, input$summary_comp)
          edge_values <- summarized2[upper.tri(summarized2)]
          edge_values <- edge_values[edge_values != 0]
          threshold <- input[[paste0("threshold_comp_", m)]]
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
            theme_minimal() +
            labs(title = paste("Healthy Edges (", m, ")", sep = ""), x = "Edge Weight", y = "Frequency")+
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
      
      if (input$diff_net_comp == "Yes") {
        output[[paste0("diff_plot_comp_", m)]] <- renderPlot({
          group1_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "IBS")
          group2_matrices <- load_matrices(input$datasets_comp, "~/MetaIBS/outputs/network-comparison/Individual/association_matrices", input$agg_level, m, group = "H")
          if (length(group1_matrices) > 0 && length(group2_matrices) > 0) {
            summarized1 <- summary_assoc_matrix(group1_matrices, input$summary_comp)
            summarized2 <- summary_assoc_matrix(group2_matrices, input$summary_comp)
            threshold <- input[[paste0("threshold_comp_", m)]]
            dataType <- if (input$summary_comp == "proportion") "proportionality" else "condDependence"
            net <- network_construct_comparison(summarized1, summarized2,dataType = dataType, thresh = threshold)
            result <- plot(net,
                 repulsion = 0.7,
                 sameLayout = TRUE,
                 nodeColor = "cluster",
                 labelScale = FALSE,
                 rmSingles = "inboth",
                 nodeSize = "eigenvector",
                 cexNodes = 0.58,
                 cexLabels = 0.7,
                 cexHubLabels = 1,
                 groupNames = c("IBS", "H"),
                 showTitle = TRUE,
                 cexTitle = 2,
                 hubBorderCol = "gray40",
                 doPlot=F)
            
            dif_mat <- diff_networks(net$input$assoMat1, net$input$assoMat2)
            net_dif <- network_construct_comparison(dif_mat$only_1, dif_mat$only_2)
            plot(net_dif,
                 repulsion = 0.7,
                 nodeColor = "cluster",
                 labelScale = FALSE,
                 layout = result$layout$layout1,
                 rmSingles = "none",
                 nodeSize = "eigenvector",
                 cexNodes = 0.58,
                 cexLabels = 0.7,
                 cexHubLabels = 1,
                 groupNames = c("IBS", "H"),
                 showTitle = TRUE,
                 cexTitle = 2,
                 hubBorderCol = "gray40")
          }
            
        })
      }
  
    })
  })
  
  # === Group Plots Update === 
  observeEvent(input$update_group, {
    output$network_ui <- renderUI(NULL)
    output$comparison_ui <- renderUI(NULL)
    
    comp_var <- input$comparison_variable

    # Build path
    path_base <- file.path("~/MetaIBS/outputs/network-comparison/Combined", comp_var, "association_matrices")

    # List all available datasets (just for naming purposes)
    files <- list.files(file.path(path_base, input$agg_level_group), pattern = "^AssocMat_.*\\.RData$")
    groups <- sort(sub("^AssocMat_(.*)\\.RData$", "\\1", files))
    
    output$group_ui <- renderUI({
      req(input$methods_group)
      
      plot_output_list <- lapply(input$methods_group, function(method) {
        method_group_ui <- lapply(groups, function(group_value) {
          
          # Load IBS and H matrices to calculate max edge value
          ibs_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "IBS")[[1]]
          h_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "H")[[1]]
            
          max_val <- 0.2
          if (length(ibs_matrices) > 0 && length(h_matrices) > 0) {
            edge_values <- c(ibs_matrices[upper.tri(ibs_matrices)], h_matrices[upper.tri(h_matrices)])
            edge_values <- edge_values[edge_values != 0]
            max_val <- ifelse(length(edge_values) > 0, round(max(abs(edge_values)), 2), 0.2)
          }
          
          output_list <- list(
            h4(paste("Method:", method, "- Group:", group_value)),
            sliderInput(
              inputId = paste0("threshold_group_", method, "_", group_value),
              label = paste("Threshold for", method, "-", group_value),
              min = 0, max = max_val, value = 0.001, step = 0.001
            ),
            plotOutput(outputId = paste0("group_plot_", method, "_", group_value), height = "600px"),
            plotOutput(outputId = paste0("hist_plot_group_ibs_", method, "_", group_value), height = "300px"),
            plotOutput(outputId = paste0("hist_plot_group_h_", method, "_", group_value), height = "300px")
          )
          
          # Conditional differential 
          if (input$diff_net_group == "Yes") {
            output_list <- append(output_list, list(
              plotOutput(paste0("diff_plot_group_", method, "_", group_value), height = "600px")
            ))
          }
          return(tagList(output_list))
        })
        
        do.call(tagList, method_group_ui)
      })
      
      do.call(tagList, plot_output_list)
    })
    
    
    lapply(input$methods_group, function(method) {
      lapply(groups, function(group_value) {
        
        output[[paste0("group_plot_", method, "_", group_value)]] <- renderPlot({
          ibs_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "IBS")[[1]]
          h_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "H")[[1]]
          
          if (length(ibs_matrices) > 0 && length(h_matrices) > 0) {
            threshold <- input[[paste0("threshold_group_", method, "_", group_value)]]
            net <- network_construct_comparison(ibs_matrices, h_matrices, thresh = threshold)
            plot(net,
                 repulsion = 0.7,
                 sameLayout = TRUE,
                 nodeColor = "cluster",
                 labelScale = FALSE,
                 rmSingles = "inboth",
                 nodeSize = "eigenvector",
                 cexNodes = 0.58,
                 cexLabels = 0.7,
                 cexHubLabels = 1,
                 groupNames = c("IBS", "H"),
                 showTitle = TRUE,
                 cexTitle = 2,
                 hubBorderCol = "gray40")
          }
        })
        
        output[[paste0("hist_plot_group_ibs_", method, "_", group_value)]] <- renderPlot({
          ibs_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "IBS")[[1]]
          if (length(ibs_matrices) > 0) {
            edge_values <- ibs_matrices[upper.tri(ibs_matrices)]
            edge_values <- edge_values[edge_values != 0]
            threshold <- input[[paste0("threshold_group_", method, "_", group_value)]]
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
              theme_minimal() +
              labs(title = paste("IBS Edges (", method, " - ", group_value,")", sep = ""), x = "Edge Weight", y = "Frequency") +
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
        
        output[[paste0("hist_plot_group_h_", method, "_", group_value)]] <- renderPlot({
          h_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "H")[[1]]
          if (length(h_matrices) > 0) {
            edge_values <- h_matrices[upper.tri(h_matrices)]
            edge_values <- edge_values[edge_values != 0]
            threshold <- input[[paste0("threshold_group_", method, "_", group_value)]]
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
              theme_minimal() +
              labs(title = paste("Healthy Edges (", method, " - ", group_value,")", sep = ""), x = "Edge Weight", y = "Frequency")+
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
        
        if (input$diff_net_group == "Yes") {
          output[[paste0("diff_plot_group_", method, "_", group_value)]] <- renderPlot({
            ibs_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "IBS")[[1]]
            h_matrices <- load_matrices(group_value, path_base, input$agg_level_group, method, group = "H")[[1]]
            
            if (length(ibs_matrices) > 0 && length(h_matrices) > 0) {
              threshold <- input[[paste0("threshold_group_", method, "_", group_value)]]
              net <- network_construct_comparison(ibs_matrices, h_matrices, thresh = threshold)
              result <- plot(net,
                             repulsion = 0.7,
                             sameLayout = TRUE,
                             nodeColor = "cluster",
                             labelScale = FALSE,
                             rmSingles = "inboth",
                             nodeSize = "eigenvector",
                             cexNodes = 0.58,
                             cexLabels = 0.7,
                             cexHubLabels = 1,
                             groupNames = c("IBS", "H"),
                             showTitle = TRUE,
                             cexTitle = 2,
                             hubBorderCol = "gray40",
                             doPlot=F)
              
              dif_mat <- diff_networks(net$input$assoMat1, net$input$assoMat2)
              net_dif <- network_construct_comparison(dif_mat$only_1, dif_mat$only_2)
              plot(net_dif,
                   repulsion = 0.7,
                   nodeColor = "cluster",
                   labelScale = FALSE,
                   layout = result$layout$layout1,
                   rmSingles = "none",
                   nodeSize = "eigenvector",
                   cexNodes = 0.58,
                   cexLabels = 0.7,
                   cexHubLabels = 1,
                   groupNames = c("IBS", "H"),
                   showTitle = TRUE,
                   cexTitle = 2,
                   hubBorderCol = "gray40")
            }
            
          })
        }
      })
      
    })    
    
    
  })
}

shinyApp(ui = ui, server = server)
