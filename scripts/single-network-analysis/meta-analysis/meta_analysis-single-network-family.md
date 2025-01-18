Meta analysis - Single Network Anaylsis
================
2024-03-01

``` r
myPaths <- .libPaths()
myPaths <- c(myPaths, "~/MetaIBS/MetaIBS-library")
myPaths <- c(myPaths[3], myPaths[1], myPaths[2])
.libPaths(myPaths)  # add new path
```

------------------------------------------------------------------------

# 1. IMPORT

------------------------------------------------------------------------

## 1.1. Libraries

``` r
library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
library(tidyverse)
library(ggplot2)
library(SpiecEasi)
library(igraph)
library(VennDiagram)
library(NetCoMi)
```

------------------------------------------------------------------------

# 2. FUNCTIONS

------------------------------------------------------------------------

``` r
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
plot_meta_analysis <- function(meta, layout = NULL, title = "Meta-Analysis", repulsion = 0.5) {
  props_asso_meta <- network_construct(meta, thresh = 0.1)
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
       cexTitle = 2,
       hubBorderCol = "gray40")
}

# Function to process individual plots
plot_individual_networks <- function(matrices, meta, layout, datasets_names) {
  for (i in seq_along(matrices)) {
    sub_adj <- matrices[[i]]
    
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
         cexLabels = 1.5,
         cexHubLabels = 1.7,
         title1 = paste0("Network Analysis ", str_to_title(datasets_names)),
         showTitle = TRUE,
         cexTitle = 2,
         hubBorderCol  = "gray40")
  }
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
```

------------------------------------------------------------------------

# 3. META-ANALYSIS

------------------------------------------------------------------------

## Meta-Analysis Plot

<img src="../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-1.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-2.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-3.png" width="33%" />

<img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-1.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-2.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-3.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-4.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-5.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-6.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-7.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-8.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-9.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-10.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-11.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-12.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-13.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-14.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-15.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-16.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-17.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-18.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-19.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-20.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-21.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-22.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-23.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-24.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-25.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-26.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-27.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-28.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-29.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-30.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-31.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-32.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-33.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-34.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-35.png" width="33%" /><img src="../../../outputs/single-network-analysis/Individual/plots/Family/single-network-36.png" width="33%" />
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* \# 2. META-ANALYSIS
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

## GLasso

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-glasso-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.53459 | 0.81818 | 0.68182 | 0.37500 | 0.84118 | 0.43478 | 0.80822 | 0.39535 | 0.65686 | 0.80597 | 0.42857 | 0.51515 |
| Clustering coefficient | 0.72151 | 0.47158 | 0.65134 | 0.00000 | 0.44557 | 0.00000 | 0.47257 | 0.45473 | 0.66895 | 0.31602 | 0.60622 | 0.42567 |
| Modularity | 0.23516 | 0.39872 | 0.26669 | 0.38281 | 0.41860 | 0.36000 | 0.40873 | 0.32609 | 0.16941 | 0.35995 | 0.30547 | 0.39683 |
| Positive edge percentage | 96.42276 | 68.83117 | 87.22628 | 37.50000 | 85.54348 | 20.00000 | 70.55556 | 52.17391 | 69.48229 | 61.71875 | 42.30769 | 52.38095 |
| Edge density | 0.17227 | 0.10762 | 0.15480 | 0.22222 | 0.09061 | 0.22222 | 0.10520 | 0.16912 | 0.16599 | 0.08945 | 0.24762 | 0.15441 |
| Natural connectivity | 0.07167 | 0.02938 | 0.04965 | 0.15807 | 0.02120 | 0.14100 | 0.02683 | 0.07977 | 0.05268 | 0.02643 | 0.09572 | 0.07877 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.66243 | 0.69249 | 0.67190 | 0.70820 | 0.66749 | 0.70637 | 0.69172 | 0.70271 | 0.68799 | 0.69386 | 0.70637 | 0.70302 |
| Average path length\*\* | 1.72017 | 2.10686 | 1.79831 | 1.79979 | 1.77729 | 1.83433 | 1.96717 | 1.82654 | 1.64028 | 2.10141 | 1.58844 | 1.93466 |

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.96226 | 0.95455 | 0.85227 | 0.20833 | 1.00000 | 0.43478 | 0.98630 | 0.67442 | 0.92157 | 0.98507 | 0.51429 | 0.60606 |
| Clustering coefficient | 0.34916 | 0.28234 | 0.43695 | 0.00000 | 0.16951 | 0.00000 | 0.19906 | 0.21936 | 0.38973 | 0.12133 | 0.18748 | 0.26034 |
| Modularity | 0.48633 | 0.44298 | 0.47806 | 0.21875 | 0.38444 | 0.36000 | 0.46660 | 0.50475 | 0.37543 | 0.44202 | 0.48125 | 0.50473 |
| Positive edge percentage | 89.53488 | 66.43357 | 87.25490 | 50.00000 | 75.91888 | 20.00000 | 68.04734 | 45.94595 | 72.99270 | 64.80000 | 40.00000 | 47.82609 |
| Edge density | 0.04438 | 0.07322 | 0.07351 | 0.40000 | 0.05493 | 0.22222 | 0.06612 | 0.09113 | 0.06269 | 0.05828 | 0.13072 | 0.12105 |
| Natural connectivity | 0.01094 | 0.02129 | 0.01968 | 0.31583 | 0.00943 | 0.14125 | 0.01849 | 0.04363 | 0.01631 | 0.01945 | 0.07253 | 0.06498 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 2.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 2.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.68578 | 0.69452 | 0.68345 | 0.70226 | 0.68176 | 0.70374 | 0.69091 | 0.70300 | 0.69568 | 0.68762 | 0.70802 | 0.70233 |
| Average path length\*\* | 2.22244 | 2.17494 | 1.96940 | 1.39736 | 1.80289 | 1.82596 | 2.12076 | 2.57994 | 2.16928 | 2.26632 | 2.15218 | 2.14937 |

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.57233 | 0.54545 | 0.56818 | 0.08333 | 0.71176 | 0.13043 | 0.67123 | 0.53488 | 0.61765 | 0.56716 | 0.08571 | 0.09091 |
| Clustering coefficient | 0.45060 | 0.36622 | 0.41125 | 0.00000 | 0.48331 | 0.00000 | 0.28166 | 0.30375 | 0.26458 | 0.34161 | 0.00000 | 0.00000 |
| Modularity | 0.26381 | 0.26344 | 0.33357 | 0.00000 | 0.15634 | -0.12500 | 0.39906 | 0.50244 | 0.33886 | 0.40920 | -0.12500 | -0.12500 |
| Positive edge percentage | 45.57957 | 50.83333 | 43.72093 | 100.00000 | 41.59021 | 50.00000 | 50.00000 | 56.25000 | 33.52941 | 55.55556 | 50.00000 | 50.00000 |
| Edge density | 0.12430 | 0.19048 | 0.17551 | 1.00000 | 0.13512 | 0.66667 | 0.10034 | 0.12648 | 0.08705 | 0.12802 | 0.66667 | 0.66667 |
| Natural connectivity | 0.02474 | 0.04215 | 0.03398 | 0.80274 | 0.04203 | 0.55598 | 0.02746 | 0.05677 | 0.02164 | 0.03569 | 0.55358 | 0.55448 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.70065 | 0.70622 | 0.70079 | 0.68419 | 0.70687 | 0.70516 | 0.70380 | 0.70698 | 0.70740 | 0.70710 | 0.71809 | 0.71351 |
| Average path length\*\* | 1.66882 | 1.46322 | 1.53411 | 0.68419 | 1.76460 | 0.94021 | 1.90551 | 2.18468 | 1.80058 | 1.74047 | 0.95746 | 0.95135 |
