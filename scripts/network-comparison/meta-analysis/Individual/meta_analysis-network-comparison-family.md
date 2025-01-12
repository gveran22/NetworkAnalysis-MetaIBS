Meta Analysis - Network Comparison (Family)
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

# 2. META-ANALYSIS

------------------------------------------------------------------------

## GLasso

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/meta-analysis-glasso-1.png)<!-- -->

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-1.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-2.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-3.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-4.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-5.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-6.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-7.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-8.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-9.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-10.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-11.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-glasso-12.png)<!-- -->

## MB

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/meta-analysis-mb-1.png)<!-- -->

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-1.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-2.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-3.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-4.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-5.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-6.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-7.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-8.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-9.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-10.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-11.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-mb-12.png)<!-- -->

## SLR

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/meta-analysis-slr-1.png)<!-- -->

![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-1.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-2.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-3.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-4.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-5.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-6.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-7.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-8.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-9.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-10.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-11.png)<!-- -->![](~/MetaIBS/outputs/network-comparison/Individual/plots/Family/single-network-slr-12.png)<!-- -->
