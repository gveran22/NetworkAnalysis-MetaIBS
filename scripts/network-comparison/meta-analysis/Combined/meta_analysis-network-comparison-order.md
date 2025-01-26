Meta Analysis - Network Comparison Combined (Order)
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

------------------------------------------------------------------------

# 3. NETWORK INFERENCE

------------------------------------------------------------------------

## 3.1. Total Merge

![](../../../../outputs/network-comparison/Combined/plots/Order/total-merge-1.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/total-merge-2.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/total-merge-3.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/total-merge-4.png)<!-- -->

## 3.2. Sample Type

## 3.3. Sequencing Technology

![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-1.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-2.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-3.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-4.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-5.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-6.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-7.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-8.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-9.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-10.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-11.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/seq-tech-12.png)<!-- -->

## 3.4. Variable Region

![](../../../../outputs/network-comparison/Combined/plots/Order/global-prop-mb-1.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/global-prop-mb-2.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/global-prop-mb-3.png)<!-- -->![](../../../../outputs/network-comparison/Combined/plots/Order/global-prop-mb-4.png)<!-- -->
