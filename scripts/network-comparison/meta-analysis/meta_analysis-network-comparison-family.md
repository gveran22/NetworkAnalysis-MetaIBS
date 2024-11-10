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

![](images/family/meta-analysis-glasso-1.png)<!-- -->

![](images/family/single-network-glasso-1.png)<!-- -->![](images/family/single-network-glasso-2.png)<!-- -->![](images/family/single-network-glasso-3.png)<!-- -->![](images/family/single-network-glasso-4.png)<!-- -->![](images/family/single-network-glasso-5.png)<!-- -->![](images/family/single-network-glasso-6.png)<!-- -->![](images/family/single-network-glasso-7.png)<!-- -->![](images/family/single-network-glasso-8.png)<!-- -->![](images/family/single-network-glasso-9.png)<!-- -->![](images/family/single-network-glasso-10.png)<!-- -->![](images/family/single-network-glasso-11.png)<!-- -->![](images/family/single-network-glasso-12.png)<!-- -->

## MB

![](images/family/meta-analysis-mb-1.png)<!-- -->

![](images/family/single-network-mb-1.png)<!-- -->![](images/family/single-network-mb-2.png)<!-- -->![](images/family/single-network-mb-3.png)<!-- -->![](images/family/single-network-mb-4.png)<!-- -->![](images/family/single-network-mb-5.png)<!-- -->![](images/family/single-network-mb-6.png)<!-- -->![](images/family/single-network-mb-7.png)<!-- -->![](images/family/single-network-mb-8.png)<!-- -->![](images/family/single-network-mb-9.png)<!-- -->![](images/family/single-network-mb-10.png)<!-- -->![](images/family/single-network-mb-11.png)<!-- -->![](images/family/single-network-mb-12.png)<!-- -->

## SLR

![](images/family/meta-analysis-slr-1.png)<!-- -->

![](images/family/single-network-slr-1.png)<!-- -->![](images/family/single-network-slr-2.png)<!-- -->![](images/family/single-network-slr-3.png)<!-- -->![](images/family/single-network-slr-4.png)<!-- -->![](images/family/single-network-slr-5.png)<!-- -->![](images/family/single-network-slr-6.png)<!-- -->![](images/family/single-network-slr-7.png)<!-- -->![](images/family/single-network-slr-8.png)<!-- -->![](images/family/single-network-slr-9.png)<!-- -->![](images/family/single-network-slr-10.png)<!-- -->![](images/family/single-network-slr-11.png)<!-- -->![](images/family/single-network-slr-12.png)<!-- -->
