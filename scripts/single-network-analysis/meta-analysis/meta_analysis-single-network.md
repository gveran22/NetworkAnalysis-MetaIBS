Your Document Title
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

![](meta_analysis-single-network_files/figure-gfm/meta-analysis-glasso-1.png)<!-- -->

![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-10.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-11.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-1-12.png)<!-- -->

## MB

![](meta_analysis-single-network_files/figure-gfm/meta-analysis-mb-1.png)<!-- -->

![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-7.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-8.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-9.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-10.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-11.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-2-12.png)<!-- -->

## SLR

![](meta_analysis-single-network_files/figure-gfm/meta-analysis-slr-1.png)<!-- -->

![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-5.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-6.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-7.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-8.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-9.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-10.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-11.png)<!-- -->![](meta_analysis-single-network_files/figure-gfm/unnamed-chunk-3-12.png)<!-- -->
