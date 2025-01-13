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

# 2. META-ANALYSIS

------------------------------------------------------------------------

## GLasso

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-glasso-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-glasso-12.png)<!-- -->
\### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.65000 | 0.94000 | 0.71000 | 0.65000 | 0.92000 | 0.82000 | 0.94000 | 0.94000 | 0.92000 | 0.93000 | 0.82000 | 0.81000 |
| Clustering coefficient | 0.64890 | 0.29814 | 0.58751 | 0.42134 | 0.42837 | 0.43542 | 0.36925 | 0.27728 | 0.43812 | 0.27395 | 0.48746 | 0.34302 |
| Modularity | 0.41848 | 0.41257 | 0.39452 | 0.72024 | 0.42520 | 0.70008 | 0.39776 | 0.57840 | 0.37074 | 0.38257 | 0.40451 | 0.59542 |
| Positive edge percentage | 80.16194 | 64.63023 | 67.90541 | 99.06542 | 68.52792 | 96.89922 | 66.56977 | 63.22870 | 62.89855 | 60.32258 | 55.59441 | 66.27907 |
| Edge density | 0.11875 | 0.07115 | 0.11911 | 0.05144 | 0.09412 | 0.03884 | 0.07870 | 0.05102 | 0.08242 | 0.07246 | 0.08612 | 0.05309 |
| Natural connectivity | 0.03145 | 0.01608 | 0.02756 | 0.02040 | 0.02335 | 0.01565 | 0.02334 | 0.01436 | 0.01901 | 0.01746 | 0.02302 | 0.01650 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.68352 | 0.69551 | 0.68206 | 0.65369 | 0.68094 | 0.65599 | 0.68062 | 0.69405 | 0.69535 | 0.69126 | 0.69700 | 0.68991 |
| Average path length\*\* | 1.86587 | 2.09263 | 1.76275 | 3.86101 | 1.87126 | 3.58853 | 2.10940 | 2.53334 | 2.08925 | 2.03216 | 2.36088 | 2.70845 |

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-12.png)<!-- -->
\### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.90000 | 0.97000 | 0.85000 | 0.98000 | 1.00000 | 0.99000 | 1.00000 | 1.00000 | 0.96000 | 1.00000 | 1.00000 | 0.98000 |
| Clustering coefficient | 0.46678 | 0.16828 | 0.41091 | 0.09816 | 0.21835 | 0.06348 | 0.16687 | 0.07563 | 0.23047 | 0.12643 | 0.14276 | 0.07865 |
| Modularity | 0.55495 | 0.45001 | 0.51908 | 0.62218 | 0.42634 | 0.56652 | 0.46210 | 0.47178 | 0.43552 | 0.38646 | 0.48541 | 0.54307 |
| Positive edge percentage | 86.32812 | 65.77947 | 76.54867 | 88.60759 | 68.38906 | 89.57055 | 67.27941 | 57.46606 | 67.63636 | 60.40268 | 61.57205 | 57.86517 |
| Edge density | 0.06392 | 0.05649 | 0.06331 | 0.03324 | 0.06646 | 0.03360 | 0.05495 | 0.04465 | 0.06031 | 0.06020 | 0.04626 | 0.03745 |
| Natural connectivity | 0.01623 | 0.01378 | 0.01667 | 0.01275 | 0.01443 | 0.01277 | 0.01345 | 0.01273 | 0.01459 | 0.01372 | 0.01300 | 0.01267 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.68829 | 0.69636 | 0.69015 | 0.66116 | 0.69154 | 0.65445 | 0.69286 | 0.69192 | 0.69493 | 0.69475 | 0.68510 | 0.68982 |
| Average path length\*\* | 2.35202 | 2.10650 | 2.25775 | 3.07428 | 1.92556 | 2.76511 | 2.08351 | 2.27104 | 2.20365 | 1.99758 | 2.41886 | 2.67187 |

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-slr-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-12.png)<!-- -->

### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.11000 | 0.79000 | 0.10000 | 0.72000 | 0.84000 | 0.99000 | 0.81000 | 0.87000 | 0.82000 | 0.87000 | 0.77000 | 0.76000 |
| Clustering coefficient | 0.00000 | 0.10746 | 0.00000 | 0.19578 | 0.16058 | 0.23643 | 0.18249 | 0.29036 | 0.36826 | 0.13746 | 0.35639 | 0.17363 |
| Modularity | 0.43500 | 0.47809 | 0.35500 | 0.31344 | 0.36716 | 0.45163 | 0.46772 | 0.47350 | 0.43014 | 0.46718 | 0.38426 | 0.37106 |
| Positive edge percentage | 90.00000 | 51.00671 | 60.00000 | 49.77376 | 48.52321 | 78.32700 | 43.60465 | 55.07246 | 61.11111 | 50.85714 | 52.01613 | 47.69231 |
| Edge density | 0.18182 | 0.04836 | 0.22222 | 0.08646 | 0.06799 | 0.05422 | 0.05309 | 0.05533 | 0.07588 | 0.04678 | 0.08476 | 0.06842 |
| Natural connectivity | 0.12695 | 0.01564 | 0.14048 | 0.01959 | 0.01601 | 0.01368 | 0.01556 | 0.01505 | 0.01748 | 0.01445 | 0.01887 | 0.01738 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.68159 | 0.70660 | 0.71355 | 0.70537 | 0.70308 | 0.68321 | 0.70697 | 0.70196 | 0.70117 | 0.69837 | 0.70307 | 0.70817 |
| Average path length\*\* | 2.38667 | 2.37080 | 1.91672 | 1.95806 | 1.91386 | 2.08453 | 2.28580 | 2.14900 | 1.84978 | 2.31050 | 1.87344 | 1.93791 |
