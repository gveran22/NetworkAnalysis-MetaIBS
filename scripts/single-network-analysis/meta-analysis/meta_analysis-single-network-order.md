Meta analysis - Single Network Anaylsis (Order)
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

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-glasso-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-glasso-12.png)<!-- -->
\### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.38824 | 0.60526 | 0.40426 | 0.20000 | 0.71429 | 0.21429 | 0.78571 | 0.45833 | 0.58491 | 0.63889 | 0.33333 | 0.50000 |
| Clustering coefficient | 0.77511 | 0.57849 | 0.74573 | 1.00000 | 0.52154 | 0.00000 | 0.41143 | 0.61989 | 0.77758 | 0.53221 | 0.60306 | 0.25146 |
| Modularity | 0.15181 | 0.34664 | 0.15695 | -0.22222 | 0.45132 | -0.12500 | 0.48996 | 0.18364 | 0.12278 | 0.21941 | 0.25000 | 0.31500 |
| Positive edge percentage | 98.49057 | 73.91304 | 91.80328 | 33.33333 | 83.60000 | 0.00000 | 73.21429 | 33.33333 | 67.15328 | 50.76923 | 50.00000 | 60.00000 |
| Edge density | 0.50189 | 0.18182 | 0.35673 | 1.00000 | 0.10352 | 0.66667 | 0.10606 | 0.32727 | 0.29462 | 0.25692 | 0.38095 | 0.22222 |
| Natural connectivity | 0.14279 | 0.06500 | 0.10588 | 0.57234 | 0.02504 | 0.55295 | 0.04009 | 0.13611 | 0.07795 | 0.07008 | 0.21850 | 0.14169 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 2.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 2.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.65889 | 0.67686 | 0.66478 | 0.71351 | 0.68399 | 0.72082 | 0.69358 | 0.70467 | 0.68535 | 0.70044 | 0.70726 | 0.70178 |
| Average path length\*\* | 1.11920 | 1.70829 | 1.24855 | 0.71351 | 1.91172 | 0.96110 | 2.14258 | 1.37969 | 1.27308 | 1.57276 | 1.23969 | 1.66335 |

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-12.png)<!-- -->

### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.90588 | 0.76316 | 0.74468 | 0.20000 | 0.98980 | 0.21429 | 0.97619 | 0.33333 | 0.69811 | 0.66667 | 0.52381 | 0.55000 |
| Clustering coefficient | 0.36240 | 0.31582 | 0.45496 | 0.00000 | 0.18828 | 0.00000 | 0.15389 | 0.40625 | 0.49871 | 0.42567 | 0.00000 | 0.24527 |
| Modularity | 0.55533 | 0.48334 | 0.52401 | -0.12500 | 0.48561 | -0.12500 | 0.52805 | 0.21605 | 0.29931 | 0.40448 | 0.41500 | 0.33471 |
| Positive edge percentage | 89.50617 | 73.17073 | 83.33333 | 50.00000 | 79.79094 | 0.00000 | 67.69231 | 44.44444 | 70.88608 | 57.44681 | 40.00000 | 54.54545 |
| Edge density | 0.05537 | 0.10099 | 0.09076 | 0.66667 | 0.06164 | 0.66667 | 0.07927 | 0.32143 | 0.11862 | 0.17029 | 0.18182 | 0.20000 |
| Natural connectivity | 0.01762 | 0.04501 | 0.03760 | 0.55386 | 0.01458 | 0.55206 | 0.03110 | 0.18654 | 0.03829 | 0.05855 | 0.12514 | 0.12691 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.67832 | 0.68773 | 0.67804 | 0.71786 | 0.67657 | 0.72565 | 0.69183 | 0.70500 | 0.69257 | 0.69627 | 0.70975 | 0.70101 |
| Average path length\*\* | 2.40812 | 2.39279 | 2.21964 | 0.95715 | 2.07846 | 0.96753 | 2.25378 | 1.45850 | 1.73131 | 1.81342 | 2.05593 | 1.68932 |

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-12.png)<!-- -->

### Global Properties

|  |  |  |  |  |  |  |  |  |  |  |  |  |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Relative LCC size | 0.52941 | 0.50000 | 0.51064 | 0.13333 | 0.65306 | 0.21429 | 0.38095 | 0.12500 | 0.15094 | 0.58333 | 0.14286 | 0.10000 |
| Clustering coefficient | 0.48308 | 0.22729 | 0.46050 | 0.00000 | 0.39911 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.16079 | 0.00000 | 0.00000 |
| Modularity | 0.29152 | 0.26698 | 0.26395 | 0.00000 | 0.20225 | -0.12500 | 0.49846 | -0.12500 | 0.33673 | 0.40496 | -0.12500 | 0.00000 |
| Positive edge percentage | 47.71574 | 58.33333 | 67.64706 | 100.00000 | 45.30744 | 0.00000 | 77.77778 | 50.00000 | 42.85714 | 54.54545 | 50.00000 | 0.00000 |
| Edge density | 0.19899 | 0.21053 | 0.24638 | 1.00000 | 0.15327 | 0.66667 | 0.15000 | 0.66667 | 0.25000 | 0.15714 | 0.66667 | 1.00000 |
| Natural connectivity | 0.04004 | 0.07369 | 0.06337 | 0.79796 | 0.03342 | 0.55253 | 0.08303 | 0.55512 | 0.18056 | 0.06338 | 0.55557 | 0.79016 |
| Vertex connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Edge connectivity | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| Average dissimilarity\* | 0.70164 | 0.70624 | 0.69840 | 0.70529 | 0.70737 | 0.72282 | 0.70261 | 0.71082 | 0.71882 | 0.70762 | 0.70725 | 0.74320 |
| Average path length\*\* | 1.45537 | 1.55168 | 1.34867 | 0.70529 | 1.77981 | 0.96376 | 2.20720 | 0.94776 | 1.75085 | 1.86486 | 0.94299 | 0.74320 |
