Meta analysis - Single Network Analysis (Order)
================
2024-03-01

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
library(reshape2)
```

------------------------------------------------------------------------

# 2. FUNCTIONS

------------------------------------------------------------------------

------------------------------------------------------------------------

# 3. META-ANALYSIS

------------------------------------------------------------------------

## 3.1. Summary Plots

<img src="../../../outputs/single-network-analysis/plots/Order/meta-analysis-1.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/meta-analysis-2.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/meta-analysis-3.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/meta-analysis-4.png" width="33%" />

## 3.2. Individual Plots

<img src="../../../outputs/single-network-analysis/plots/Order/single-network-1.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-2.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-3.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-4.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-5.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-6.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-7.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-8.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-9.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-10.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-11.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-12.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-13.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-14.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-15.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-16.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-17.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-18.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-19.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-20.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-21.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-22.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-23.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-24.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-25.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-26.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-27.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-28.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-29.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-30.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-31.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-32.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-33.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-34.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-35.png" width="33%" /><img src="../../../outputs/single-network-analysis/plots/Order/single-network-36.png" width="33%" />

## 3.3. Global Properties

### 3.3.1. Summary

| Study | Method | Rel_LCC | Clustering | Modularity | Positive_Edges | Edge_Density | Natural_Conn | Avg_Path_Length | Hubs |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|:---|
| agp | gl | 0.38824 | 0.77511 | 0.15181 | 98.49057 | 0.50189 | 0.14279 | 1.11920 | 12_Acidimicrobiia(C), 22_Actinobacteria(C), Acetobacterales, Brevibacillales, Bryobacterales, Cardiobacteriales, Frankiales, Propionibacteriales, Syntrophomonadales |
| fukui | gl | 0.60526 | 0.57849 | 0.34664 | 73.91304 | 0.18182 | 0.06500 | 1.70829 | Bacteroidales, Micrococcales, Saccharimonadales, Sphingomonadales |
| hugerth | gl | 0.40426 | 0.74573 | 0.15695 | 91.80328 | 0.35673 | 0.10588 | 1.24855 | 2_Firmicutes(P), Chitinophagales, Eubacteriales, Obscuribacterales, Sphingomonadales |
| labus | gl | 0.20000 | 1.00000 | -0.22222 | 33.33333 | 1.00000 | 0.57234 | 0.71351 | Erysipelotrichales, Lachnospirales |
| liu | gl | 0.71429 | 0.52154 | 0.45132 | 83.60000 | 0.10352 | 0.02504 | 1.91172 | Acetobacterales, Bacillales, Campylobacterales, Exiguobacterales, Flavobacteriales, Frankiales, Microtrichales, Rhizobiales, Sphingomonadales, Synechococcales |
| lopresti | gl | 0.21429 | 0.00000 | -0.12500 | 0.00000 | 0.66667 | 0.55295 | 0.96110 | Acidaminococcales, Veillonellales-Selenomonadales |
| mars | gl | 0.78571 | 0.41143 | 0.48996 | 73.21429 | 0.10606 | 0.04009 | 2.14258 | Flavobacteriales, Fusobacteriales, Micrococcales, Oscillospirales, Staphylococcales |
| nagel | gl | 0.45833 | 0.61989 | 0.18364 | 33.33333 | 0.32727 | 0.13611 | 1.37969 | Burkholderiales, Coriobacteriales, Erysipelotrichales |
| pozuelo | gl | 0.58491 | 0.77758 | 0.12278 | 67.15328 | 0.29462 | 0.07795 | 1.27308 | 1_Proteobacteria(P), 4_Bacteroidota(P), 5_Gammaproteobacteria(C), Corynebacteriales, Micrococcales, Rhizobiales |
| zeber | gl | 0.63889 | 0.53221 | 0.21941 | 50.76923 | 0.25692 | 0.07008 | 1.57276 | Bacteroidales, Clostridia vadinBB60 group, Izemoplasmatales, Saccharimonadales |
| zhu | gl | 0.33333 | 0.60306 | 0.25000 | 50.00000 | 0.38095 | 0.21850 | 1.23969 | Clostridiales, Erysipelotrichales |
| zhuang | gl | 0.50000 | 0.25146 | 0.31500 | 60.00000 | 0.22222 | 0.14169 | 1.66335 | Bacteroidales, Christensenellales |
| agp | mb | 0.90588 | 0.36240 | 0.55533 | 89.50617 | 0.05537 | 0.01762 | 2.40812 | 22_Actinobacteria(C), 7_Bacilli(C), Brevibacillales, Cardiobacteriales, Deinococcales, Frankiales, Micromonosporales, Streptosporangiales, Syntrophomonadales |
| fukui | mb | 0.76316 | 0.31582 | 0.48334 | 73.17073 | 0.10099 | 0.04501 | 2.39279 | Bacteroidales, Campylobacterales, Micrococcales, Saccharimonadales |
| hugerth | mb | 0.74468 | 0.45496 | 0.52401 | 83.33333 | 0.09076 | 0.03760 | 2.21964 | 2_Firmicutes(P), Chitinophagales, Eubacteriales, RF39, Sphingomonadales |
| labus | mb | 0.20000 | 0.00000 | -0.12500 | 50.00000 | 0.66667 | 0.55386 | 0.95715 | Erysipelotrichales, Lachnospirales |
| liu | mb | 0.98980 | 0.18828 | 0.48561 | 79.79094 | 0.06164 | 0.01458 | 2.07846 | 55_Clostridia(C), Campylobacterales, Candidatus Kaiserbacteria, Flavobacteriales, Frankiales, Hydrogenispora, Lachnospirales, Microtrichales, Rhizobiales, Sphingomonadales |
| lopresti | mb | 0.21429 | 0.00000 | -0.12500 | 0.00000 | 0.66667 | 0.55206 | 0.96753 | Acidaminococcales, Veillonellales-Selenomonadales |
| mars | mb | 0.97619 | 0.15389 | 0.52805 | 67.69231 | 0.07927 | 0.03110 | 2.25378 | Corynebacteriales, Micrococcales, Oscillospirales, Propionibacteriales, Staphylococcales |
| nagel | mb | 0.33333 | 0.40625 | 0.21605 | 44.44444 | 0.32143 | 0.18654 | 1.45850 | Bifidobacteriales, Coriobacteriales, Erysipelotrichales |
| pozuelo | mb | 0.69811 | 0.49871 | 0.29931 | 70.88608 | 0.11862 | 0.03829 | 1.73131 | 1_Proteobacteria(P), 5_Gammaproteobacteria(C), 59_Negativicutes(C), Clostridia vadinBB60 group, Rhizobiales, Sphingobacteriales |
| zeber | mb | 0.66667 | 0.42567 | 0.40448 | 57.44681 | 0.17029 | 0.05855 | 1.81342 | Bacteroidales, Clostridia vadinBB60 group, Saccharimonadales, Staphylococcales |
| zhu | mb | 0.52381 | 0.00000 | 0.41500 | 40.00000 | 0.18182 | 0.12514 | 2.05593 | Clostridiales, Rhodobacterales |
| zhuang | mb | 0.55000 | 0.24527 | 0.33471 | 54.54545 | 0.20000 | 0.12691 | 1.68932 | Bacteroidales, Christensenellales |
| agp | slr | 0.07059 | 0.00000 | 0.26000 | 100.00000 | 0.33333 | 0.25474 | 1.64108 | Actinomycetales, Bifidobacteriales, Campylobacterales, Caulobacterales, Clostridia UCG-014, Corynebacteriales, DTU014, Fusobacteriales, Veillonellales-Selenomonadales |
| fukui | slr | 0.07895 | 1.00000 | -0.22222 | 100.00000 | 1.00000 | 0.57615 | 0.70149 | Clostridia UCG-014, Clostridia vadinBB60 group, RF39, Verrucomicrobiales |
| hugerth | slr | 0.14894 | 0.00000 | 0.29167 | 16.66667 | 0.28571 | 0.21214 | 1.61531 | Acholeplasmatales, Methanobacteriales, Staphylococcales, Synergistales, Veillonellales-Selenomonadales |
| labus | slr | 0.13333 | 0.00000 | 0.00000 | 100.00000 | 1.00000 | 0.79756 | 0.70711 | Acidaminococcales, Burkholderiales |
| liu | slr | 0.67347 | 0.55376 | 0.48040 | 89.35185 | 0.10070 | 0.02484 | 2.09274 | Acetobacterales, Bacillales, Campylobacterales, Exiguobacterales, Flavobacteriales, Frankiales, Microtrichales, Rhizobiales, Sphingomonadales, Synechococcales |
| lopresti | slr | 0.14286 | 0.00000 | 0.00000 | 0.00000 | 1.00000 | 0.79496 | 0.71927 | Enterobacterales, Oscillospirales |
| mars | slr | 0.11905 | 0.00000 | 0.21875 | 50.00000 | 0.40000 | 0.31501 | 1.41195 | Bacteroidales, Desulfovibrionales, Erysipelotrichales, RF39, Rhodospirillales |
| nagel | slr | 0.08333 | 0.00000 | 0.00000 | 100.00000 | 1.00000 | 0.79887 | 0.70113 | Christensenellales, RF39 |
| pozuelo | slr | 0.22642 | 0.30000 | 0.46450 | 53.84615 | 0.19697 | 0.11515 | 2.38741 | Bifidobacteriales, Clostridiales, Gastranaerophilales, Peptostreptococcales-Tissierellales, Rhodospirillales, Synergistales |
| zeber | slr | 0.05556 | 0.00000 | 0.00000 | 100.00000 | 1.00000 | 0.79766 | 0.70665 | 3_Bacteroidia(C), Flavobacteriales |
| zhu | slr | 0.09524 | 0.00000 | 0.00000 | 100.00000 | 1.00000 | 0.80741 | 0.66476 | Rhizobiales, Rhodobacterales |
| zhuang | slr | 0.10000 | 0.00000 | 0.00000 | 100.00000 | 1.00000 | 0.79756 | 0.70710 | Burkholderiales, Desulfovibrionales |

### Some Plots

![](../../../outputs/single-network-analysis/plots/Order/global-prop-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/plots/Order/global-prop-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/plots/Order/global-prop-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/plots/Order/global-prop-mb-4.png)<!-- -->

### 3.3.3. SLR

![](../../../outputs/single-network-analysis/plots/Order/global-prop-slr-1.png)<!-- -->

| Taxon                               | HubCount |
|:------------------------------------|---------:|
| Bacteroidales                       |        7 |
| Sphingomonadales                    |        6 |
| Erysipelotrichales                  |        6 |
| Rhizobiales                         |        6 |
| Frankiales                          |        5 |
| Micrococcales                       |        5 |
| Campylobacterales                   |        5 |
| Flavobacteriales                    |        5 |
| Saccharimonadales                   |        4 |
| Veillonellales-Selenomonadales      |        4 |
| Staphylococcales                    |        4 |
| Clostridia vadinBB60 group          |        4 |
| RF39                                |        4 |
| Acetobacterales                     |        3 |
| Lachnospirales                      |        3 |
| Microtrichales                      |        3 |
| Acidaminococcales                   |        3 |
| Oscillospirales                     |        3 |
| Burkholderiales                     |        3 |
| Corynebacteriales                   |        3 |
| Clostridiales                       |        3 |
| Christensenellales                  |        3 |
| Bifidobacteriales                   |        3 |
| 22_Actinobacteria(C)                |        2 |
| Brevibacillales                     |        2 |
| Cardiobacteriales                   |        2 |
| Propionibacteriales                 |        2 |
| Syntrophomonadales                  |        2 |
| 2_Firmicutes(P)                     |        2 |
| Chitinophagales                     |        2 |
| Eubacteriales                       |        2 |
| Bacillales                          |        2 |
| Exiguobacterales                    |        2 |
| Synechococcales                     |        2 |
| Fusobacteriales                     |        2 |
| Coriobacteriales                    |        2 |
| 1_Proteobacteria(P)                 |        2 |
| 5_Gammaproteobacteria(C)            |        2 |
| Rhodobacterales                     |        2 |
| Clostridia UCG-014                  |        2 |
| Synergistales                       |        2 |
| Desulfovibrionales                  |        2 |
| Rhodospirillales                    |        2 |
| 12_Acidimicrobiia(C)                |        1 |
| Bryobacterales                      |        1 |
| Obscuribacterales                   |        1 |
| 4_Bacteroidota(P)                   |        1 |
| Izemoplasmatales                    |        1 |
| 7_Bacilli(C)                        |        1 |
| Deinococcales                       |        1 |
| Micromonosporales                   |        1 |
| Streptosporangiales                 |        1 |
| 55_Clostridia(C)                    |        1 |
| Candidatus Kaiserbacteria           |        1 |
| Hydrogenispora                      |        1 |
| 59_Negativicutes(C)                 |        1 |
| Sphingobacteriales                  |        1 |
| Actinomycetales                     |        1 |
| Caulobacterales                     |        1 |
| DTU014                              |        1 |
| Verrucomicrobiales                  |        1 |
| Acholeplasmatales                   |        1 |
| Methanobacteriales                  |        1 |
| Enterobacterales                    |        1 |
| Gastranaerophilales                 |        1 |
| Peptostreptococcales-Tissierellales |        1 |
| 3_Bacteroidia(C)                    |        1 |
