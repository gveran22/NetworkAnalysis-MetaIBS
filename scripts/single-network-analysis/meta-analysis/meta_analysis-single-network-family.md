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
library(kableExtra)
```

------------------------------------------------------------------------

# 2. META-ANALYSIS

------------------------------------------------------------------------

## GLasso

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-glasso-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 85 | 54 | 60 | 9 | 143 | 10 | 59 | 17 | 67 | 54 | 15 | 17 |
| lccSizeRel1 | 0.534591194968553 | 0.818181818181818 | 0.681818181818182 | 0.375 | 0.841176470588235 | 0.434782608695652 | 0.808219178082192 | 0.395348837209302 | 0.656862745098039 | 0.805970149253731 | 0.428571428571429 | 0.515151515151515 |
| avDiss1 | 0.662429034251132 | 0.692494591560237 | 0.671895561710602 | 0.708198236795545 | 0.66749040755003 | 0.706374761281987 | 0.691722641843888 | 0.702713935561142 | 0.687986060494955 | 0.693862712183184 | 0.706369318453593 | 0.703023834208629 |
| avPath1 | 1.72017043504308 | 2.10685888899396 | 1.79831068343951 | 1.79978550906889 | 1.77728535009589 | 1.83432832488777 | 1.96717226419015 | 1.82654257666179 | 1.64028274881498 | 2.10141091816064 | 1.58844025688646 | 1.93466146439973 |
| clustCoef1 | 0.721511150549166 | 0.471582303592243 | 0.651338153788387 | 0 | 0.445570878918826 | 0 | 0.472569736935148 | 0.454733917629925 | 0.668951397786687 | 0.316016711489906 | 0.606220582306464 | 0.425670005921011 |
| modularity1 | 0.235157644259369 | 0.398718164951931 | 0.266689754382226 | 0.3828125 | 0.418598180529301 | 0.36 | 0.408734567901234 | 0.326086956521739 | 0.169408786166651 | 0.359954833984375 | 0.305473372781065 | 0.396825396825397 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0716654270381417 | 0.029377217826254 | 0.0496476099641446 | 0.158065212257023 | 0.0212006734961743 | 0.141001578843248 | 0.0268289742084101 | 0.0797735132520083 | 0.0526755092440833 | 0.0264297785493569 | 0.0957167305666659 | 0.0787656386080988 |
| density1 | 0.172268907563025 | 0.107617051013277 | 0.154802259887006 | 0.222222222222222 | 0.0906136117403723 | 0.222222222222222 | 0.105201636469901 | 0.169117647058824 | 0.165988240615106 | 0.0894479385045423 | 0.247619047619048 | 0.154411764705882 |
| pep1 | 96.4227642276423 | 68.8311688311688 | 87.2262773722628 | 37.5 | 85.5434782608696 | 20 | 70.5555555555556 | 52.1739130434783 | 69.4822888283379 | 61.71875 | 42.3076923076923 | 52.3809523809524 |

Global network properties

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 153 | 63 | 75 | 5 | 170 | 10 | 72 | 29 | 94 | 66 | 18 | 20 |
| lccSizeRel1 | 0.962264150943396 | 0.954545454545455 | 0.852272727272727 | 0.208333333333333 | 1 | 0.434782608695652 | 0.986301369863014 | 0.674418604651163 | 0.92156862745098 | 0.985074626865672 | 0.514285714285714 | 0.606060606060606 |
| avDiss1 | 0.685781032251058 | 0.694524359267811 | 0.683452444046151 | 0.702259014049142 | 0.681764691635611 | 0.703738883525277 | 0.690910041666369 | 0.702997593321346 | 0.69567647028295 | 0.687615475915117 | 0.708023804826509 | 0.702334552747312 |
| avPath1 | 2.22243680539441 | 2.17494109843637 | 1.9693967910723 | 1.3973626207529 | 1.80289346472507 | 1.82596401323556 | 2.12075597015007 | 2.57993792694107 | 2.16927939822225 | 2.26632367052457 | 2.15218350313103 | 2.14937441673878 |
| clustCoef1 | 0.349161712560609 | 0.282339813868713 | 0.436946627149203 | 0 | 0.169513455894334 | 0 | 0.199057687147153 | 0.219359388367274 | 0.389733210997839 | 0.121333784523161 | 0.187475118274439 | 0.260341819040766 |
| modularity1 | 0.486325190793822 | 0.442980096826251 | 0.478061322568243 | 0.21875 | 0.384439239800746 | 0.36 | 0.466597808199993 | 0.504747991234478 | 0.375432894666738 | 0.442016 | 0.48125 | 0.504725897920605 |
| vertConnect1 | 1 | 1 | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.010935580494926 | 0.0212881265379349 | 0.0196760017407197 | 0.315833362267956 | 0.00942974321133175 | 0.141248198904508 | 0.0184937162192857 | 0.0436279230921132 | 0.0163099388054905 | 0.0194456605371336 | 0.072532704394851 | 0.0649822746489298 |
| density1 | 0.04437564499484 | 0.0732206861239119 | 0.0735135135135135 | 0.4 | 0.0549251653324052 | 0.222222222222222 | 0.0661189358372457 | 0.0911330049261084 | 0.0626858842370167 | 0.0582750582750583 | 0.130718954248366 | 0.121052631578947 |
| pep1 | 89.5348837209302 | 66.4335664335664 | 87.2549019607843 | 50 | 75.9188846641318 | 20 | 68.0473372781065 | 45.945945945946 | 72.992700729927 | 64.8 | 40 | 47.8260869565217 |

Global network properties

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 91 | 36 | 50 | 2 | 121 | 3 | 49 | 23 | 63 | 38 | 3 | 3 |
| lccSizeRel1 | 0.572327044025157 | 0.545454545454545 | 0.568181818181818 | 0.0833333333333333 | 0.711764705882353 | 0.130434782608696 | 0.671232876712329 | 0.534883720930233 | 0.617647058823529 | 0.567164179104478 | 0.0857142857142857 | 0.0909090909090909 |
| avDiss1 | 0.700646306387569 | 0.706223468133027 | 0.700787947783672 | 0.684185110732245 | 0.706874633734318 | 0.705158983118493 | 0.703804720633022 | 0.706979742338351 | 0.707395046022618 | 0.707100053135472 | 0.718091928245468 | 0.713509996719303 |
| avPath1 | 1.66882155578144 | 1.46321663356725 | 1.53410961723649 | 0.684185110732245 | 1.76460145669788 | 0.940211977491324 | 1.90550504260161 | 2.18467959128068 | 1.80058321033371 | 1.74047082408843 | 0.95745590432729 | 0.951346662292404 |
| clustCoef1 | 0.450602526224225 | 0.366220882100833 | 0.411254850156885 | 0 | 0.48331494401136 | 0 | 0.281664361416238 | 0.303749015841667 | 0.264584792542678 | 0.341606685971297 | 0 | 0 |
| modularity1 | 0.263809387797639 | 0.2634375 | 0.333574905354246 | 0 | 0.156336343638198 | -0.125 | 0.399059178397012 | 0.50244140625 | 0.338858131487889 | 0.409197530864198 | -0.125 | -0.125 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0247445975422372 | 0.0421537841815226 | 0.0339835736952382 | 0.802739159903281 | 0.0420341149229 | 0.555983770134316 | 0.0274613545527469 | 0.0567733709712566 | 0.0216389617695608 | 0.0356920915212157 | 0.553581387047936 | 0.554481096319041 |
| density1 | 0.124297924297924 | 0.19047619047619 | 0.175510204081633 | 1 | 0.135123966942149 | 0.666666666666667 | 0.100340136054422 | 0.126482213438735 | 0.0870455709165387 | 0.128022759601707 | 0.666666666666667 | 0.666666666666667 |
| pep1 | 45.5795677799607 | 50.8333333333333 | 43.7209302325581 | 100 | 41.5902140672783 | 50 | 50 | 56.25 | 33.5294117647059 | 55.5555555555556 | 50 | 50 |

Global network properties
