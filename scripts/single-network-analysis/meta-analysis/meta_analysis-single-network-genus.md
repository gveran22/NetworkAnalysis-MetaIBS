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

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 65 | 94 | 71 | 65 | 92 | 82 | 94 | 94 | 92 | 93 | 82 | 81 |
| lccSizeRel1 | 0.65 | 0.94 | 0.71 | 0.65 | 0.92 | 0.82 | 0.94 | 0.94 | 0.92 | 0.93 | 0.82 | 0.81 |
| avDiss1 | 0.683520859608678 | 0.695508958207824 | 0.68206444248382 | 0.653689375761609 | 0.680937573812101 | 0.655994715539515 | 0.680619325838493 | 0.694047027191546 | 0.695354824908217 | 0.691262302973964 | 0.696997659855179 | 0.689905330676439 |
| avPath1 | 1.8658688992353 | 2.09263148230427 | 1.76275114569351 | 3.86100631637969 | 1.87126351548386 | 3.58853085972095 | 2.10939632435737 | 2.53334307884992 | 2.08924821214753 | 2.03216154630871 | 2.36088208704014 | 2.70845114140831 |
| clustCoef1 | 0.648897988880081 | 0.29814285210272 | 0.587505946085482 | 0.421335793421067 | 0.428368370726821 | 0.435424033998171 | 0.369245361888669 | 0.277283714682487 | 0.438121591626599 | 0.273953745607432 | 0.487456372852235 | 0.343017845027415 |
| modularity1 | 0.418479240767756 | 0.412568108270179 | 0.394516983199416 | 0.72023757533409 | 0.425200984307764 | 0.70007812030527 | 0.397761458896701 | 0.578404954855316 | 0.370737240075614 | 0.382570239334027 | 0.404506332828011 | 0.595423201730665 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0314477467781196 | 0.0160800533169269 | 0.0275646602563352 | 0.0204023886400028 | 0.0233495264829871 | 0.0156543438491922 | 0.0233376875578032 | 0.0143585743085393 | 0.0190054071469611 | 0.0174614179500067 | 0.0230231642160516 | 0.0165038458780681 |
| density1 | 0.11875 | 0.071150766415008 | 0.119114688128773 | 0.0514423076923077 | 0.0941232680363115 | 0.038843721770551 | 0.0787005261953786 | 0.051018073667353 | 0.0824175824175824 | 0.072463768115942 | 0.0861186389641674 | 0.0530864197530864 |
| pep1 | 80.1619433198381 | 64.6302250803859 | 67.9054054054054 | 99.0654205607477 | 68.5279187817259 | 96.8992248062015 | 66.5697674418605 | 63.2286995515695 | 62.8985507246377 | 60.3225806451613 | 55.5944055944056 | 66.2790697674419 |

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-12.png)<!-- -->
\### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 90 | 97 | 85 | 98 | 100 | 99 | 100 | 100 | 96 | 100 | 100 | 98 |
| lccSizeRel1 | 0.9 | 0.97 | 0.85 | 0.98 | 1 | 0.99 | 1 | 1 | 0.96 | 1 | 1 | 0.98 |
| avDiss1 | 0.688294660714036 | 0.69636233618507 | 0.690152918004346 | 0.661163262393192 | 0.691543800795555 | 0.654454075308491 | 0.692858561995277 | 0.691915770242329 | 0.694934836718616 | 0.694747866983756 | 0.685098877174403 | 0.689816462545654 |
| avPath1 | 2.35202431796687 | 2.10649502553295 | 2.25774573929018 | 3.0742772402051 | 1.92555515372572 | 2.76510811992778 | 2.08350803504889 | 2.2710416612422 | 2.20364920276538 | 1.99758346912443 | 2.41886487492929 | 2.67186740895412 |
| clustCoef1 | 0.466776030375109 | 0.168275962389785 | 0.410906397828782 | 0.0981628500374906 | 0.218353934990481 | 0.0634818833688867 | 0.166866109163709 | 0.0756325253781296 | 0.230469707060141 | 0.126432638965612 | 0.142755841491832 | 0.0786477955902612 |
| modularity1 | 0.554946899414062 | 0.450013734476427 | 0.519079411073694 | 0.622175933344015 | 0.426344915512606 | 0.566524897436863 | 0.46209991349481 | 0.471775762167032 | 0.435517355371901 | 0.386463447592451 | 0.485412177494708 | 0.543065900770105 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0162312507904918 | 0.0137813106697014 | 0.0166694129563964 | 0.0127517766303322 | 0.0144291685230007 | 0.0127666486275706 | 0.0134500988688268 | 0.0127340717470919 | 0.0145859467829369 | 0.0137237553244396 | 0.0130020306175739 | 0.0126653505261135 |
| density1 | 0.0639200998751561 | 0.0564862542955326 | 0.0633053221288515 | 0.0332421628445193 | 0.0664646464646465 | 0.033601319315605 | 0.054949494949495 | 0.0446464646464646 | 0.0603070175438596 | 0.0602020202020202 | 0.0462626262626263 | 0.0374500315590154 |
| pep1 | 86.328125 | 65.7794676806084 | 76.5486725663717 | 88.6075949367088 | 68.3890577507599 | 89.5705521472393 | 67.2794117647059 | 57.4660633484163 | 67.6363636363636 | 60.4026845637584 | 61.5720524017467 | 57.8651685393258 |

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-slr-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-12.png)<!-- -->

### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 11 | 79 | 10 | 72 | 84 | 99 | 81 | 87 | 82 | 87 | 77 | 76 |
| lccSizeRel1 | 0.11 | 0.79 | 0.1 | 0.72 | 0.84 | 0.99 | 0.81 | 0.87 | 0.82 | 0.87 | 0.77 | 0.76 |
| avDiss1 | 0.681585057688382 | 0.706598124248728 | 0.713545511206251 | 0.705367982292231 | 0.703079050746049 | 0.683211580982192 | 0.706970930525911 | 0.701963535982266 | 0.701165722086699 | 0.698370312781416 | 0.703074652487027 | 0.708165866846396 |
| avPath1 | 2.38667414521821 | 2.3707992880601 | 1.91672261373546 | 1.95806385801618 | 1.91385993558944 | 2.08452899947848 | 2.28580057069957 | 2.14899862014077 | 1.84978013081304 | 2.31049806197597 | 1.87343550910137 | 1.93791096960998 |
| clustCoef1 | 0 | 0.10745889091561 | 0 | 0.195775250073537 | 0.160581078725564 | 0.236430623477709 | 0.182492488217619 | 0.290359975848526 | 0.368257598877322 | 0.137459923754634 | 0.356391105158531 | 0.173629417773071 |
| modularity1 | 0.435 | 0.478086572676906 | 0.355 | 0.313435433344936 | 0.367159821253716 | 0.451625728288684 | 0.467719037317469 | 0.473499964993349 | 0.430138259007307 | 0.467183673469388 | 0.384259560353798 | 0.371058514135437 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.126947241003703 | 0.0156360464020043 | 0.140478356159199 | 0.0195914810689348 | 0.016009591133165 | 0.0136847800079947 | 0.0155614897988292 | 0.0150457177081146 | 0.0174845837367627 | 0.0144542106887611 | 0.0188745137450467 | 0.017375927839805 |
| density1 | 0.181818181818182 | 0.0483609217786433 | 0.222222222222222 | 0.0864632237871674 | 0.0679862306368331 | 0.0542156256441971 | 0.0530864197530864 | 0.0553327987169206 | 0.0758807588075881 | 0.0467789361133387 | 0.0847573479152427 | 0.068421052631579 |
| pep1 | 90 | 51.006711409396 | 60 | 49.7737556561086 | 48.5232067510549 | 78.3269961977186 | 43.6046511627907 | 55.0724637681159 | 61.1111111111111 | 50.8571428571429 | 52.0161290322581 | 47.6923076923077 |
