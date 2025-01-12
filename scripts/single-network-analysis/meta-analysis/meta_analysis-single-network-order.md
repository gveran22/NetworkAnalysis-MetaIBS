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

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 33 | 23 | 19 | 3 | 70 | 3 | 33 | 11 | 31 | 23 | 7 | 10 |
| lccSizeRel1 | 0.388235294117647 | 0.605263157894737 | 0.404255319148936 | 0.2 | 0.714285714285714 | 0.214285714285714 | 0.785714285714286 | 0.458333333333333 | 0.584905660377358 | 0.638888888888889 | 0.333333333333333 | 0.5 |
| avDiss1 | 0.658886208568376 | 0.676858047943762 | 0.664784668281467 | 0.713513024109345 | 0.683990727778986 | 0.720822968214867 | 0.693582918940938 | 0.704671485080845 | 0.685352067500124 | 0.700437105632864 | 0.707255348497786 | 0.701780423848017 |
| avPath1 | 1.11919501145032 | 1.70828760059253 | 1.24854838207781 | 0.713513024109345 | 1.91172404921003 | 0.961097290953156 | 2.14257838137211 | 1.37968854100236 | 1.27307821936622 | 1.57275830134327 | 1.23968935666481 | 1.66335022772168 |
| clustCoef1 | 0.775109454886084 | 0.578489647664439 | 0.745733975935263 | 1 | 0.521536607440231 | 0 | 0.411427850093052 | 0.619888138227922 | 0.777583143736999 | 0.532212561154 | 0.603063612580471 | 0.25146015380179 |
| modularity1 | 0.151812032751869 | 0.34664461247637 | 0.156947057242677 | -0.222222222222222 | 0.45132 | -0.125 | 0.489955357142857 | 0.183641975308642 | 0.122782247322713 | 0.219408284023669 | 0.25 | 0.315 |
| vertConnect1 | 1 | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 2 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.142793903096842 | 0.0649973962619939 | 0.10588019212637 | 0.572341307959345 | 0.0250397400511789 | 0.552946950827164 | 0.0400934599249989 | 0.136113095874708 | 0.0779452444112079 | 0.0700788675023799 | 0.218495961698751 | 0.141685941384823 |
| density1 | 0.501893939393939 | 0.181818181818182 | 0.35672514619883 | 1 | 0.10351966873706 | 0.666666666666667 | 0.106060606060606 | 0.327272727272727 | 0.294623655913978 | 0.256916996047431 | 0.380952380952381 | 0.222222222222222 |
| pep1 | 98.4905660377358 | 73.9130434782609 | 91.8032786885246 | 33.3333333333333 | 83.6 | 0 | 73.2142857142857 | 33.3333333333333 | 67.1532846715328 | 50.7692307692308 | 50 | 60 |

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-12.png)<!-- -->

### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 77 | 29 | 35 | 3 | 97 | 3 | 41 | 8 | 37 | 24 | 11 | 11 |
| lccSizeRel1 | 0.905882352941176 | 0.763157894736842 | 0.74468085106383 | 0.2 | 0.989795918367347 | 0.214285714285714 | 0.976190476190476 | 0.333333333333333 | 0.69811320754717 | 0.666666666666667 | 0.523809523809524 | 0.55 |
| avDiss1 | 0.6783200704136 | 0.687727230344644 | 0.678044899569827 | 0.717861706925037 | 0.676573965602452 | 0.725649043816422 | 0.691829731403922 | 0.70499700376799 | 0.692574436373632 | 0.696267516283932 | 0.709754488797441 | 0.701012936353421 |
| avPath1 | 2.40811906614453 | 2.39279249909985 | 2.21963502383403 | 0.957148942566716 | 2.0784593590407 | 0.967532058421895 | 2.25378016054739 | 1.45850126229308 | 1.73130625127215 | 1.81342085106891 | 2.05593239315743 | 1.68931952035005 |
| clustCoef1 | 0.362401601437677 | 0.315822929650036 | 0.454961715995446 | 0 | 0.188275091150283 | 0 | 0.153893100495699 | 0.406254526928855 | 0.498707341033337 | 0.425672489066256 | 0 | 0.245273867495073 |
| modularity1 | 0.555326931870142 | 0.483343248066627 | 0.52400548696845 | -0.125 | 0.485613519649383 | -0.125 | 0.528047337278107 | 0.216049382716049 | 0.299311007851306 | 0.404481665912177 | 0.415 | 0.334710743801653 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0176171813943167 | 0.0450090952856019 | 0.0376038240799855 | 0.553863801986142 | 0.0145821420288644 | 0.552062621771074 | 0.0310996180644101 | 0.186540559305829 | 0.0382904991715624 | 0.0585517157413078 | 0.125141463954739 | 0.126908265234667 |
| density1 | 0.0553656869446343 | 0.100985221674877 | 0.0907563025210084 | 0.666666666666667 | 0.0616408934707904 | 0.666666666666667 | 0.0792682926829268 | 0.321428571428571 | 0.118618618618619 | 0.170289855072464 | 0.181818181818182 | 0.2 |
| pep1 | 89.5061728395062 | 73.1707317073171 | 83.3333333333333 | 50 | 79.7909407665505 | 0 | 67.6923076923077 | 44.4444444444444 | 70.8860759493671 | 57.4468085106383 | 40 | 54.5454545454545 |

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-12.png)<!-- -->

### Global Properties

|  | agp | fukui | hugerth | labus | liu | lopresti | mars | nagel | pozuelo | zeber | zhu | zhuang |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| lccSize1 | 45 | 19 | 24 | 2 | 64 | 3 | 16 | 3 | 8 | 21 | 3 | 2 |
| lccSizeRel1 | 0.529411764705882 | 0.5 | 0.51063829787234 | 0.133333333333333 | 0.653061224489796 | 0.214285714285714 | 0.380952380952381 | 0.125 | 0.150943396226415 | 0.583333333333333 | 0.142857142857143 | 0.1 |
| avDiss1 | 0.701642779425181 | 0.706235098631646 | 0.698402122221439 | 0.705289332852998 | 0.707374270525729 | 0.722823126142369 | 0.702613744378417 | 0.710821129806322 | 0.718822725268692 | 0.707623809917162 | 0.707245072618375 | 0.743202523845768 |
| avPath1 | 1.45537112716362 | 1.55168264760283 | 1.34866848074102 | 0.705289332852998 | 1.77980630577923 | 0.963764168189826 | 2.20720130210823 | 0.94776150640843 | 1.75085070683059 | 1.86485613682646 | 0.942993430157833 | 0.743202523845768 |
| clustCoef1 | 0.483078445070843 | 0.227288203175652 | 0.460495857998052 | 0 | 0.399105235440164 | 0 | 0 | 0 | 0 | 0.160787918439202 | 0 | 0 |
| modularity1 | 0.291517431523616 | 0.266975308641975 | 0.263948961937716 | 0 | 0.202249662236466 | -0.125 | 0.498456790123457 | -0.125 | 0.336734693877551 | 0.40495867768595 | -0.125 | 0 |
| vertConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| edgeConnect1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| natConnect1 | 0.0400381409590422 | 0.0736873961849749 | 0.0633748251529095 | 0.79795729064521 | 0.0334233118963411 | 0.552534048658085 | 0.0830344697417884 | 0.555117045373099 | 0.180563227212174 | 0.0633791690637251 | 0.555568416925851 | 0.790155027253564 |
| density1 | 0.198989898989899 | 0.210526315789474 | 0.246376811594203 | 1 | 0.15327380952381 | 0.666666666666667 | 0.15 | 0.666666666666667 | 0.25 | 0.157142857142857 | 0.666666666666667 | 1 |
| pep1 | 47.7157360406091 | 58.3333333333333 | 67.6470588235294 | 100 | 45.3074433656958 | 0 | 77.7777777777778 | 50 | 42.8571428571429 | 54.5454545454545 | 50 | 0 |
