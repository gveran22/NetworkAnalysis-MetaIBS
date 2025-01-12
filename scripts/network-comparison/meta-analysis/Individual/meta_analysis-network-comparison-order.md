Meta Analysis - Network Comparison (Order)
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

# 2. META-ANAYLSIS

------------------------------------------------------------------------

## GLasso

### Meta-Analysis Plot

![](../../../../outputs/network-comparison/Individual/plots/Order/meta-analysis-glasso-1.png)<!-- -->

### Individual Plots

![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-1.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-2.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-3.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-4.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-5.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-6.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-7.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-8.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-9.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-10.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-glasso-11.png)<!-- -->

### Global Properties

    ##              agp        fukui      hugerth    labus     liu       lopresti 
    ## lccSize1     41         6          28         2         9         2        
    ## lccSize2     43         14         32         1         14        3        
    ## lccSizeRel1  0.640625   0.24       0.7368421  0.2       0.3214286 0.1818182
    ## lccSizeRel2  0.671875   0.56       0.8421053  0.1       0.5       0.2727273
    ## avDiss1      0.6708164  0.7171649  0.6847399  0.7166576 0.705986  0.7581994
    ## avDiss2      0.6697365  0.688924   0.6934112  1         0.7007523 0.7181149
    ## avPath1      1.78276    1.390441   2.021653   0.7166576 1.227052  0.7581994
    ## avPath2      1.697609   1.927786   1.729562   1         2.125009  0.9574865
    ## clustCoef1   0.6229015  0.4202739  0.4834911  0         0.6930643 0        
    ## clustCoef2   0.6758735  0.2106017  0.4258867  0         0.4090039 0        
    ## modularity1  0.288424   0.2083333  0.4826153  0         0.1403061 0        
    ## modularity2  0.3969223  0.4199219  0.4111604  0         0.26125   -0.125   
    ## vertConnect1 1          1          1          1         1         1        
    ## vertConnect2 1          1          1          0         1         1        
    ## edgeConnect1 1          1          1          1         1         1        
    ## edgeConnect2 1          1          1          0         1         1        
    ## natConnect1  0.04810539 0.2574509  0.05347776 0.7955108 0.1699395 0.7873521
    ## natConnect2  0.04783053 0.09828566 0.04542623 0         0.1011277 0.5534462
    ## density1     0.1609756  0.4        0.1666667  1         0.3888889 1        
    ## density2     0.1727575  0.1758242  0.1572581  0         0.2197802 0.6666667
    ## pep1         84.09091   33.33333   77.77778   0         50        0        
    ## pep2         83.33333   62.5       52.5641    0         60        0        
    ##              mars      nagel     pozuelo    zeber      zhuang    
    ## lccSize1     10        3         26         6          1         
    ## lccSize2     16        2         27         15         2         
    ## lccSizeRel1  0.3333333 0.15      0.6046512  0.2307692  0.08333333
    ## lccSizeRel2  0.5333333 0.1       0.627907   0.5769231  0.1666667 
    ## avDiss1      0.6957265 0.7014486 0.694892   0.70936    1         
    ## avDiss2      0.6859261 0.6858671 0.6999954  0.6991636  0.7570607 
    ## avPath1      1.686351  0.9352648 1.421588   1.085017   1         
    ## avPath2      2.30686   0.6858671 1.412003   1.554631   0.7570607 
    ## clustCoef1   0.2913293 0         0.5619863  0.6958953  0         
    ## clustCoef2   0.4472637 0         0.5848027  0.4731653  0         
    ## modularity1  0.3760331 -0.125    0.1851502  0.1171875  0         
    ## modularity2  0.51875   0         0.1453804  0.2939509  0         
    ## vertConnect1 1         1         1          1          0         
    ## vertConnect2 1         1         1          1          1         
    ## edgeConnect1 1         1         1          1          0         
    ## edgeConnect2 1         1         1          1          1         
    ## natConnect1  0.1433535 0.5568351 0.06621613 0.2667359  0         
    ## natConnect2  0.0855368 0.8023467 0.06552479 0.09414401 0.7875592 
    ## density1     0.2444444 0.6666667 0.2553846  0.5333333  0         
    ## density2     0.1666667 1         0.2621083  0.2190476  1         
    ## pep1         63.63636  50        59.03614   50         0         
    ## pep2         85        100       55.43478   56.52174   0

## MB

### Meta-Analysis Plot

![](../../../../outputs/network-comparison/Individual/plots/Order/meta-analysis-mb-1.png)<!-- -->

### Individual Plots

![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-1.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-2.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-3.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-4.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-5.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-6.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-7.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-8.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-9.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-10.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-mb-11.png)<!-- -->

### Global Properties

    ##              agp        fukui      hugerth    labus     liu        lopresti 
    ## lccSize1     52         4          17         2         9          2        
    ## lccSize2     51         17         34         2         14         3        
    ## lccSizeRel1  0.8125     0.16       0.4473684  0.2       0.3214286  0.1818182
    ## lccSizeRel2  0.796875   0.68       0.8947368  0.2       0.5        0.2727273
    ## avDiss1      0.6848351  0.7160374  0.6848059  0.7223595 0.7043328  0.7779132
    ## avDiss2      0.6809628  0.6906705  0.6858086  0.7215003 0.6908534  0.7223955
    ## avPath1      2.062771   1.193986   2.091708   0.7223595 1.398375   0.7779132
    ## avPath2      2.234596   2.590235   2.046216   0.7215003 2.351755   0.963194 
    ## clustCoef1   0.438984   0          0.2071235  0         0.4851666  0        
    ## clustCoef2   0.4343792  0.1014402  0.3390864  0         0.2568916  0        
    ## modularity1  0.4577897  0.1666667  0.4307851  0         0.3140496  0        
    ## modularity2  0.5273993  0.5519031  0.4793163  0         0.4101562  -0.125   
    ## vertConnect1 1          1          1          1         1          1        
    ## vertConnect2 1          1          1          1         1          1        
    ## edgeConnect1 1          1          1          1         1          1        
    ## edgeConnect2 1          1          1          1         1          1        
    ## natConnect1  0.02643817 0.4061114  0.08011932 0.7943182 0.1633434  0.7839147
    ## natConnect2  0.0267048  0.07761042 0.03909001 0.7944964 0.09787065 0.5526567
    ## density1     0.08220211 0.5        0.1617647  1         0.3055556  1        
    ## density2     0.07921569 0.125      0.1051693  1         0.1758242  0.6666667
    ## pep1         84.40367   33.33333   63.63636   0         54.54545   0        
    ## pep2         90.09901   58.82353   64.40678   0         62.5       0        
    ##              mars       nagel     pozuelo    zeber      zhuang    
    ## lccSize1     12         6         30         7          1         
    ## lccSize2     21         2         33         17         2         
    ## lccSizeRel1  0.4        0.3       0.6976744  0.2692308  0.08333333
    ## lccSizeRel2  0.7        0.1       0.7674419  0.6538462  0.1666667 
    ## avDiss1      0.6892061  0.7000745 0.6954247  0.7111036  1         
    ## avDiss2      0.6882127  0.6716302 0.6968291  0.6937194  0.7942105 
    ## avPath1      1.763979   1.312052  1.778081   1.428522   1         
    ## avPath2      2.719111   0.6716302 1.878488   1.839303   0.7942105 
    ## clustCoef1   0.2112182  0         0.13991    0          0         
    ## clustCoef2   0.3103452  0         0.260367   0.2736541  0         
    ## modularity1  0.4289941  0.22      0.3126594  0.3194444  0         
    ## modularity2  0.5789931  0         0.3331299  0.4445983  0         
    ## vertConnect1 1          1         1          1          0         
    ## vertConnect2 1          1         1          1          1         
    ## edgeConnect1 1          1         1          1          0         
    ## edgeConnect2 1          1         1          1          1         
    ## natConnect1  0.1163039  0.25539   0.04492032 0.2120084  0         
    ## natConnect2  0.06208362 0.8057307 0.04094101 0.07832051 0.7812868 
    ## density1     0.1969697  0.3333333 0.1287356  0.2857143  0         
    ## density2     0.1142857  1         0.1212121  0.1397059  1         
    ## pep1         61.53846   40        66.07143   50         0         
    ## pep2         70.83333   100       56.25      57.89474   0

## SLR

### Meta-Analysis Plot

![](../../../../outputs/network-comparison/Individual/plots/Order/meta-analysis-slr-1.png)<!-- -->

### Individual Plots

![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-1.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-2.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-3.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-4.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-5.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-6.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-7.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-8.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-9.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-10.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-11.png)<!-- -->![](../../../../outputs/network-comparison/Individual/plots/Order/single-network-slr-12.png)<!-- -->

### Global Properties

    ##              agp        fukui     hugerth    labus     liu        lopresti 
    ## lccSize1     32         2         21         4         2          3        
    ## lccSize2     36         12        22         2         18         3        
    ## lccSizeRel1  0.5        0.08      0.5526316  0.4       0.07142857 0.2727273
    ## lccSizeRel2  0.5625     0.48      0.5789474  0.2       0.6428571  0.2727273
    ## avDiss1      0.7024364  0.776751  0.6986779  0.7151544 0.8140755  0.7286807
    ## avDiss2      0.6977606  0.7051905 0.6990338  0.7011821 0.7089299  0.7290229
    ## avPath1      1.626841   0.776751  3.057362   1.192643  0.8140755  0.9715743
    ## avPath2      2.153518   1.737985  2.921447   0.7011821 2.029441   0.9720305
    ## clustCoef1   0.3572371  0         0          0         0          0        
    ## clustCoef2   0.4413905  0.1663617 0          0         0.4243104  0        
    ## modularity1  0.2817519  0         0.59375    0.1666667 0          -0.125   
    ## modularity2  0.5317397  0.3577778 0.5295139  0         0.4356509  -0.125   
    ## vertConnect1 1          1         1          1         1          1        
    ## vertConnect2 1          1         1          1         1          1        
    ## edgeConnect1 1          1         1          1         1          1        
    ## edgeConnect2 1          1         1          1         1          1        
    ## natConnect1  0.04491827 0.7841095 0.06049053 0.4055486 0.7783472  0.5515419
    ## natConnect2  0.03717547 0.1167278 0.05807865 0.7988635 0.07483593 0.5513703
    ## density1     0.1653226  1         0.0952381  0.5       1          0.6666667
    ## density2     0.1126984  0.2272727 0.1038961  1         0.1699346  0.6666667
    ## pep1         67.07317   0         70         0         0          0        
    ## pep2         64.78873   60        62.5       100       42.30769   0        
    ##              mars       nagel     pozuelo    zeber     zhuang   
    ## lccSize1     5          3         23         3         2        
    ## lccSize2     2          2         19         4         2        
    ## lccSizeRel1  0.1666667  0.15      0.5348837  0.1153846 0.1666667
    ## lccSizeRel2  0.06666667 0.1       0.4418605  0.1538462 0.1666667
    ## avDiss1      0.6972057  0.6965062 0.7003339  0.7055987 0.6989653
    ## avDiss2      0.7340062  0.6950517 0.7073827  0.6812237 0.7686543
    ## avPath1      1.191387   0.9286749 2.228933   0.9407983 0.6989653
    ## avPath2      0.7340062  0.6950517 1.683193   1.139669  0.7686543
    ## clustCoef1   0.5858378  0         0.07896239 0         0        
    ## clustCoef2   0          0         0.3790025  0         0        
    ## modularity1  0.22       -0.125    0.4472222  -0.125    0        
    ## modularity2  0          0         0.3343426  0.1666667 0        
    ## vertConnect1 1          1         1          1         1        
    ## vertConnect2 1          1         1          1         1        
    ## edgeConnect1 1          1         1          1         1        
    ## edgeConnect2 1          1         1          1         1        
    ## natConnect1  0.3231086  0.5577613 0.05628357 0.5558926 0.7993576
    ## natConnect2  0.7919537  0.8002381 0.07258143 0.4119523 0.785494 
    ## density1     0.5        0.6666667 0.1185771  0.6666667 1        
    ## density2     1          1         0.1988304  0.5       1        
    ## pep1         80         100       63.33333   100       100      
    ## pep2         0          100       58.82353   66.66667  0
