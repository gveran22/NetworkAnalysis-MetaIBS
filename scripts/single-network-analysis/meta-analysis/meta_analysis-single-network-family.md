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

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-glasso-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-glasso-12.png)<!-- -->
\### Global Properties

    ##              agp        fukui      hugerth    labus     liu        lopresti 
    ## lccSize1     85         54         60         9         143        10       
    ## lccSizeRel1  0.5345912  0.8181818  0.6818182  0.375     0.8411765  0.4347826
    ## avDiss1      0.662429   0.6924946  0.6718956  0.7081982 0.6674904  0.7063748
    ## avPath1      1.72017    2.106859   1.798311   1.799786  1.777285   1.834328 
    ## clustCoef1   0.7215112  0.4715823  0.6513382  0         0.4455709  0        
    ## modularity1  0.2351576  0.3987182  0.2666898  0.3828125 0.4185982  0.36     
    ## vertConnect1 1          1          1          1         1          1        
    ## edgeConnect1 1          1          1          1         1          1        
    ## natConnect1  0.07166543 0.02937722 0.04964761 0.1580652 0.02120067 0.1410016
    ## density1     0.1722689  0.1076171  0.1548023  0.2222222 0.09061361 0.2222222
    ## pep1         96.42276   68.83117   87.22628   37.5      85.54348   20       
    ##              mars       nagel      pozuelo    zeber      zhu        zhuang    
    ## lccSize1     59         17         67         54         15         17        
    ## lccSizeRel1  0.8082192  0.3953488  0.6568627  0.8059701  0.4285714  0.5151515 
    ## avDiss1      0.6917226  0.7027139  0.6879861  0.6938627  0.7063693  0.7030238 
    ## avPath1      1.967172   1.826543   1.640283   2.101411   1.58844    1.934661  
    ## clustCoef1   0.4725697  0.4547339  0.6689514  0.3160167  0.6062206  0.42567   
    ## modularity1  0.4087346  0.326087   0.1694088  0.3599548  0.3054734  0.3968254 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.02682897 0.07977351 0.05267551 0.02642978 0.09571673 0.07876564
    ## density1     0.1052016  0.1691176  0.1659882  0.08944794 0.247619   0.1544118 
    ## pep1         70.55556   52.17391   69.48229   61.71875   42.30769   52.38095

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-mb-12.png)<!-- -->
\### Global Properties

    ##              agp        fukui      hugerth    labus     liu         lopresti 
    ## lccSize1     153        63         75         5         170         10       
    ## lccSizeRel1  0.9622642  0.9545455  0.8522727  0.2083333 1           0.4347826
    ## avDiss1      0.685781   0.6945244  0.6834524  0.702259  0.6817647   0.7037389
    ## avPath1      2.222437   2.174941   1.969397   1.397363  1.802893    1.825964 
    ## clustCoef1   0.3491617  0.2823398  0.4369466  0         0.1695135   0        
    ## modularity1  0.4863252  0.4429801  0.4780613  0.21875   0.3844392   0.36     
    ## vertConnect1 1          1          1          1         2           1        
    ## edgeConnect1 1          1          1          1         2           1        
    ## natConnect1  0.01093558 0.02128813 0.019676   0.3158334 0.009429743 0.1412482
    ## density1     0.04437564 0.07322069 0.07351351 0.4       0.05492517  0.2222222
    ## pep1         89.53488   66.43357   87.2549    50        75.91888    20       
    ##              mars       nagel      pozuelo    zeber      zhu       zhuang    
    ## lccSize1     72         29         94         66         18        20        
    ## lccSizeRel1  0.9863014  0.6744186  0.9215686  0.9850746  0.5142857 0.6060606 
    ## avDiss1      0.69091    0.7029976  0.6956765  0.6876155  0.7080238 0.7023346 
    ## avPath1      2.120756   2.579938   2.169279   2.266324   2.152184  2.149374  
    ## clustCoef1   0.1990577  0.2193594  0.3897332  0.1213338  0.1874751 0.2603418 
    ## modularity1  0.4665978  0.504748   0.3754329  0.442016   0.48125   0.5047259 
    ## vertConnect1 1          1          1          1          1         1         
    ## edgeConnect1 1          1          1          1          1         1         
    ## natConnect1  0.01849372 0.04362792 0.01630994 0.01944566 0.0725327 0.06498227
    ## density1     0.06611894 0.091133   0.06268588 0.05827506 0.130719  0.1210526 
    ## pep1         68.04734   45.94595   72.9927    64.8       40        47.82609

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Family/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Family/single-network-slr-12.png)<!-- -->
\### Global Properties

    ##              agp       fukui      hugerth    labus      liu        lopresti 
    ## lccSize1     91        36         50         2          121        3        
    ## lccSizeRel1  0.572327  0.5454545  0.5681818  0.08333333 0.7117647  0.1304348
    ## avDiss1      0.7006463 0.7062235  0.7007879  0.6841851  0.7068746  0.705159 
    ## avPath1      1.668822  1.463217   1.53411    0.6841851  1.764601   0.940212 
    ## clustCoef1   0.4506025 0.3662209  0.4112549  0          0.4833149  0        
    ## modularity1  0.2638094 0.2634375  0.3335749  0          0.1563363  -0.125   
    ## vertConnect1 1         1          1          1          1          1        
    ## edgeConnect1 1         1          1          1          1          1        
    ## natConnect1  0.0247446 0.04215378 0.03398357 0.8027392  0.04203411 0.5559838
    ## density1     0.1242979 0.1904762  0.1755102  1          0.135124   0.6666667
    ## pep1         45.57957  50.83333   43.72093   100        41.59021   50       
    ##              mars       nagel      pozuelo    zeber      zhu        zhuang    
    ## lccSize1     49         23         63         38         3          3         
    ## lccSizeRel1  0.6712329  0.5348837  0.6176471  0.5671642  0.08571429 0.09090909
    ## avDiss1      0.7038047  0.7069797  0.707395   0.7071001  0.7180919  0.71351   
    ## avPath1      1.905505   2.18468    1.800583   1.740471   0.9574559  0.9513467 
    ## clustCoef1   0.2816644  0.303749   0.2645848  0.3416067  0          0         
    ## modularity1  0.3990592  0.5024414  0.3388581  0.4091975  -0.125     -0.125    
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.02746135 0.05677337 0.02163896 0.03569209 0.5535814  0.5544811 
    ## density1     0.1003401  0.1264822  0.08704557 0.1280228  0.6666667  0.6666667 
    ## pep1         50         56.25      33.52941   55.55556   50         50
