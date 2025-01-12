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

    ##              agp       fukui     hugerth   labus      liu        lopresti 
    ## lccSize1     33        23        19        3          70         3        
    ## lccSizeRel1  0.3882353 0.6052632 0.4042553 0.2        0.7142857  0.2142857
    ## avDiss1      0.6588862 0.676858  0.6647847 0.713513   0.6839907  0.720823 
    ## avPath1      1.119195  1.708288  1.248548  0.713513   1.911724   0.9610973
    ## clustCoef1   0.7751095 0.5784896 0.745734  1          0.5215366  0        
    ## modularity1  0.151812  0.3466446 0.1569471 -0.2222222 0.45132    -0.125   
    ## vertConnect1 1         1         1         2          1          1        
    ## edgeConnect1 1         1         1         2          1          1        
    ## natConnect1  0.1427939 0.0649974 0.1058802 0.5723413  0.02503974 0.552947 
    ## density1     0.5018939 0.1818182 0.3567251 1          0.1035197  0.6666667
    ## pep1         98.49057  73.91304  91.80328  33.33333   83.6       0        
    ##              mars       nagel     pozuelo    zeber      zhu       zhuang   
    ## lccSize1     33         11        31         23         7         10       
    ## lccSizeRel1  0.7857143  0.4583333 0.5849057  0.6388889  0.3333333 0.5      
    ## avDiss1      0.6935829  0.7046715 0.6853521  0.7004371  0.7072553 0.7017804
    ## avPath1      2.142578   1.379689  1.273078   1.572758   1.239689  1.66335  
    ## clustCoef1   0.4114279  0.6198881 0.7775831  0.5322126  0.6030636 0.2514602
    ## modularity1  0.4899554  0.183642  0.1227822  0.2194083  0.25      0.315    
    ## vertConnect1 1          1         1          1          1         1        
    ## edgeConnect1 1          1         1          1          1         1        
    ## natConnect1  0.04009346 0.1361131 0.07794524 0.07007887 0.218496  0.1416859
    ## density1     0.1060606  0.3272727 0.2946237  0.256917   0.3809524 0.2222222
    ## pep1         73.21429   33.33333  67.15328   50.76923   50        60

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-mb-12.png)<!-- -->

### Global Properties

    ##              agp        fukui     hugerth    labus     liu        lopresti 
    ## lccSize1     77         29        35         3         97         3        
    ## lccSizeRel1  0.9058824  0.7631579 0.7446809  0.2       0.9897959  0.2142857
    ## avDiss1      0.6783201  0.6877272 0.6780449  0.7178617 0.676574   0.725649 
    ## avPath1      2.408119   2.392792  2.219635   0.9571489 2.078459   0.9675321
    ## clustCoef1   0.3624016  0.3158229 0.4549617  0         0.1882751  0        
    ## modularity1  0.5553269  0.4833432 0.5240055  -0.125    0.4856135  -0.125   
    ## vertConnect1 1          1         1          1         1          1        
    ## edgeConnect1 1          1         1          1         1          1        
    ## natConnect1  0.01761718 0.0450091 0.03760382 0.5538638 0.01458214 0.5520626
    ## density1     0.05536569 0.1009852 0.0907563  0.6666667 0.06164089 0.6666667
    ## pep1         89.50617   73.17073  83.33333   50        79.79094   0        
    ##              mars       nagel     pozuelo   zeber      zhu       zhuang   
    ## lccSize1     41         8         37        24         11        11       
    ## lccSizeRel1  0.9761905  0.3333333 0.6981132 0.6666667  0.5238095 0.55     
    ## avDiss1      0.6918297  0.704997  0.6925744 0.6962675  0.7097545 0.7010129
    ## avPath1      2.25378    1.458501  1.731306  1.813421   2.055932  1.68932  
    ## clustCoef1   0.1538931  0.4062545 0.4987073 0.4256725  0         0.2452739
    ## modularity1  0.5280473  0.2160494 0.299311  0.4044817  0.415     0.3347107
    ## vertConnect1 1          1         1         1          1         1        
    ## edgeConnect1 1          1         1         1          1         1        
    ## natConnect1  0.03109962 0.1865406 0.0382905 0.05855172 0.1251415 0.1269083
    ## density1     0.07926829 0.3214286 0.1186186 0.1702899  0.1818182 0.2      
    ## pep1         67.69231   44.44444  70.88608  57.44681   40        54.54545

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Order/meta-analysis-slr-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Order/single-network-slr-12.png)<!-- -->

### Global Properties

    ##              agp        fukui     hugerth    labus     liu        lopresti 
    ## lccSize1     45         19        24         2         64         3        
    ## lccSizeRel1  0.5294118  0.5       0.5106383  0.1333333 0.6530612  0.2142857
    ## avDiss1      0.7016428  0.7062351 0.6984021  0.7052893 0.7073743  0.7228231
    ## avPath1      1.455371   1.551683  1.348668   0.7052893 1.779806   0.9637642
    ## clustCoef1   0.4830784  0.2272882 0.4604959  0         0.3991052  0        
    ## modularity1  0.2915174  0.2669753 0.263949   0         0.2022497  -0.125   
    ## vertConnect1 1          1         1          1         1          1        
    ## edgeConnect1 1          1         1          1         1          1        
    ## natConnect1  0.04003814 0.0736874 0.06337483 0.7979573 0.03342331 0.552534 
    ## density1     0.1989899  0.2105263 0.2463768  1         0.1532738  0.6666667
    ## pep1         47.71574   58.33333  67.64706   100       45.30744   0        
    ##              mars       nagel     pozuelo   zeber      zhu       zhuang   
    ## lccSize1     16         3         8         21         3         2        
    ## lccSizeRel1  0.3809524  0.125     0.1509434 0.5833333  0.1428571 0.1      
    ## avDiss1      0.7026137  0.7108211 0.7188227 0.7076238  0.7072451 0.7432025
    ## avPath1      2.207201   0.9477615 1.750851  1.864856   0.9429934 0.7432025
    ## clustCoef1   0          0         0         0.1607879  0         0        
    ## modularity1  0.4984568  -0.125    0.3367347 0.4049587  -0.125    0        
    ## vertConnect1 1          1         1         1          1         1        
    ## edgeConnect1 1          1         1         1          1         1        
    ## natConnect1  0.08303447 0.555117  0.1805632 0.06337917 0.5555684 0.790155 
    ## density1     0.15       0.6666667 0.25      0.1571429  0.6666667 1        
    ## pep1         77.77778   50        42.85714  54.54545   50        0
