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

    ##              agp        fukui      hugerth    labus      liu        lopresti  
    ## lccSize1     65         94         71         65         92         82        
    ## lccSizeRel1  0.65       0.94       0.71       0.65       0.92       0.82      
    ## avDiss1      0.6835209  0.695509   0.6820644  0.6536894  0.6809376  0.6559947 
    ## avPath1      1.865869   2.092631   1.762751   3.861006   1.871264   3.588531  
    ## clustCoef1   0.648898   0.2981429  0.5875059  0.4213358  0.4283684  0.435424  
    ## modularity1  0.4184792  0.4125681  0.394517   0.7202376  0.425201   0.7000781 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.03144775 0.01608005 0.02756466 0.02040239 0.02334953 0.01565434
    ## density1     0.11875    0.07115077 0.1191147  0.05144231 0.09412327 0.03884372
    ## pep1         80.16194   64.63023   67.90541   99.06542   68.52792   96.89922  
    ##              mars       nagel      pozuelo    zeber      zhu        zhuang    
    ## lccSize1     94         94         92         93         82         81        
    ## lccSizeRel1  0.94       0.94       0.92       0.93       0.82       0.81      
    ## avDiss1      0.6806193  0.694047   0.6953548  0.6912623  0.6969977  0.6899053 
    ## avPath1      2.109396   2.533343   2.089248   2.032162   2.360882   2.708451  
    ## clustCoef1   0.3692454  0.2772837  0.4381216  0.2739537  0.4874564  0.3430178 
    ## modularity1  0.3977615  0.578405   0.3707372  0.3825702  0.4045063  0.5954232 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.02333769 0.01435857 0.01900541 0.01746142 0.02302316 0.01650385
    ## density1     0.07870053 0.05101807 0.08241758 0.07246377 0.08611864 0.05308642
    ## pep1         66.56977   63.2287    62.89855   60.32258   55.59441   66.27907

## MB

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-mb-1.png)<!-- -->
\### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-mb-12.png)<!-- -->
\### Global Properties

    ##              agp        fukui      hugerth    labus      liu        lopresti  
    ## lccSize1     90         97         85         98         100        99        
    ## lccSizeRel1  0.9        0.97       0.85       0.98       1          0.99      
    ## avDiss1      0.6882947  0.6963623  0.6901529  0.6611633  0.6915438  0.6544541 
    ## avPath1      2.352024   2.106495   2.257746   3.074277   1.925555   2.765108  
    ## clustCoef1   0.466776   0.168276   0.4109064  0.09816285 0.2183539  0.06348188
    ## modularity1  0.5549469  0.4500137  0.5190794  0.6221759  0.4263449  0.5665249 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.01623125 0.01378131 0.01666941 0.01275178 0.01442917 0.01276665
    ## density1     0.0639201  0.05648625 0.06330532 0.03324216 0.06646465 0.03360132
    ## pep1         86.32812   65.77947   76.54867   88.60759   68.38906   89.57055  
    ##              mars       nagel      pozuelo    zeber      zhu        zhuang    
    ## lccSize1     100        100        96         100        100        98        
    ## lccSizeRel1  1          1          0.96       1          1          0.98      
    ## avDiss1      0.6928586  0.6919158  0.6949348  0.6947479  0.6850989  0.6898165 
    ## avPath1      2.083508   2.271042   2.203649   1.997583   2.418865   2.671867  
    ## clustCoef1   0.1668661  0.07563253 0.2304697  0.1264326  0.1427558  0.0786478 
    ## modularity1  0.4620999  0.4717758  0.4355174  0.3864634  0.4854122  0.5430659 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.0134501  0.01273407 0.01458595 0.01372376 0.01300203 0.01266535
    ## density1     0.05494949 0.04464646 0.06030702 0.06020202 0.04626263 0.03745003
    ## pep1         67.27941   57.46606   67.63636   60.40268   61.57205   57.86517

## SLR

### Meta-Analysis Plot

![](../../../outputs/single-network-analysis/Individual/plots/Genus/meta-analysis-slr-1.png)<!-- -->

### Individual Plots

![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-1.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-2.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-3.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-4.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-5.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-6.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-7.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-8.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-9.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-10.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-11.png)<!-- -->![](../../../outputs/single-network-analysis/Individual/plots/Genus/single-network-slr-12.png)<!-- -->

### Global Properties

    ##              agp       fukui      hugerth   labus      liu        lopresti  
    ## lccSize1     11        79         10        72         84         99        
    ## lccSizeRel1  0.11      0.79       0.1       0.72       0.84       0.99      
    ## avDiss1      0.6815851 0.7065981  0.7135455 0.705368   0.7030791  0.6832116 
    ## avPath1      2.386674  2.370799   1.916723  1.958064   1.91386    2.084529  
    ## clustCoef1   0         0.1074589  0         0.1957753  0.1605811  0.2364306 
    ## modularity1  0.435     0.4780866  0.355     0.3134354  0.3671598  0.4516257 
    ## vertConnect1 1         1          1         1          1          1         
    ## edgeConnect1 1         1          1         1          1          1         
    ## natConnect1  0.1269472 0.01563605 0.1404784 0.01959148 0.01600959 0.01368478
    ## density1     0.1818182 0.04836092 0.2222222 0.08646322 0.06798623 0.05421563
    ## pep1         90        51.00671   60        49.77376   48.52321   78.327    
    ##              mars       nagel      pozuelo    zeber      zhu        zhuang    
    ## lccSize1     81         87         82         87         77         76        
    ## lccSizeRel1  0.81       0.87       0.82       0.87       0.77       0.76      
    ## avDiss1      0.7069709  0.7019635  0.7011657  0.6983703  0.7030747  0.7081659 
    ## avPath1      2.285801   2.148999   1.84978    2.310498   1.873436   1.937911  
    ## clustCoef1   0.1824925  0.29036    0.3682576  0.1374599  0.3563911  0.1736294 
    ## modularity1  0.467719   0.4735     0.4301383  0.4671837  0.3842596  0.3710585 
    ## vertConnect1 1          1          1          1          1          1         
    ## edgeConnect1 1          1          1          1          1          1         
    ## natConnect1  0.01556149 0.01504572 0.01748458 0.01445421 0.01887451 0.01737593
    ## density1     0.05308642 0.0553328  0.07588076 0.04677894 0.08475735 0.06842105
    ## pep1         43.60465   55.07246   61.11111   50.85714   52.01613   47.69231
