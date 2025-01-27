# Network Analysis - MetaIBS

This repository contains the data, scripts, and results for the **Network Analysis - MetaIBS** project. The study focuses on analyzing microbial networks using taxonomy-level aggregations and comparing network structures across datasets and experimental conditions.

---

## **Repository Structure**

```plaintext
.
├── data/                     # Raw input data
│   └── Individual/           # Original input files
├── build/                    # Intermediate processed data
│   ├── Agglomeration/        # Aggregated data by taxonomy level
│   └── Combined/             # Preprocessed data for analyses
├── outputs/                  # Key results and exploratory outputs
│   ├── investigation/        # Filtering and exploratory results
│   ├── network-comparison/   # Final network comparison outputs
│   │   ├── Combined/
│   │   └── Individual/
│   └── single-network-analysis/ # Outputs for single network analysis
├── scripts/                  # Code for analyses
│   ├── network-comparison/   # Pipeline for network comparison
│   └── single-network-analysis/  # Pipeline for single-network analysis
├── docs/                     # Documentation and reports
│   ├── meta-analysis/        # `.Rmd` and `.md` files with plots and analysis
│   └── methodology.md        # Filtering amd methodology details
├── tools/                    # Reusable functions and utilities
│   └── functions.R
├── LICENSE                   # License file
└── README.md                 # Project overview

```

<br/>

## **Installation**
To replicate the analysis, you will need:

1. R (version >= 4.1.0)
2. Required R packages:
  - tidyverse
  - igraph
  - vegan
  - ggplot2

Clone the repository:  
`git clone https://github.com/gveran22/NetworkAnalysis-MetaIBS.git  
cd NetworkAnalysis-MetaIBS`

## **Usage**
1. Preprocessing Data
   - Place your raw input files in the `data/Individual/` folder.
   - Run the preprocessing script in `scripts/`:  
        [Contribution guidelines for this project] (scripts/single-network-analysis/preprocessing.R)

2. Single-Network Analysis
   - Analyze individual datasets using the `single-network-analysis` scripts:  
Rscript scripts/single-network-analysis/run_analysis.R

3. Network Comparison
   - Compare networks using the `network-comparison` pipeline:
Rscript scripts/network-comparison/run_comparison.R
4. Meta-Analysis
   - Review the final meta-analysis results in `docs/meta-analysis/` 

## **Folder Details**
**[data/](data/)**
- `Individual/`: Raw input files for the project.  
**[build/](build/)**
- `Agglomeration/`: Aggregated datasets by taxonomy level.
- `Combined/`: Merged and preprocessed datasets.  
**[outputs/](outputs/)**
- `investigation/`: Results from exploratory filtering and testing.
- `network-comparison/`: Final results from network comparison analyses.
   - `Combined/`: Results from merged datasets.
   - `Individual/`: Results from individual datasets.
- `single-network-analysis/`: Results from analyzing individual datasets.  
**[scripts/](scripts/)**
- `network-comparison/`: Scripts to process and compare networks across datasets.
- `single-network-analysis/`: Scripts to analyze individual datasets.  
**[docs/](docs/)**
- `meta-analysis/`: Reports and visualizations summarizing the findings.
- `methodology.md`: Details on data filtering and processing steps.  

## **Contributing**
Contributions are welcome! Please:
1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request with a detailed description of changes.

## **License**
This project is licensed under the MIT License. See the LICENSE file for details.

## **Contact**
For questions or suggestions, contact Your Name.




