# Network Inference and Comparative Analysis of Gut Microbiomes in Healthy and IBS-Patients

This repository is an extension of the **MetaIBS project** ([link to project](https://github.com/bio-datascience/MetaIBS)). It contains the data, scripts, and results for the Inference and Network Analysis of the MetaIBs datasets. The study focuses on analyzing microbial networks using taxonomy-level aggregations and comparing network structures divided by the `healthy`status across datasets as well as the merged of the datasets based on some important variables.

---

## **Repository Structure**

```plaintext
.
├── data/                     # Raw input data
├── build/                    # Intermediate processed data
│   ├── Agglomeration/        # Aggregated data by taxonomy level
│   └── ...             
├── outputs/                  # Key results and exploratory outputs
│   ├── investigation/        # Filtering and exploratory results
│   ├── network-comparison/   # Final network comparison outputs
│   │   ├── Combined/
│   │   └── Individual/
│   └── single-network-analysis/ # Outputs for single network analysis
├── scripts/                  # Code for analyses
│   ├── network-comparison/   # Pipeline for network comparison
│   ├── single-network-analysis/  # Pipeline for single-network analysis
│   ├── filtering_investigation.R 
│   ├── preprocessing.R 
│   └── README.md
├── docs/                     # Documentation and reports
│   ├── meta-analysis/        # `.Rmd` and `.md` files with plots and analysis
│   └── methodology.md        # Filtering amd methodology details
├── tools/                    # Reusable functions and utilities
│   ├── functions.R
│   ├── analysis_configs.R 
│   └── anaylsis_variables.R
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

<br/>

## **How to use this repository**
1. Preprocessing Data
   - Place your raw input files in the [data/](data/) folder.
   - Use the preprocessing script in `scripts/` as a reference for preprocessing:
   everything you need to know is in the [scripts/preprocessing.R](scripts/preprocessing.R).
   - Save all your preprocessed data in the [build/](build/) folder.

2. Single-Network Analysis
   - Define your filtering parameters in the [tools/analysis_configs.R](tools/analysis_configs.R) file.  
   - Run the individual analysis on each dataset using the [scripts/single-network-analysis/](scripts/single-network-analysis/) scripts:   `scripts/single-network-analysis/run_analysis.R`  

3. Network Comparison
   - Define your variables parameters in the [tools/analysis_variables.R](tools/analysis_variables.R) file.  
   - Run the network comparison analysis for each dataset, the total merge and the variables you defined in the step before, using the [scripts/network-comparison/](scripts/network-comparison/) scripts:  scripts/network-comparison/run_comparison.R

4. Meta-Analysis
   - Review the final meta-analysis results in [docs/meta-analysis/](docs/meta-analysis/) adapting the files with your variables.

<br/>

## **Folder Details**
**[data/](data/)**: Raw input files for the project.

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

<br/>

## **Contributing**
Contributions are welcome! Please:
1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request with a detailed description of changes.

<br/>

## **License**
See the LICENSE file for details.

<br/>

## **Contact**
For questions or suggestions, contact .




