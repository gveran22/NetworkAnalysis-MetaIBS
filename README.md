# 🧬 Network Inference and Comparative Analysis of Gut Microbiomes in Healthy and IBS-Patients

This repository is an extension of the **MetaIBS project** ([link to project](https://github.com/bio-datascience/MetaIBS)). It contains the data, scripts, and results for the inference and network Analysis of the MetaIBs datasets. The study focuses on microbial network inference using taxonomy-level aggregations and compares network structures across `health` statuses (IBS vs. Healthy) within and across datasets.

---

## 📁 Repository Structure

```plaintext
.
├── data/                     # Raw input data
├── build/                    # Intermediate processed data
│   ├── Agglomeration/        # Aggregated data by taxonomy level
│   └── Merge/                # Merged datasets
├── outputs/                  # Key results and exploratory outputs
│   ├── investigation/        # Filtering and exploratory results
│   ├── network-comparison/   # Results of network comparison analysis
│   │   ├── Combined/
│   │   └── Individual/
│   └── single-network-analysis/ # Outputs from single network analysis
├── scripts/                  # Code for analyses
│   ├── network-comparison/   # Pipeline for network comparison
│   ├── single-network-analysis/  # Pipeline for single-network analysis
│   ├── filtering_investigation.R 
│   ├── preprocessing.R 
│   └── README.md
├── docs/                     # Documentation and reports
│   ├── meta-analysis/        # `.Rmd` and `.md` files with plots and summaries
│   └── methodology.md        # Filtering amd methodology details
├── tools/                    # Reusable functions and utilities
│   ├── functions.R
│   ├── analysis_configs.R 
│   └── anaylsis_variables.R
└── README.md                 # Project overview (this file)

```

<br/>

##⚙️ Installation 
To replicate the analysis, you will need base R, SpiecEasi and NetComi. See the [SpiecEasi github](https://github.com/zdk123/SpiecEasi) and the [NetComi github](https://github.com/stefpeschel/NetCoMi) pages for more details. This involves a few auxiliary R packeges like: devtools and phyloseq

Clone the repository:  
`git clone https://github.com/gveran22/NetworkAnalysis-MetaIBS.git`  
`cd NetworkAnalysis-MetaIBS`

<br/>

## 🚀 How to use this repository 
**1. Preprocess Data**
   - Place your raw input files in the [data/](data/) folder.
   - Use [scripts/preprocessing.R](scripts/preprocessing.R) to preprocess and agglomerate data .
   - Save your outputs to the [build/](build/) directory.
   
**2. Defining filter parameters**
   - Edit filtering settings in [tools/analysis_configs.R](tools/analysis_configs.R).
   - If unsure of optimal values, run [scripts/filtering_investigation.R](scripts/filtering_investigation.R) to explore them.
   
**3. Single-Network Analysis**
   - Run the pipeline on each dataset using scripts in [scripts/single-network-analysis/](scripts/single-network-analysis/)   
      - Entry point: `scripts/single-network-analysis/process.R`
**4. Network Comparison**
   - Set anaylsis parameters in [tools/analysis_variables.R](tools/analysis_variables.R).  
   - Run comparisons across merged and individual datasets using scripts in [scripts/network-comparison/](scripts/network-comparison/) scripts:  scripts/network-comparison/run_comparison.R
      - Entry point: `scripts/network-comparison/process_comparison.R`
**5. Meta-Analysis**
   - Explore and adapt summary plots and insights in [docs/meta-analysis/](docs/meta-analysis/).

<br/>

## **Folder Details**
**[data/](data/)**: Raw input files for the project.

**[build/](build/)** 
- `Agglomeration/`: Datasets aggregated at various taxonomy levels.
- `Merge/`: Merge phyloseq dataset.  

**[outputs/](outputs/)**
- `investigation/`: Results from exploratory filtering and testing. 
- `network-comparison/`: Results of IBS vs. Healthy comparisons.
   - `Combined/`: Across merged datasets.
   - `Individual/`: Separate per datasets.
- `single-network-analysis/`: Results from individual dataset analysis. 

**[scripts/](scripts/)**
- `network-comparison/`: Scripts to process and compare networks across merged and individual datasets.
- `single-network-analysis/`: Code for individual dataset analysis.  

**[docs/](docs/)**
- `meta-analysis/`: Reports and visualizations of findings.
- `methodology.md`: Filtering strategy and pipeline overview.  

<br/>

## 📬 Contact
For questions or suggestions, contact [Gilary Vera Nunez](mailto:gilary.vera22@gmail.com).




