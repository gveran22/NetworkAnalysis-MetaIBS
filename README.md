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
├── tools/                    # Reusable functions and configs
│   ├── functions.R
│   ├── analysis_configs.R 
│   └── anaylsis_variables.R
└── README.md                 # Project overview (this file)

```

<br/>

## ⚙️ Installation
To replicate the analysis, you will need **R**, along with the packages **SpiecEasi** and **NetCoMi**. For installation and details, refer to their GitHub repositories:

- [SpiecEasi GitHub](https://github.com/zdk123/SpiecEasi)
- [NetCoMi GitHub](https://github.com/stefpeschel/NetCoMi)


These packages depend on additional tools such as `devtools`, `phyloseq`, and others, which can be installed from CRAN or GitHub.

To clone this repository:  

`git clone https://github.com/gveran22/NetworkAnalysis-MetaIBS.git`  
`cd NetworkAnalysis-MetaIBS`

<br/>

## 🚀 How to Use

### 1. **Preprocess Data**
- Place your raw input files in the [`data/`](data/) folder.
- Use [`scripts/preprocessing.R`](scripts/preprocessing.R) to:
  - Merge datasets.
  - Aggregate by taxonomy level.
  - Split agglomerated data.
- Save the resulting files in the [`build/`](build/) directory.
   
### 2. **Define Filtering Parameters**
- Set your filtering thresholds and options in [`tools/analysis_configs.R`](tools/analysis_configs.R).
- If you're unsure of the best parameters, run the exploratory script [`scripts/filtering_investigation.R`](scripts/filtering_investigation.R) to test different configurations.
   
### 3. **Single-Network Analysis**
- Analyze each dataset individually using the scripts in [`scripts/single-network-analysis/`](scripts/single-network-analysis/).
- Execute the full pipeline with `process_analysis.R`.
      
### 4. **Network Comparison**
- Set comparison variables in [`tools/analysis_variables.R`](tools/analysis_variables.R).
- Run the comparative analysis using the scripts in [`scripts/network-comparison/`](scripts/network-comparison/), particularly `process_comparison.R`.
- This includes comparisons of IBS vs Healthy across datasets or merged data according the comparison variables.
      
### 5. **Meta-Analysis and Visualization**
- Review the meta-analysis outputs in [`docs/meta-analysis/`](docs/meta-analysis/).
- Customize `.Rmd` reports with your results and variable definitions for final reporting and publication-ready visualizations.


<br/>

## 📂 Folder Details

### [`data/`](data/)
Contains raw input data such as phyloseq objects.

### [`build/`](build/)
- `Agglomeration/`: Taxonomy-level aggregated phyloseq data.
- `Merge/`: Combined data from all datasets.


### [`outputs/`](outputs/)
- `investigation/`: Results from filtering exploration.
- `network-comparison/`: Network comparison results.
  - `Combined/`: Merged dataset results.
  - `Individual/`: Dataset-specific results.
- `single-network-analysis/`: Outputs from individual network inference pipelines.

### [`scripts/`](scripts/)
- `network-comparison/`: All scripts for comparative network inference.
- `single-network-analysis/`: Scripts for analyzing each dataset independently.
- `filtering_investigation.R`: Explore optimal filtering strategies.
- `preprocessing.R`: Handles merging, agglomeration, and splitting of datasets.

### [`docs/`](docs/)
- `meta-analysis/`: Final reports and visualizations.
- `methodology.md`: Filtering strategy and pipeline overview.

### [`tools/`](tools/)
- `functions.R`: Reusable functions for the pipeline and meta-analysis.
- `analysis_configs.R`: Filtering parameters for phyloseq objects.
- `analysis_variables.R`: Variables for network comparison setup.

<br/>

## 📬 Contact

For questions, feedback, or contributions, please reach out to:

**Gilary Vera Nuñez**  
📧 [gilary.vera22@gmail.com](mailto:gilary.vera22@gmail.com)
