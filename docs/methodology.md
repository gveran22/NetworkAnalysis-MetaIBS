# 📘 Methodology

This document outlines the methodology for preprocessing, analyzing, and comparing microbial networks in the **NetworkAnalysis-MetaIBS** project.

---

## 1. 🧹 Data Preprocessing

### Structure
- Raw data is stored in the `data/` folder. Each file contains a `phyloseq` object for one study (ASV table + metadata).
- For dataset origins, refer to the [MetaIBS project](https://github.com/bio-datascience/MetaIBS).

### Steps

1. **Merge**
   - All study-specific datasets are merged into a single `phyloseq` object.

2. **Agglomeration**
   - The merged dataset is agglomerated at various taxonomy levels (e.g., Class, Family, Order).
   - Missing taxonomy labels are cleaned and renamed during this step.

3. **Splitting Agglomerated Data**
   - The agglomerated dataset is split into subsets for downstream analysis:
     - **a. Individual datasets:** Split by study (`author` variable).
     - **b. Combined datasets:** Split by selected comparison variables (e.g., `sample_type`, `sequencing_tech`, etc.).

   These outputs serve as input for both the single-network and comparison-network inference pipelines.

🧾 Script:  
- `scripts/preprocessing.R`

📂 Outputs:  
- Saved in `build/Merge/` and `build/Agglomeration/`

---

## 2. 🔍 Filtering Exploration

- Various filter combination were tested and visualized to identify the most consistent approach across all study-specific datasets.
- Bar plots were used to evaluate outcomes.

📂 Outputs:  
- `outputs/investigation/`
- `outputs/filtering-investigation/`

🧾 Script:  
- `scripts/filtering_investigation.R`

---

## 3. 🧠 Network Inference

### 3.1 Single-Network Analysis
- Microbial networks are inferred separately for each dataset.
- Methods: SpiecEasi using three approaches: `glasso`, `mb` and `slr`.
- Metrics: Relative LCC size, clustering coefficient, modularity, and more.

📂 Outputs:  
- `outputs/single-network-analysis/` — includes association matrices, SpiecEasi results, plots, and network metrics

🧾 Script:  
- `scripts/single-network-analysis/run_analysis.R`

### 3.2 Network Comparison
- Preprocessing: The split agglomerated datasets are split into IBS and Healthy samples.
- Comparative networks are inferred for each individual dataset study and combined dataset.
- Methods: SpiecEasi with `glasso`, `mb` and `slr`.
- Metrics: Relative LCC size, clustering coefficient, modularity, etc.
  
📂 Outputs:  
- `outputs/network-comparison/Combined/`  
- `outputs/network-comparison/Individual/`

🧾 Script:  
- `scripts/network-comparison/run_comparison.R`

---

## 4. 📊 Meta-Analysis

- Aggregates and summarizes network results.
- Visual outputs include network plots.

📂 Reports:  
- Located in `docs/meta-analysis/`

🧾 Source Files:  
- Includes `.Rmd` and `.md` reports

---

## 5. 🧰 Tools

### Functions
- Custom functions used across the analysis are stored in:
  - `tools/functions.R`

### Configuration
- Filtering parameters: `tools/analysis_configs.R`  
- Grouping and comparison variables: `tools/analysis_variables.R`

---

## 📚 Packages Used

- `tidyverse` – data wrangling & plotting  
- `SpiecEasi` – network inference  
- `igraph`, `NetCoMi` – network metrics & visualization  
- `phyloseq`, `vegan` – ecological analyses

---

## 📖 References

1. Zachary D. Kurtz, Richard Bonneau, Christian L. Müller (2019). *Disentangling microbial associations from hidden environmental and technical factors via latent graphical models*. bioRxiv. [https://www.biorxiv.org/content/10.1101/2019.12.21.885889v1](https://www.biorxiv.org/content/10.1101/2019.12.21.885889v1)
2. Zachary D. Kurtz et al. (2015). *Sparse and compositionally robust inference of microbial ecological networks*. arXiv. [https://arxiv.org/abs/1408.4158](https://arxiv.org/abs/1408.4158)
3. Salomé Carcy et al. (2024). *MetaIBS - large-scale amplicon-based meta-analysis of irritable bowel syndrome*. bioRxiv. [https://www.biorxiv.org/content/10.1101/2024.01.22.575775v1](https://www.biorxiv.org/content/10.1101/2024.01.22.575775v1)

---


## 📬 Contact

For questions or contributions: [Gilary Vera Nuñez](mailto:gilary.vera22@gmail.com)