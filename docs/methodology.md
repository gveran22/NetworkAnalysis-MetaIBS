# Methodology

This document outlines the methodologies used for preprocessing, analyzing, and comparing microbial networks in the **Network Analysis - MetaIBS** project.

---

## **1. Data Preprocessing**

### 1.1 Raw Data Structure
- Raw input data is stored in the `data/Individual/` folder.
- Each dataset contains:
  - Taxonomic abundance tables.
  - Metadata associated with experimental conditions.

### 1.2 Preprocessing Steps
1. **Normalization:**
   - Abundance values are normalized using relative abundance.
   - Scripts: `scripts/single-network-analysis/preprocessing.R`
2. **Filtering:**
   - Operational Taxonomic Units (OTUs) with low prevalence (<5% occurrence across samples) are removed.
   - OTUs with zero variance are excluded.
3. **Agglomeration:**
   - Taxonomic levels aggregated into Class, Family, Order, etc., using custom functions in `tools/functions.R`.

Outputs:
- Preprocessed datasets are stored in `build/Agglomeration/`.

---

## **2. Network Analysis**

### 2.1 Single-Network Analysis
- Individual datasets are analyzed to construct microbial association networks.
- Methodology:
  - Correlation-based network construction (e.g., Spearman correlation).
  - Thresholding based on correlation strength (e.g., |r| > 0.6).
  - Network metrics calculated: Degree, Betweenness Centrality, Clustering Coefficient.
- Outputs:
  - Association matrices: `outputs/single-network-analysis/Individual/association_matrices/`
  - Filtered networks: `outputs/single-network-analysis/Individual/filtered_outputs/`
  - Final plots: `outputs/single-network-analysis/Individual/final_plots/`

Scripts:
- `scripts/single-network-analysis/run_analysis.R`

---

### 2.2 Network Comparison
- Comparison of networks constructed from combined datasets or across conditions.
- Methodology:
  - Pairwise comparison of network metrics (e.g., average degree, modularity).
  - Statistical significance tested using permutation tests.
- Outputs:
  - Combined and individual results: `outputs/network-comparison/`
  - Final plots: `outputs/network-comparison/final_plots/`

Scripts:
- `scripts/network-comparison/run_comparison.R`

---

## **3. Filtering and Exploration**

### 3.1 Investigation of Filtering
- Exploratory filtering methods were applied to test their impact on OTU tables.
- Methodology:
  - Comparison of retained OTUs across different filtering thresholds.
  - Visualization of filtering effects using diversity metrics.
- Outputs:
  - Exploratory results: `outputs/investigation/`

Scripts:
- `scripts/filter_exploration.R`

---

## **4. Meta-Analysis**

- Summary reports and visualizations of final results.
- Methodology:
  - Aggregation of network metrics across datasets.
  - Comparison across experimental conditions using visualizations (e.g., heatmaps, barplots).
- Reports:
  - Located in `docs/meta-analysis/`
  - Includes `.Rmd` and `.md` files for reproducibility.

---

## **5. Tools and Utilities**

### 5.1 Functions
- Custom R functions for preprocessing and network analysis.
- Key files:
  - `tools/functions.R`

### 5.2 Libraries Used
- `tidyverse`: Data manipulation and visualization.
- `igraph`: Network construction and metrics calculation.
- `vegan`: Diversity analysis.
- `ggplot2`: Plotting results.

---

## **References**

1. **Network Construction Methodology:** Barabási, A.-L. (2016). *Network Science*. Cambridge University Press.
2. **Microbial Network Analysis:** Faust, K., & Raes, J. (2012). Microbial interactions: From networks to models. *Nature Reviews Microbiology*.

---

## **Contact**

For additional details or questions, contact [Your Name](mailto:youremail@example.com).
