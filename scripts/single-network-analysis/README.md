# 🧰 Project single-network-analysis

This folder contains all the scripts for the network inference pipeline used in the single network analysis. 

---

## 📄 Files

| **File**                              | **Short Description**                                     |
|---------------------------------------|------------------------------------------------------------|
| [`filtering.R`](./filtering.R) | Filters agglomerated phyloseq objects |
| [`fit_network.R`](./fit_network.R) | Infers networks using SpiecEasi pipeline and saves the resulting objects |
| [`assoc_mat.R`](./assoc_mat.R)  | Extracts association matrices from SpiecEasi objects |
| [`net_properties.R`](./net_properties.R)  | Computes network properties from association matrices |
| [`process.R`](./process.R)  | Runs the complete single-network inference pipeline  |
| [`venn-diagramm.R`](./venn-diagramm.R)  | Generates Venn diagrams to compare overlaps between methods   |
| [`README.md`](./README.md)           | This documentation file    

---

## 📝 File Details

- `filtering.R`: Applies filtering steps to a phyloseq object and saves the output to `outputs/single-network-analysis/filtered_otus/`.
- `fit_network.R`: Runs the SpiecEasi pipeline using three methods (`glasso`, `mb`, `slr`). Outputs are saved in `outputs/single-network-analysis/spiec-easi-results/`. Parameters for SpiecEasi can be modified in this script.
- `process.R`: Executes the full pipeline, including filtering and network inference. Filtering parameters are loaded from `tools/analysis_configs`.
- `venn-diagramm.R`: Generates Venn diagrams to compare edge overlaps across the different network inference methods. This script is standalone and should be run **after the inference is complete**, using the association matrices as input.