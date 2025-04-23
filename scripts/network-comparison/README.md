# 🧰 Project network-comparison

This folder contains all the scripts for the network inference pipeline used in the Network Comparison analysis. 

---

## 📄 Files

| **File**                              | **Short Description**                                     |
|---------------------------------------|------------------------------------------------------------|
| [`preprocessing_comparison.R`](./preprocessing_comparison.R) | Splits agglomerated phyloseq objects by the `host_disease` variable (`IBS` vs `Healthy`) |
| [`filtering_comparison.R`](./filtering_comparison.R) | Applies filtering to the split phyloseq objects |
| [`fit_network_comparison.R`](./fit_network_comparison.R) | Infers networks using SpiecEasi pipeline and saves the resulting objects |
| [`assoc_mat_comparison.R`](./assoc_mat_comparison.R)  | Extracts association matrices from SpiecEasi objects |
| [`net_properties_comparison.R`](./net_properties_comparison.R)  | Computes network properties from association matrices |
| [`process_comparison.R`](./process_comparison.R)  | Runs the complete single-network inference pipeline  |
| [`venn-diagramm_comparison.R`](./venn-diagramm_comparison.R)  | Generates Venn diagrams to compare overlaps between methods   |
| [`README.md`](./README.md)           | This documentation file    

---

## 📝 File Details

- `filtering.R`: Applies filtering steps to a phyloseq object and saves the output to `outputs/network-comparison/filtered_otus/`.
- `fit_network.R`: Runs the SpiecEasi pipeline using three methods (`glasso`, `mb`, `slr`). Outputs are saved in `outputs/network-comparison/spiec-easi-results/`. Parameters for SpiecEasi can be modified in this script.
- `process.R`: Executes the full pipeline, including filtering and network inference. Filtering parameters are loaded from `tools/analysis_configs`.
- `venn-diagramm.R`: Generates Venn diagrams to compare edge overlaps across the different network inference methods. This script is standalone and should be run **after the inference is complete**, using the association matrices as input.