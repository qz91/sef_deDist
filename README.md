# Repository for code review

## Data Dictionary for *Differential Density Analysis in Single-Cell Genomics Using Specially Designed Exponential Families*

We provide information on the structure of downstream analysis and reproducibility of figures for the manuscript here.

--- 

### 1. Data Overview
#### Single Cell Data

The data used in our analysis for each cell type can be accessed via Zenodo using the following [link](https://zenodo.org/records/17402494?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjM0MzNkNDE1LWY4ZTctNDVhYi1hODk5LWJmNzhjNzg4MDUxNyIsImRhdGEiOnt9LCJyYW5kb20iOiI1OTVhOGVjZTBkYmZkZjBjMDA2ZTY4ZTBmNmVjN2Q3NiJ9.ONISAR5Zgx5GZ0odRZKmfSKmKTzBUTRyZ250S-hCc18EzXopSVeq12rdOqvJt_VgHZaHObG8x909Sya_aV9CVQ).

---

### 2. File Overview
We briefly explain the file purpose in our data analysis and simulations, which are split into the `RDA` and `simulations` directories in the repository.

#### Simulations
The simulations section is split into separate `R` files  for distinct functions. First, we have a main file that defines all functions used in simulations. We then have files generate the expression data. In addition to these core files, we include auxiliary files that run the SEF tests for the mean, variance, and modality shift settings. We include the simulation schemes for competing methods for pseudobulk and cell-level methods. Lastly, we include supplementary codes as well, including tests involving discrete carriers, bin size sensitivity, and comparisons with a direct moment test.

| File Name | Format     | Description |
|---------------|----------|------------|
| `lpc_simulation_main_5_8_2026`  | `.R`   | `main` file defining functions used in simulations |
| `lpc_revisions_simulation_discrete_6_7_2026`  | `.R`   | `main` file defining discrete-specific functions used in simulations |
| `lpc_main_NB_4_30_2026`  | `.R`   | `main` file used to define data  generation for `gen_sim_obj_EFFECT` scripts |
| `gen_sim_obj_EFFECT`  | `.R`   | generate simulated expression data based on effect size (modal_NB for modalitiy, strong_effect for mean, small_effect for variance) |
| `lpc_generalized_sef_omnicarrier`  | `.R`   | Workflow used on computing cluster to conduct SEF framework in empirical FDR test for mean and variance shift |
| `lpc_pval_combination_generalized_(ss)_final`  | `.R`   | P-value combination test at varying sample sizes |
| `lpc_NB_modality_analysis_final`  | `.R`   | `main` file defining modality-specific functions used in simulations  |
| `lpc_NB_modal_scaled_SEF_analysis_final`  | `.R`   | Workflow used on computing cluster to conduct SEF framework in empirical FDR test for modality shift |
---

#### Real Data Analysis
In the real data analysis section, we provide code to demonstrate the main real data contributions using the proposed SEF methodology as well as replicate the key figures with respect to these data.
This includes the regression pipeline and enrichment figures.


| File Name | Format     | Description |
|---------------|----------|------------|
| `merged_newCov_revisions_RDA_final`      | `.R`  | Runs cell-type specific SEF regression modeling and testing; provides general enrichment analyses across all cell types as well as unique CD8+ enrichment |
| `revisions_pseudobulk_comparison_5_5_2026`  | `.R`   | Procedural pipeline for pseudobulk methods |
| `merged_self_contained_RDA_final`  | `.R`   | self-contained `main` file defining functions used in real data analysis and downstream analysis |
| `run_pseudobulk_final`  | `.R`   | Script to generate results from comparison pseudobulk analysis |
| `prepare_pseudobulk_final`  | `.R`   | Prepare pseudobulk data |
| `bootstrap_RDA`  | `.R`   | Workflow to reproduce bootstrap analysis of CD8+ T-cells. Results can be found in the subdirectory `RDA/aux_response/cd8_bootstrap_stability` |
| `b_cells_explore`  | `.R`   | Basic QC metrics for B cells as requested in revisions|
| ``  | `.R`   | Workflow to reproduce Figure 3; Figures S.5, S.6, S.7 |

---

### 3. Object Overview
In this section, we overview the formatted data objects used in our simulations and real data analysis that are found on [Zenodo](https://zenodo.org/records/17402494?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjM0MzNkNDE1LWY4ZTctNDVhYi1hODk5LWJmNzhjNzg4MDUxNyIsImRhdGEiOnt9LCJyYW5kb20iOiI1OTVhOGVjZTBkYmZkZjBjMDA2ZTY4ZTBmNmVjN2Q3NiJ9.ONISAR5Zgx5GZ0odRZKmfSKmKTzBUTRyZ250S-hCc18EzXopSVeq12rdOqvJt_VgHZaHObG8x909Sya_aV9CVQ).
`CELLTYPE` can be exchanged for the specific cell types used in our analysis.


| Object Name | Format     | Description |
|---------------|----------|------------|
| `covariate_data`           | `.csv`     | Covariate data with donor ID and disease status specification |
| `CELLTYPE_lcf_matrix` | `.RDS`  | Library-size corrected counts matrix (genes x cell) |
| `CELLTYPE_new_metadata`| `.RDS`   | Contains metadata from original data (cell-level rows) |
| `CELL_TYPE_donors_used`         | `.RDS`   | Donors used in each cell type-specific analyses |
| `CELLTYPE_enrichment_07_01_2026`       | `.csv`   | `clusterProfiler` enrichment analysis results; used to re-construct figures |
| `COMPARISON_CELLTYPE_sig_pvals`         | `.csv`   | Significant genes identified by a compared method (DESeq2, edgeR, Wilcoxon, MAST) after Bonferroni-correction |
| `pb_CELLTYPE_matrix`         | `.RDS`   | Pseudobulked expression matrix after using standard `AggregateExpression()` |
| `pb_CELLTYPE_metadata`         | `.RDS`   | Pseudobulk metadata |

---
