# Repository for code review

## Data Dictionary for *Differential Density Analysis in Single-Cell Genomics Using Specially Designed Exponential Families*

We provide information on the structure of downstream analysis and reproducibility of figures for the manuscript here.

--- 

### 1. Data Overview
#### Single Cell Data

The data used in our analysis for each cell type can be accessed via Zenodo using the following [link](https://zenodo.org/records/17402494?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjM0MzNkNDE1LWY4ZTctNDVhYi1hODk5LWJmNzhjNzg4MDUxNyIsImRhdGEiOnt9LCJyYW5kb20iOiI1OTVhOGVjZTBkYmZkZjBjMDA2ZTY4ZTBmNmVjN2Q3NiJ9.ONISAR5Zgx5GZ0odRZKmfSKmKTzBUTRyZ250S-hCc18EzXopSVeq12rdOqvJt_VgHZaHObG8x909Sya_aV9CVQ).

---

### 2. File Overview
We briefly explain the file purpose in our data analysis and simulations, which are split into the `RDA` and `simulations` directories in the repository. In each file, we provide a brief description for each of the functions used.

#### Simulations
The simulations section is split into separate `R` files  for distinct functions. First, we have a main file that defines all functions used in simulations. We then have files generate the expression data. In addition to these core files, we include auxiliary files that run the SEF tests for the mean, variance, and modality shift settings. We include the simulation schemes for competing methods for pseudobulk and cell-level methods. Lastly, we include supplementary codes as well, including tests involving discrete carriers, bin size sensitivity, and comparisons with a direct moment test.

| File Name | Format     | Description |
|---------------|----------|------------|
| `simulations_main_10_13_2025`  | `.R`   | `main` file defining functions used in simulations |
| `refined_simulations_10_13_2025`  | `.R`   | Workflow for Figure 1: QQ-Plots and power analysis in main contents for Poisson-Gamma and ZINB |
| `gen_sim_obj_EFFECT`  | `.R`   | generate simulated expression data based on effect size |
| `refined_supp_simulations_10_13_2025`  | `.R`   | Supplementary tests and figure construction code (Figures S.1-S.4) |
---

#### Real Data Analysis
In the real data analysis section, we provide code to demonstrate the main real data contributions using the proposed SEF methodology as well as replicate the key figures with respect to these data.
This includes the regression pipeline and enrichment figures.


| File Name | Format     | Description |
|---------------|----------|------------|
| `merged_newCov_revisions_RDA_final`      | `.R`  | Runs cell-type specific SEF regression modeling and testing; provides enrichment analyses |
| `revisions_pseudobulk_comparison_5_5_2026`  | `.R`   | Procedural pipeline for pseudobulk methods |
| ``  | `.R`   | Workflow to reproduce Table 1 numbers and Figure 2 |
| ``  | `.R`   | Workflow to reproduce proportions in Table 2 |
| ``  | `.R`   | Workflow to reproduce Figure 3; Figures S.5, S.6, S.7 |
| `merged_self_contained_RDA_final`  | `.R`   | self-contained `main` file defining functions used in real data analysis and downstream analysis |

---

### 3. Object Overview
In this section, we overview the formatted data objects used in our simulations and real data analysis that are found on [Zenodo](https://zenodo.org/records/17402494?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjM0MzNkNDE1LWY4ZTctNDVhYi1hODk5LWJmNzhjNzg4MDUxNyIsImRhdGEiOnt9LCJyYW5kb20iOiI1OTVhOGVjZTBkYmZkZjBjMDA2ZTY4ZTBmNmVjN2Q3NiJ9.ONISAR5Zgx5GZ0odRZKmfSKmKTzBUTRyZ250S-hCc18EzXopSVeq12rdOqvJt_VgHZaHObG8x909Sya_aV9CVQ).
`CELLTYPE` can be exchanged for the specific cell types used in our analysis.


| Object Name | Format     | Description |
|---------------|----------|------------|
| `covariate_data`           | `.csv`     | Covariate data with donor ID and disease status specification |
| `CELLTYPE_lcf_matrix` | `.RDS`  | Library-size corrected counts matrix (genes x cell) |
| `CELLTYPE_new_metadata`| `.RDS`   | Contains metadata from original data (cell-level rows) |
| `pb_CELLTYPE_sObj`         | `.RDS`   | Pseudobulked Seurat object after using standard `AggregateExpression()` |
| ``         | `.csv`   | Cell type-specific differentially distributed genes with adjusted p-values (Bonferroni) using SEF regression |
| ``         | `.csv`   | Cell type specific DEGs with adjusted p-values (Bonferroni) pseudobulked data and log fold change values |
| ``       | `.csv`   | Includes `clusterProfiler` enrichment analysis results and used to re-construct barplots from original analysis |

---
