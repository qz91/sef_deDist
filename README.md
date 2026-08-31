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
The simulations section is split into three `R` files. We have a main file that defines all functions used in simulations. The remaining two files contain the main simulation results and that from the supplementary materials.

| File Name | Format     | Description |
|---------------|----------|------------|
| `simulations_main_10_13_2025`  | `.R`   | `main` file defining functions used in simulations |
| `refined_simulations_10_13_2025`  | `.R`   | Workflow for Figure 1: QQ-Plots and power analysis in main contents for Poisson-Gamma and ZINB |
| `refined_supp_simulations_10_13_2025`  | `.R`   | Supplementary tests and figure construction code (Figures S.1-S.4) |
---

#### Real Data Analysis
In the real data analysis section, we provide code to demonstrate the main real data contributions using the proposed SEF methodology as well as replicate the key figures with respect to these data.
This includes the regression pipeline and enrichment figures.

Within the `walkthrough` directory of the repository, we include the main function file, denoted as `main_github_10_14_2025.R`, which defines all the functions needed to conduct our analysis.
We apply these in practice in the walkthrough titled `SLE_SEF_regression_walkthrough_10_24_2025.R`, where we test for distributional differences between SLE and control populations with the CD8+ T-cell type.
This walkthrough also includes code snippets to reproduce enrichment analysis. Results can be replicated for the two other cell types. 

In `SLE_fig3_10_21_2025.R`, we provide code snippets that reproduce figures 3, S.5, S.6, and S.7,
which illustrate carrier densities, individual-level densities, and group-wise density comparisons for uniquely-identified genes (S100A4, TSC22D3, MT2A, and LTB.)

In the remaining files, we provide scripts run on the command line for the pseudobulk data and related comparisons. 


| File Name | Format     | Description |
|---------------|----------|------------|
| `SLE_SEF_regression_walkthrough_10_24_2025`      | `.R`  | Runs cell-type specific SEF regression modeling and testing; code to reproduce Figure 5 for a given cell-type (CD8+ T-cells as example) |
| ``  | `.R`   | Procedural pipeline for pseudobulk methods; run on command line but also can be run on any `R`-compatible environments |
| ``  | `.R`   | Workflow to reproduce Table 1 numbers and Figure 2 |
| ``  | `.R`   | Workflow to reproduce proportions in Table 2 |
| ``  | `.R`   | Workflow to reproduce Figure 3; Figures S.5, S.6, S.7 |
| `gen_sim_obj_EFFECT`  | `.R`   | generate simulated expression data based on effect size |
| `main_github_10_14_2025`  | `.R`   | `main` file defining functions used in real data analysis and downstream analysis |

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
