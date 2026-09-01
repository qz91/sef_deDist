####### load relevant packages ####### 
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(parallel)
library(doSNOW)
library(foreach) 
library(doParallel)
library(pbmcapply)  
library(Ake)
source("./density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 

####### load data ####### 
cat("[1] loading covariate data...\n")
covariate_df = read.csv("covariate_data.csv")
cat("[1] covariate data instantiated\n")

# Parse command line arguments
args = commandArgs(trailingOnly = TRUE)
n_perm = if (length(args) >= 1) suppressWarnings(as.integer(args[[1]])) else 10
if (is.na(n_perm) || n_perm <= 0) {
  cat("Invalid or missing permutation count. Using default n_perm = 10\n")
  n_perm = 10
}
choice = if (length(args) >= 2) suppressWarnings(as.integer(args[[2]])) else 1
if (!choice %in% c(1, 2, 3)) {
  stop("Invalid cell type choice. Must be 1 (Monocytes), 2 (CD4), or 3 (CD8).")
}

if (choice == 1) {
  cat("You chose Monocytes\n")
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_exprMatReal.RDS")
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_genes_processed.RDS")
} else if (choice == 2) {
  cat("You chose CD4\n")
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_exprMat.RDS")
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_genes_processed.RDS")
} else if (choice == 3) {
  cat("You chose CD8\n")
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_exprMat.RDS")
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_genes_processed.RDS")
} else {
  cat("Invalid choice\n")
}

genes_of_interest = names(donors_to_use_per_gene)
print(length(genes_of_interest))
print("loaded up data.")


####### permutation test ####### 
base_seed = 100
all_perm_results = vector("list", n_perm)
tested_genes = names(donors_to_use_per_gene)


# ---- n_perm loops ---- 
for (i in seq_len(n_perm)) { # run sequentially for each permutation iteration
  seed = base_seed + i
  message(paste0("Running permutation %d with seed %d", i, seed))
  
  # run in parallel for each gene
  perm_result = setNames(pbmclapply(tested_genes, function(name) {
    tryCatch({
      exprMat1 = exprMatReal[name, , drop = F]
      stablePermutationTest(
        exprMat = exprMat1,
        covariate_df = covariate_df,
        sObj_meta = sObj_metadata,
        donor_list = donors_to_use_per_gene[[name]],
        seed = seed, p = 2, plot_flag = F
      )
    }, error = function(e) {
      message(paste0("Error in processing: ", name, " - ", e$message))
      return(NULL)
    })
  }, mc.cores = 5, ignore.interactive = TRUE), tested_genes)
  all_perm_results[[i]] = perm_result
  gc() #save memory
}
cat("done\n")

if(choice == 1){
  cat("saving monocyte permutation test results...\n")
  saveRDS(all_perm_results, file = "./cM_perm_res_new.rds")
  cat("saved\n")
} else if(choice == 2){
  cat("saving CD4+ permutation test results...\n")
  saveRDS(all_perm_results, file = "./cd4_perm_res_new.rds")
  cat("saved\n")
} else if(choice == 3){
  cat("saving CD8+ permutation test results...\n")
  saveRDS(all_perm_results, file = "./cd8_perm_res_new.rds")
  cat("saved\n")
}
cat(paste0("permutation test with ", n_perm, " total tests complete.\n"))

# ---- make qqplots ---- 
pval_vector = unlist( # extract from permutation object, from each iteration and for each gene the p-value
  lapply(all_perm_results, function(inner_list) {
    vapply(inner_list, function(res) {
      if (!is.null(res) && "pval_1" %in% names(res)) res$pval_1 else NA_real_
    }, numeric(1))
  }),
  recursive = T,
  use.names = F)

print(mean(pval_vector < 0.05)) # Proportion of false positives

permutation_plot = gg_qqplot(pval_vector) + ggtitle("QQ Plot: Permutation Test") +
  theme(plot.title = element_text(size = 18))
GraphName = paste0("permutation_plot_new.pdf")
pdf(GraphName, width = 5, height = 5)
print(permutation_plot)
dev.off()
