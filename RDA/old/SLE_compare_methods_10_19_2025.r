# comparison with pseudobulk method, run on command line
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pbmcapply)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 
covariate_df = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript fair_pb.R <cell_type>\n  <cell_type> must be one of: cM, cd8, cd4", call. = FALSE)
}
cell_type = tolower(args[[1]])
valid_ct = c("cM", "cd8", "cd4") 
if (cell_type == "cm"){ # re-correct back to cM format
  cell_type = "cM"
}
if (!cell_type %in% c("cM", "cd8", "cd4")) {
  stop("cell_type must be one of: cM, cd8, cd4", call. = FALSE)
}

if (cell_type == "cd8") {
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_genes_processed.RDS") #renamed to cd8 for consistency, original was "updated_cd8"
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  pb_sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/pb_cd8_sObj.RDS")
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd8_fixed_pvals.csv", row.names = 1)
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
  
  
}else if(cell_type == "cd4"){
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_genes_processed.RDS") #renamed to cd4 for consistency, original was "TEST_cd4"
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  pb_sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/pb_cd4_sObj.RDS")
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd4_fixed_pvals.csv", row.names = 1)
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
  
}else if(cell_type == "cM"){
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_genes_processed.RDS") #renamed to cM for consistency, original was "monocytes"
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  pb_sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/pb_cM_sObj.RDS")
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cM_fixed_pvals.csv", row.names = 1) # renamed to just "sef_cM", original was "filter_before_" 
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
}

# pseudobulk
fair_pseudo_wilcoxon = function(pb_sObj, donors_i, gene_i){
  # function: runs pseudobulked wilcoxon rank-sum test for each gene using every donor used in SEF regression test
  # input: pseudobulked seurat object; donor IDs; gene of interest
  # output: avg log FC and p-value
  single_gene_sObj = subset(pb_sObj,
                            subset  = donor_id %in% donors_i,
                            features = gene_i)
  design_matrix = cbind(single_gene_sObj@meta.data[,c("donor_id","disease")],
                        expression = as.vector(single_gene_sObj@assays$SCT@layers$data)) #use on pseudobulk values
  
  x = design_matrix$expression[design_matrix$disease == "normal"]
  y = design_matrix$expression[design_matrix$disease == "systemic lupus erythematosus"]
  wilcox_res = wilcox.test(x = x, y = y,alternative = "two.sided")
  logFC = mean(log2(y + 1)) - mean(log2(x + 1)) 
  return(
    data.frame(
      gene = gene_i,
      p_value = wilcox_res$p.value,
      statistic = unname(wilcox_res$statistic),
      logFC = logFC, # compute average logFC
      n_normal = sum(design_matrix$disease == "normal"),
      n_sle = sum(design_matrix$disease == "systemic lupus erythematosus"),
      stringsAsFactors = FALSE
    )
  )
}

results_list = pbmclapply( # run in parallel
  genes_of_interest,
  mc.cores = 3,
  function(g) {
    donors_i = donors_to_use_per_gene[[g]]
    try(
      fair_pseudo_wilcoxon(pb_sObj, donors_i, g),
      silent = TRUE
    )
  }
)
results_df = do.call(rbind, Filter(Negate(is.null), results_list))
results_df$p_val_adj = p.adjust(results_df$p_value, method = "bonferroni")
revised_pseudobulk = results_df %>%
  arrange(p_val_adj)
print(dim(subset(revised_pseudobulk, p_val_adj < 0.05)))

# saved as pb_CELLTYPE_10_15_2025.csv
