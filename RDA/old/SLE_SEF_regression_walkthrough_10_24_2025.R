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
cell_type = "cd8" # provide workflow with CD8+ T-cells
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")
sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_metadata.RDS") 
exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_exprMat.RDS") 
donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_genes_processed.RDS") 
genes_of_interest = names(donors_to_use_per_gene)
print(length(genes_of_interest))

####### sef regression ####### 
p = 2
res_list = setNames(pbmclapply(genes_of_interest, function(name) { # run in parallel
  tryCatch({
    exprMat1 = exprMatReal[name, , drop = F]
    stablerunSeuratCounts_Concise(exprMat = exprMat1, donor_list = donors_to_use_per_gene[[name]], sObj_meta = sObj_metadata, p = p, plot_flag = F) 
  }, error = function(e) {
    message(paste0("Error in processing: ", name, " - ", e$message))
    return(NULL)
  })
}, mc.cores = 3, ignore.interactive = TRUE), genes_of_interest)

# extract p-values, adjust for multiple comparisons
res_list = Filter(Negate(is.null), res_list)
pvals = unlist(sapply(res_list, function(x) x$pval_1))
pval_df = data.frame(p_value = pvals, row.names = names(pvals))
pval_df$bonferroni_pval = p.adjust(pval_df$p_value, method = "bonferroni")
dim(subset(pval_df, bonferroni_pval < 0.05)) 
pval_df = pval_df %>% arrange(bonferroni_pval)
sig = subset(pval_df, bonferroni_pval < 0.05)

# enrichment Figure 5
res = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/results/enrichment_res/cd8_enriched_pathways.csv", row.names = 1)
res$Count = sapply(strsplit(res$geneID, "/"), length)
ego_reconstructed = new("enrichResult")
ego_reconstructed@result = res
ego_reconstructed@organism = "Homo sapiens"
ego_reconstructed@keytype = "ENTREZID"
ego_reconstructed@ontology = "BP"
ego_reconstructed@readable = TRUE
barplot(height = ego_reconstructed) + ggtitle("CD8+") + guides(size = "none")

