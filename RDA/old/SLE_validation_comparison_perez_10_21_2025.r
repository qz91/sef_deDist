# validation with perez et al supplementary genes
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 

perez = read.csv("./density_estimation/sef_lupus_results/exploration/perez_validation_DEGs.csv", row.names = 1)
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")
donors_cd8 = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_genes_processed.RDS") 
donors_cd4 = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_genes_processed.RDS") 
donors_cM = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_genes_processed.RDS") 
cd8_validated_genes = intersect(rownames(perez)[perez$T8 == 1], names(donors_cd8))
cd4_validated_genes= intersect(rownames(perez)[perez$T4 == 1], names(donors_cd4))
cM_validated_genes = intersect(rownames(perez)[perez$cM == 1], names(donors_cM))

# load results
sef_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cM_fixed_pvals.csv", row.names = 1)
sef_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd8_fixed_pvals.csv", row.names = 1)
sef_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd4_fixed_pvals.csv", row.names = 1)
SCT_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/pb_cM_10_15_2025.csv", row.names = 1)
SCT_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/pb_cd8_10_15_2025.csv", row.names = 1)
SCT_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/pb_cd4_10_15_2025.csv", row.names = 1)

# cd8
sum(row.names(sef_cd8) %in% cd8_validated_genes)/length(cd8_validated_genes)
sum(SCT_cd8$gene %in% cd8_validated_genes)/length(cd8_validated_genes)

#cd4
sum(row.names(sef_cd4) %in% cd4_validated_genes)/length(cd4_validated_genes) 
sum(SCT_cd4$gene %in% cd4_validated_genes)/length(cd4_validated_genes)

# cM
sum(row.names(sef_cM) %in% cM_validated_genes)/length(cM_validated_genes)
sum(SCT_cM$gene %in% cM_validated_genes)/length(cM_validated_genes)

