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
# cell_type = "cM"
# cell_type = "cd8"
# cell_type = "cd4"
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")

if(cell_type == "cd4"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_metadata.RDS") 
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_exprMat.RDS") 
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_genes_processed.RDS")
} else if (cell_type == "cd8"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_metadata.RDS") 
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_exprMat.RDS") 
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_genes_processed.RDS") 
} else if (cell_type == "cM"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_exprMatReal.RDS")
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_genes_processed.RDS") 
}
# specific_genes = c("S100A9", "ISG15","CYBA", "S100A4", "S100A8") # cM; S.7: S100A4
# specific_genes = c("HLA-A","NKG7", "LGALS1", "MT2A", "LTB") # cd8; S.7: MT2A, LTB
# specific_genes = c("B2M", "IL32","LTB", "TSC22D3", "KLF6") # cd4; S.7: TSC22D3

don = "1101"
genes_of_interest = specific_genes
p = 2
res_list = setNames(pbmclapply(genes_of_interest, function(name) { #run sef regression for the genes of interest
  tryCatch({
    exprMat1 = exprMatReal[name, , drop = F]
    stablerunSeuratCounts_Concise(exprMat = exprMat1, donor_list = donors_to_use_per_gene[[name]], sObj_meta = sObj_metadata, p = p, plot_flag = F) 
  }, error = function(e) {
    message(paste0("Error in processing: ", name, " - ", e$message))
    return(NULL)
  })
}, mc.cores = 3, ignore.interactive = TRUE), genes_of_interest)


for (g in specific_genes) {
  print(g)
  carrier_test = plot_carrier_with_hist(exprMat = exprMatReal, gene = g, sObj_meta = sObj_metadata, donors_to_use_list = donors_to_use_per_gene)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-carrier-",g, ".pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(carrier_test)
  dev.off()
  
  comparison_test = group_model_comparison_plot(sef_object = res_list, gene = g)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-model-comparison-",g, ".pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(comparison_test)
  dev.off()
  
  ind_test = get_individual_sef_counts_model(donor_id = don, gene = g, exprMat = exprMatReal, sObj_meta = sObj_metadata, donor_list = donors_to_use_per_gene)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-",g,"-donor-",don,"-ind-model.pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(ind_test)
  dev.off()
}

