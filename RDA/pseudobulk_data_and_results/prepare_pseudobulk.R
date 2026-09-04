rm(list = ls())  
library(dplyr)
library(ggplot2)
library(Seurat)
library(edgeR)
library(DESeq2)
library(Matrix)
library(future)
library(MAST)
library(dplyr)
library(tibble)
source("/Users/zaqian/Desktop/merge_sim_RDA/merged_self_contained_RDA.R")
# notes: we don't use SCT, setting to null to save memory
# ------ get genes tested ------
sef_res_dir = file.path("/Users/zaqian/Desktop/merge_sim_RDA/sef_results")
flag_gene = T
if(flag_gene == F){
  p = 2
  cell_types = c("cd8","cd4","cM")
  for(i in 1:length(cell_types)){
    new_sef = readRDS(paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/with_revised_covariance_RDA_5_1_2026/res_",
                             cell_types[i], "_p_",p,"_sef_5_1_2026.RDS"))
    genes_tested = names(new_sef)
    genes_tested_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/with_revised_covariance_RDA_5_1_2026/pseudobulk/",
                               cell_types[i], "_genes_tested_5_4_2026.RDS") # based at least 300 cells, 3% non-zeroness
    saveRDS(genes_tested, genes_tested_name)
  }
}


####### load data ####### 
cell_type = "cM"
genes_tested = readRDS(paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/with_revised_covariance_RDA_5_1_2026/pseudobulk/",
                              cell_type,"_genes_tested_5_4_2026.RDS"))
donors_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/comparisons_4_20_2026/", cell_type,"_donors_used.RDS")
donors_tested = readRDS(donors_name)
if(cell_type == "cd4"){
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd4_sObj_new.RDS") 
  DefaultAssay(sObj) = "RNA"
  sObj[["SCT"]] = NULL
  sObj = subset(sObj, subset = donor_id != "ICC_control")
  sObj = subset(sObj, subset = disease_state %in% c("managed", "na"))  # keep managed SLE + healthy controls, drop flare + treated
  sObj = subset(sObj, subset = percent.mt < 10) # remove high % mitochondrial counts
  sObj = filter_extreme_coverage(sObj = sObj, log_scale = F)
  genes_to_remove = c(
    grep("^RPL",  rownames(sObj), value = TRUE),
    grep("^RPS",  rownames(sObj), value = TRUE),
    grep("^MRPL", rownames(sObj), value = TRUE),
    grep("^MRPS", rownames(sObj), value = TRUE),
    "MALAT1", "NEAT1", "XIST"
  )
  cat("Removing", length(genes_to_remove), "ribosomal genes and non-coding genes (MALAT1, NEAT1, XIST) \n")
  genes_keep = setdiff(rownames(sObj), genes_to_remove)
  sObj = subset(sObj, features = genes_keep)
  
} else if (cell_type == "cd8"){
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd8_sObj_new.RDS")
  DefaultAssay(sObj) = "RNA"
  sObj[["SCT"]] = NULL
  sObj = subset(sObj, subset = donor_id != "ICC_control")
  sObj = subset(sObj, subset = disease_state %in% c("managed", "na"))  # keep managed SLE + healthy controls, drop flare + treated
  sObj = subset(sObj, subset = percent.mt < 10) # remove high % mitochondrial counts
  sObj = filter_extreme_coverage(sObj = sObj, log_scale = F)
  genes_to_remove = c(
    grep("^RPL",  rownames(sObj), value = TRUE),
    grep("^RPS",  rownames(sObj), value = TRUE),
    grep("^MRPL", rownames(sObj), value = TRUE),
    grep("^MRPS", rownames(sObj), value = TRUE),
    "MALAT1", "NEAT1", "XIST"
  )
  cat("Removing", length(genes_to_remove), "ribosomal genes and non-coding genes (MALAT1, NEAT1, XIST) \n")
  genes_keep = setdiff(rownames(sObj), genes_to_remove)
  sObj = subset(sObj, features = genes_keep)
  
} else if (cell_type == "cM"){
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cM_sObj_new.RDS")
  DefaultAssay(sObj) = "RNA"
  sObj[["SCT"]] = NULL
  sObj = subset(sObj, subset = donor_id != "ICC_control")
  sObj = subset(sObj, subset = disease_state %in% c("managed", "na"))  # keep managed SLE + healthy controls, drop flare + treated
  sObj = subset(sObj, subset = percent.mt < 10) # remove high % mitochondrial counts
  sObj = filter_extreme_coverage(sObj = sObj, log_scale = F)
  genes_to_remove = c(
    grep("^RPL",  rownames(sObj), value = TRUE),
    grep("^RPS",  rownames(sObj), value = TRUE),
    grep("^MRPL", rownames(sObj), value = TRUE),
    grep("^MRPS", rownames(sObj), value = TRUE),
    "MALAT1", "NEAT1", "XIST"
  )
  cat("Removing", length(genes_to_remove), "ribosomal genes and non-coding genes (MALAT1, NEAT1, XIST) \n")
  genes_keep = setdiff(rownames(sObj), genes_to_remove)
  sObj = subset(sObj, features = genes_keep)
}

# ----- pseudobulk preparation -----
options(future.globals.maxSize = 8 * 1024^3)   # 8 GiB
pb = AggregateExpression(
  sObj,
  assays   = "RNA",
  group.by = "donor_id",
  return.seurat = F,
  slot     = "counts"
)$RNA
colnames(pb) = sub("^g", "", colnames(pb))
# colnames(pb) <- gsub("-", "_", sub("^g(?=[0-9])", "", colnames(pb), perl = TRUE))

md = sObj@meta.data
md = data.frame(unique(md[, c("donor_id", "disease")]))
rownames(md) = md$donor_id
md = md[colnames(pb), , drop = FALSE]
setdiff(rownames(md), colnames(pb))
stopifnot(all(rownames(md) == colnames(pb)))
stopifnot(!anyNA(md$disease))
any(duplicated(md$donor_id))   # must be FALSE
md$disease = factor(md$disease, levels = c("normal", "systemic lupus erythematosus"))   # adjust to match your data
pb = pb[, colnames(pb) %in% donors_tested]
md = md[rownames(md) %in% donors_tested, , drop = FALSE]

dim(pb)
dim(md)
# NOTE: we don't need to run filter_SLE() because we already know the donors that have at least 300 cells, as they are removed from the data. 
# we also don't even test genes with < 3% non-zero expression

# ----- save RDS files -----
out.dir = file.path("/Users/zaqian/Desktop/finalSims/submit_to_github/application/pseudobulk_data_and_results/pb_dat")
saveRDS(pb, file.path(out.dir,paste0("pb_",cell_type,"_matrix.RDS")))
saveRDS(md, file.path(out.dir,paste0("pb_",cell_type,"_metadata.RDS")))



