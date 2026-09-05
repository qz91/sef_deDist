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
source("merged_self_contained_RDA_final.R")
# notes: we don't use SCT, setting to null to save memory

####### load data ####### 
cell_type = "cM"
genes_tested = readRDS(paste0(cell_type,"_genes_tested_5_4_2026.RDS"))
donors_name = paste0(cell_type,"_donors_used.RDS")
donors_tested = readRDS(donors_name)
if(cell_type == "cd4"){
  sObj = readRDS("cd4_sObj_new.RDS") 
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
  sObj = readRDS("cd8_sObj_new.RDS")
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
  sObj = readRDS("cM_sObj_new.RDS")
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
out.dir = file.path("pseudobulk_data_and_results")
saveRDS(pb, file.path(out.dir,paste0("pb_",cell_type,"_matrix.RDS")))
saveRDS(md, file.path(out.dir,paste0("pb_",cell_type,"_metadata.RDS")))



