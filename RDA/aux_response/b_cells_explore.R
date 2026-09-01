# b-cells 
rm(list = ls())  
library(dplyr)
library(ggplot2)
library(Seurat)
source_dir = file.path("/Users/zaqian/Desktop/merge_sim_RDA")
source(file.path(source_dir, "merged_self_contained_RDA.R"))
cell_type = "b_cells"

subset_flag = T
if(subset_flag == F){
  og.data.dir = file.path("/Users/zaqian/Desktop/lupus")
  sObj = readRDS(file.path(og.data.dir, "lupus_seurat.rds"))
  b_cells = subset(sObj, cell_type == "B cell")
  b_cells[["percent.mt"]] = PercentageFeatureSet(
    b_cells,
    pattern = "^MT-"
  )
  saveRDS(b_cells, file.path(source_dir, "b_cells","b_sObj.RDS"))
}


filtered_flag = F
if(filtered_flag == F){
  sObj = readRDS(file.path(source_dir, "b_cells","b_sObj.RDS")) 
  DefaultAssay(sObj) = "RNA"
  print(dim(sObj)) # 20476 x 149840 (G x C)
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
  
  rna_matrix = LayerData(sObj, assay = "RNA", layer = "counts")
  lcf_matrix = lib_size_correction_sparse(rna_matrix)
  # file.path(source_dir, "b_cells","b_sObj.RDS")
  saveRDS(lcf_matrix,file.path(source_dir, "b_cells","b_lcf_matrix.RDS")) #file.path(source_dir, "b_cells","b_lcf_matrix.RDS")
  meta_name = file.path(source_dir,"b_cells",paste0(cell_type, "_new_metadata.RDS"))
  # paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  saveRDS(sObj@meta.data, file = meta_name)
}

lcf_matrix = readRDS(file.path(source_dir, "b_cells","b_lcf_matrix.RDS"))
meta_name = file.path(source_dir,"b_cells",paste0(cell_type, "_new_metadata.RDS"))
sObj_meta = readRDS(meta_name)

cell_counts = sObj_meta %>%
  dplyr::count(donor_id, disease, disease_state, ind_cov, name = "n_cells") %>%
  dplyr::arrange(desc(n_cells))
cell_counts = cell_counts[order(cell_counts$n_cells, decreasing = F), ]
write.csv(cell_counts, file.path(source_dir,"b_cells", paste0(cell_type,"_cell_counts_by_donor.csv")))

head(cell_counts,10)

dup_donors = names(which(table(cell_counts$ind_cov) > 1))
stopifnot(
  "Duplicate ind_cov values found across donors; please remove flare samples." = length(dup_donors) == 0
)

donor_ids = sObj_meta$donor_id
length(donor_ids)
dim(lcf_matrix)
stopifnot(
  "Column names of library-corrected count matrix do not match rownames of sObj_meta" = 
    all(colnames(lcf_matrix) == rownames(sObj_meta)),
  "Length of donor_ids does not match number of columns in library-corrected count matrix" = 
    length(donor_ids) == ncol(lcf_matrix)
)


# Within each cell-type group, we filter out donors with fewer than 300 cells, 
# and filter out genes with donor-level cell-type-specific nonzero expression rates below 3%. 
# note: at 2%, cM: 2665, cd4: 1687, cd8: 1516 total gene pools
lcf_filtered = filter_SLE(lcf_matrix, donor_ids, min_cells = 300, min_prop = 0.03)
genes_of_interest = lcf_filtered$genes
donors_keep = lcf_filtered$donors
donors_name = file.path(source_dir,"b_cells",paste0(cell_type,"_donors_used.RDS"))
saveRDS(donors_keep, donors_name)
sObj_meta = sObj_meta[sObj_meta$donor_id %in% lcf_filtered$donors, ]
lcf2= lcf_filtered$matrix
dim(lcf2)
rm(donor_ids) # save memory
rm(lcf_matrix)
rm(lcf_filtered)
gc()
saveRDS(lcf2, file.path(source_dir,"b_cells",""))

# ----- final summary counts -----
print(length(genes_of_interest))

final_disease_summary <- sObj_meta %>%
  group_by(disease) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    n_cells = n(),
    .groups = "drop"
  )

final_disease_summary
