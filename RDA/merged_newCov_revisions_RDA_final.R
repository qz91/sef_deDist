# Wider pool of genes tested, more donors, criteria dropped to nonzero counts total of at least 5%?
rm(list = ls())  
library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
source_dir = file.path("")
source(file.path(source_dir, "merged_self_contained_RDA_final.R"))

#----- cell type ----- 
cell_type = "cd8"

## ----- otherwise, load data -----
if(cell_type == "cd4"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd4_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
  
} else if (cell_type == "cd8"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd8_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
  
} else if (cell_type == "cM"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cM_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
}
cell_counts = sObj_meta %>%
  dplyr::count(donor_id, disease, disease_state, ind_cov, name = "n_cells") %>%
  dplyr::arrange(desc(n_cells))
cell_counts = cell_counts[order(cell_counts$n_cells, decreasing = F), ]
counts_dir = file.path(source_dir, "counts_by_donor")
# dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
# write.csv(cell_counts, file.path(counts_dir, paste0(cell_type,"_cell_counts_by_donor.csv")))
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

lcf_filtered = filter_SLE(lcf_matrix, donor_ids, min_cells = 300, min_prop = 0.03)
genes_of_interest = lcf_filtered$genes
donors_keep = lcf_filtered$donors
# donors_name = paste0(cell_type,"_donors_used.RDS")
# saveRDS(donors_keep, donors_name)
sObj_meta = sObj_meta[sObj_meta$donor_id %in% lcf_filtered$donors, ]
lcf2= lcf_filtered$matrix
dim(lcf2)
rm(donor_ids) # save memory
rm(lcf_matrix)
rm(lcf_filtered)
gc()

p = 2
bw = 0.5
options(future.globals.maxSize = 5 * 1024^3)  # 3 GiB
plan(multisession, workers = 2)
handlers(global = TRUE)  # enables progress bars

new_sef = with_progress({
  p_bar = progressor(along = genes_of_interest)
  future_lapply(genes_of_interest, function(name) {
    p_bar(name)
    tryCatch({
      exprVec = lcf2[name, , drop = FALSE]
      run_sef_procedure(
        exprVec = exprVec,
        sObj_meta = sObj_meta,
        p = p,
        bandwidth = bw,
        plot_flag = FALSE
      )
    }, error = function(e) {
      message(paste0("FAILED: ", name, " - ", e$message))
      return(list(gene = name, error = e$message))
    })
  }, future.seed = TRUE)
})
new_sef = setNames(new_sef, genes_of_interest)
sef_dir = file.path(source_dir, "sef_results")
dir.create(sef_dir, recursive = TRUE, showWarnings = FALSE)

date_tag = format(Sys.Date(), "%m_%d_%Y")
res_name = file.path(sef_dir, paste0("res_", cell_type,"_p_",p,"_sef_",date_tag,".RDS"))
saveRDS(new_sef, res_name)

pvals = unlist(sapply(new_sef, function(x) x$pval_1))
pval_df = data.frame(p_value = pvals, row.names = names(pvals))
pval_df$bonferroni_pval = p.adjust(pval_df$p_value, method = "bonferroni")
pval_df$fdr =  p.adjust(pval_df$p_value, method = "fdr")
dim(subset(pval_df, bonferroni_pval < 0.05)) 
dim(subset(pval_df, fdr < 0.05)) 
pval_df = pval_df %>% arrange(bonferroni_pval)
sig = subset(pval_df, bonferroni_pval < 0.05)

pval_dir = file.path(source_dir, "")
dir.create(pval_dir, recursive = TRUE, showWarnings = FALSE)

date_tag = format(Sys.Date(), "%m_%d_%Y")
write.csv(sig, file.path(pval_dir, paste0(cell_type,"_p_", p, "_sef_pvals_",date_tag,".csv")))


# ----- enrichment analysis -----
cell_type = "cd8" # options include cd8, cd4, cM
p = 2
sig_genes = rownames(sig)

# enrichment
orgDb = org.Hs.eg.db; ontology = "BP"
gene_df = bitr(sig_genes, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)
ego = enrichGO(gene = gene_df$ENTREZID, OrgDb = orgDb,
               keyType = "ENTREZID", ont = ontology,
               pAdjustMethod = "BH", readable = TRUE)
ego = clusterProfiler::simplify(ego, by = "p.adjust")
viewer = c("ID","Description", "RichFactor","FoldEnrichment","p.adjust","qvalue", "geneID","Count")
sig_ego = ego@result[ego@result$p.adjust < 0.05,viewer]
ego_name = paste0(cell_type,"_enrichment_", date_tag,".csv")
write.csv(sig_ego, file.path(enrich_dir, ego_name))


# ----- cd8 unique SEF genes higher p enrichment -----
cell_type = "cd8"
pb_name = paste0("deseq_",cell_type,"_sig_pvals_4_20_2026.csv")
pb_deseq = read.csv(pb_name, row.names = 1)
pb_genes = rownames(pb_deseq)

unique_genes = setdiff(sig_genes, pb_genes)

# enrichment
orgDb = org.Hs.eg.db; ontology = "BP"
gene_df = bitr(unique_genes, 
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
ego = enrichGO(gene = gene_df$ENTREZID, OrgDb = orgDb,
               keyType = "ENTREZID", ont = ontology,
               pAdjustMethod = "BH", readable = TRUE)

viewer = c("ID","Description", "RichFactor","FoldEnrichment","p.adjust","qvalue", "geneID","Count")
sig_ego = ego@result[ego@result$p.adjust < 0.1,viewer]
date_tag = format(Sys.Date(), "%m_%d_%Y")
ego_name = paste0(cell_type,"_p", p,"_sef_unique_enrichment_", date_tag,".csv")
write.csv(sig_ego, file.path(enrich_dir, ego_name))


pathways_to_show = c(
  "granzyme-mediated programmed cell death signaling pathway",
  "MHC class II protein complex assembly",
  "peptide antigen assembly with MHC class II protein complex",
  "negative regulation of lymphocyte migration",
  "negative regulation of T-helper cell differentiation",
  "MHC protein complex assembly",
  "antigen processing and presentation of endogenous peptide antigen",
  "negative regulation of CD4-positive, alpha-beta T cell differentiation",
  "regulation of T-helper cell differentiation",
  "regulation of CD4-positive, alpha-beta T cell differentiation",
  "regulation of alpha-beta T cell differentiation",
  "regulation of CD4-positive, alpha-beta T cell activation",
  "T-helper cell differentiation",
  "CD4-positive, alpha-beta T cell differentiation involved in immune response",
  "alpha-beta T cell activation involved in immune response",
  "alpha-beta T cell differentiation involved in immune response",
  "regulation of leukocyte cell-cell adhesion",
  "positive regulation of leukocyte cell-cell adhesion",
  "leukocyte cell-cell adhesion",
  "regulation of T cell activation"
)

enrich_fdr_10 = sig_ego %>%
  filter(p.adjust < 0.10, Description %in% pathways_to_show) %>%
  mutate(
    neg_log_adj_p = -log10(p.adjust),
    Description = factor(Description, levels = rev(pathways_to_show))
  )

p = enrich_fdr_10 %>%
  ggplot(aes(x = FoldEnrichment, y = Description, fill = neg_log_adj_p)) +
  geom_col() +
  scale_fill_viridis_c(
    option = "plasma",
    name = expression(-log[10]~adj.~p)
  ) +
  ggtitle("SEF-Unique Enriched Pathways (CD8+)") +
  xlab("Fold Enrichment") +
  ylab(NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p


# ----- argument: unique PB are useless -----
orgDb = org.Hs.eg.db; ontology = "BP"
gene_df = bitr(setdiff(pb_genes, sig_genes), 
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
ego = enrichGO(gene = gene_df$ENTREZID, OrgDb = orgDb,
               keyType = "ENTREZID", ont = ontology,
               pAdjustMethod = "BH", readable = TRUE)

View(ego@result)
viewer = c("ID","Description", "RichFactor","FoldEnrichment","p.adjust","qvalue", "geneID","Count")
pb_ego = ego@result[ego@result$p.adjust < 0.1,viewer]
ego_name = paste0("pb_vs_p",p,"_",cell_type,"_enrichment_", date_tag,".csv")
write.csv(pb_ego, file.path(enrich_dir, ego_name))

# ----- table counts for DD genes between p = 2, p = 3
sef_dir = file.path(source_dir, "sef_results")
pval_dir = file.path(source_dir, "pvals")
p = 3
cell_types = c("cd8", "cd4", "cM")
intersect_results = data.frame()
for (cell_type in cell_types) {
  sig = read.csv(file.path(pval_dir, paste0(cell_type, "_p_", p, "_sef_pvals_07_01_2026.csv")), row.names = 1)
  p2_sig = read.csv(file.path(paste0(cell_type, "_sef_pvals_5_1_2026.csv")), row.names = 1)
  sig_genes = rownames(sig)
  p2_genes = rownames(p2_sig)
  
  intersect_results = rbind(
    intersect_results,
    data.frame(
      cell_type = cell_type,
      n_p3 = length(sig_genes),
      n_p2 = length(p2_genes),
      pct_intersect_with_p2 = length(intersect(sig_genes, p2_genes)) / length(p2_genes)
    )
  )
}
intersect_results

overlap_results = data.frame()

p_choice = 4
overlap_results = data.frame()
for (cell_type in cell_types) {
  p2_sig = read.csv(file.path(paste0(cell_type, "_sef_pvals_5_1_2026.csv")), row.names = 1)
  
  p_choice_sig = read.csv(file.path(pval_dir, paste0(cell_type, "_p_", p_choice, "_sef_pvals_07_01_2026.csv")), row.names = 1)
  
  p2_genes = unique(rownames(p2_sig))
  p_choice_genes = unique(rownames(p_choice_sig))
  
  n_overlap = length(intersect(p_choice_genes, p2_genes))
  overlap_coefficient = n_overlap / min(length(p_choice_genes), length(p2_genes))
  
  overlap_results = rbind(
    overlap_results,
    data.frame(
      cell_type = cell_type,
      p_choice = p_choice,
      n_p_choice = length(p_choice_genes),
      n_p2 = length(p2_genes),
      n_overlap = n_overlap,
      overlap_coefficient = overlap_coefficient
    )
  )
}
overlap_results

# ----- get counts per donor -----
cell_types = c("cd8", "cd4", "cM")
for (ct in cell_types) {
  cts_donors = read.csv(file.path(source_dir,"counts_by_donor",paste0(ct, "_cell_counts_by_donor.csv")),row.names = 1)
  
  eligible_donors = subset(
    cts_donors,
    donor_id != "ICC_control" &
      n_cells >= 300
  )
  
  cat("\n", ct, ":\n", sep = "")
  print(table(eligible_donors$disease))
  write.csv(eligible_donors, file.path(source_dir,"counts_by_donor" ,paste0(ct, "_eligible_counts_by_donor_cell_counts_by_donor.csv")))
}


