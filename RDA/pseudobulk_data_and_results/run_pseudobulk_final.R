# run_pseudobulk_final.R runs pseudobulk tests
rm(list = ls())  
library(dplyr)
library(Seurat)
library(edgeR)
library(DESeq2)
library(Matrix)
library(future)
library(tibble)
cell_type = "cd8" # options: cd8, cd4, cM
out.dir = file.path("pb_dat")
pb = readRDS(file.path(out.dir,paste0("pb_",cell_type,"_matrix.RDS")))
md = readRDS(file.path(out.dir,paste0("pb_",cell_type,"_metadata.RDS")))

genes_tested = readRDS(paste0(cell_type,"_genes_tested_5_4_2026.RDS")) # stored in /RDA/pseudobulk_data_and_results/

# ----- edgeR -----
design = model.matrix(~ disease, data = md)
dge  = DGEList(counts = pb, samples = md)
keep = filterByExpr(dge, design = design)      # filter on full set
dge  = dge[keep, , keep.lib.sizes = FALSE]
dge  = calcNormFactors(dge, method = "TMM")    # TMM on full filtered set
dge = dge[rownames(dge) %in% genes_tested, , keep.lib.sizes = TRUE]

dge = estimateDisp(dge, design)
fit = glmQLFit(dge, design)
qlf = glmQLFTest(fit, coef = 2)

res_edgeR = topTags(qlf, n = Inf, adjust.method = "none")$table
res_edgeR$bonferroni = p.adjust(res_edgeR$PValue, method = "bonferroni")
sum(res_edgeR$bonferroni < 0.05)

edgeR_name1 = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/with_revised_covariance_RDA_5_1_2026/pseudobulk/all_pvals/edgeR_",
                     cell_type,"_all_pvals_5_5_2026.csv")
write.csv(res_edgeR, edgeR_name1)
sig_edgeR = subset(res_edgeR, bonferroni < 0.05)


# ----- DESeq2 -----
dds = DESeqDataSetFromMatrix(
  countData = pb,
  colData   = md,
  design    = ~ disease
)
# dds <- dds[rownames(dds) %in% genes_tested, ]
dds = estimateSizeFactors(dds)
# Now restrict to genes_tested; size factors stay attached to the object
dds = dds[rownames(dds) %in% genes_tested, ]
# Continue with dispersion + Wald on the restricted set
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)

res_deseq = results(dds,
                    name = resultsNames(dds)[2],
                    independentFiltering = FALSE,
                    cooksCutoff = FALSE)
res_deseq = as.data.frame(res_deseq)

res_deseq$bonferroni = p.adjust(res_deseq$pvalue, method = "bonferroni")
cat("DESeq2:", nrow(res_deseq), "genes tested;",
    sum(res_deseq$bonferroni < 0.05, na.rm = TRUE), "significant at Bonferroni < 0.05\n")
sig_deseq = subset(res_deseq, bonferroni < 0.05)


# ----- Wilcoxon -----
# follows: https://rpubs.com/LiYumei/806213
conditions = data.frame(conditions = md$disease)
rownames(conditions) = rownames(md)
y = edgeR::DGEList(counts = pb, group = md$disease)
keep = edgeR::filterByExpr(y)
y = y[keep, , keep.lib.sizes = FALSE]
y = edgeR::calcNormFactors(y, method = "TMM")
count_norm = edgeR::cpm(y)
count_norm = as.data.frame(count_norm)
count_norm = count_norm[rownames(count_norm) %in% genes_tested, , drop = FALSE] # wilcoxon on genes tested
g1 = colnames(count_norm)[md$disease == "normal"] 
g2 = colnames(count_norm)[md$disease == "systemic lupus erythematosus"]

pvalues = sapply(seq_len(nrow(count_norm)), function(i) {
  data = cbind.data.frame(
    gene = as.numeric(count_norm[i, ]),
    conditions
  )
  if (length(unique(data$gene)) <= 1) return(1)
  suppressWarnings(wilcox.test(gene ~ conditions, data = data, exact = FALSE)$p.value)
})
wilcox_res = data.frame(
  gene = rownames(count_norm),
  logFC = log2(rowMeans(count_norm[, g2, drop = FALSE])) - # log fold change of average pseudobulked expressions per gene
    log2(rowMeans(count_norm[, g1, drop = FALSE])),
  p_val = pvalues,
  bonferroni = p.adjust(pvalues, method = "bonferroni"),
  padj = p.adjust(pvalues, method = "fdr"),
  row.names = NULL,
  stringsAsFactors = FALSE
)

head(wilcox_res)
sum(wilcox_res$bonferroni < 0.05)

# ----- pseudobulk-specific -----
p = 2
cell_type = "cd8"
if(p == 2){
  df_name = paste0(cell_type,"_sef_pvals_5_1_2026.csv")
} else if (p == 3){
  df_name = paste0(cell_type,"_p_",p,"_sef_pvals_5_1_2026.csv")
}
sig = read.csv(df_name, row.names = 1)

deseq_name = paste0("deseq_",cell_type,"_sig_pvals_4_20_2026.csv")
sig_deseq = read.csv(deseq_name, row.names = 1)
edgeR_name = paste0("edgeR_",cell_type,"_sig_pvals_4_20_2026.csv")
sig_edgeR = read.csv(edgeR_name, row.names = 1)
wilcox_name = paste0("wilcox_",cell_type,"_sig_pvals_5_4_2026.csv")
sig_wilcox = read.csv(wilcox_name, row.names = 1)

c(dim(sig_edgeR), dim(sig_deseq), dim(sig_wilcox), dim(sig))

setdiff(row.names(sig), row.names(sig_edgeR))
length(setdiff(row.names(sig), row.names(sig_edgeR)))
length(intersect(row.names(sig), row.names(sig_edgeR)))
length(setdiff(row.names(sig_edgeR), row.names(sig)))


setdiff(row.names(sig), row.names(sig_deseq)) 
length(intersect(row.names(sig), row.names(sig_deseq)))
length(setdiff(row.names(sig), row.names(sig_deseq)))
length(setdiff(row.names(sig_deseq), row.names(sig)))
setdiff(row.names(sig_deseq), row.names(sig))

setdiff(row.names(sig), sig_wilcox$gene) 
setdiff(sig_wilcox$gene, row.names(sig))
length(intersect(row.names(sig), sig_wilcox$gene))
length(setdiff(row.names(sig), sig_wilcox$gene))
length(setdiff(sig_wilcox$gene, row.names(sig)))

u_edgeR = setdiff(row.names(sig_edgeR), row.names(sig))
unique_edgeR = sig_edgeR[u_edgeR,]
u_deseq = setdiff(row.names(sig_deseq), row.names(sig))
unique_deseq = sig_deseq[u_deseq,]
unique_deseq$logFC = unique_deseq$log2FoldChange
u_wilcox = setdiff(sig_wilcox$gene, row.names(sig))
unique_wilcox = subset(sig_wilcox, gene %in% u_wilcox)

# ----- make volcano plots -----
get_volcano_plot = function(df, pval_col = "bonferroni", fc_thresh = 0.5, pval_thresh = 0.05, p_title = "Volcano Plot") {
  df$neg_log10_pval = -log10(df[,pval_col])
  ggplot(df, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(color = "cornflowerblue", alpha = 0.8, size = 2) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      title = p_title
    ) +
    theme(plot.title = element_text(hjust = 0.5))
}

plt1 = get_volcano_plot(df = unique_edgeR, pval_col = "bonferroni", p_title = "CD8+ edgeR-specific")
plt2 = get_volcano_plot(df = unique_deseq, pval_col = "bonferroni", p_title = "CD8+ DESeq2-specific")
plt3 = get_volcano_plot(df = unique_wilcox, pval_col = "bonferroni", p_title = "CD8+ Wilcoxon-specific")
