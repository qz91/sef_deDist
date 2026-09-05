rm(list = ls())
library(pbapply)
library(moments)
base_dir = file.path("")
source(file.path(base_dir, "lpc_simulation_main_5_8_2026_final.R"))
source(file.path(base_dir, "lpc_main_NB_4_30_2026_final.R"))

# ----- param instantiation -----
n_sim = 50
n1 = 100; n2 = 100
lower = 500
upper = 1000
n_nonDE_unimodal = 1000
n_nonDE_mixture = 1000
n_DE = 200
nonDE_unimodal_mu_range = c(5, 20)
nonDE_unimodal_theta_range = c(3, 10)
nonDE_mixture_mu_range = c(5, 15)
nonDE_mixture_theta_range = c(3, 10)
nonDE_mixture_fc_range = c(2, 4)
nonDE_mixture_prob_range = c(0.30, 0.70)
DE_mu_range = c(5, 15)
DE_theta_range = c(3, 10)
DE_mixture_fc_range = c(1.25, 2.25)
sigma_sq_nonDE = 0.05
sigma_sq = 0.05
sigma_sq_DE = 0.05
de_unimodal_group = "group1"
shift = "variance"
DE_mixture_prob_range = c(0.5, 0.5)
share_mode = FALSE # midpoint matches the equally weighted mixture mean

out.dir = file.path(base_dir, "sef_comparisons", "results", "NB_small_var")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
message("beginning small variance shift Wilcoxon variance tests...")

# ----- testing donor-specific variance: original run_one_vartest -----
run_one_vartest = function(seed = 1) {
  # wilcox test of sample variances
  sim_obj = generate_NB_gene_data(
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    n_nonDE_unimodal = n_nonDE_unimodal,
    n_nonDE_mixture = n_nonDE_mixture,
    n_DE = n_DE,
    nonDE_unimodal_mu_range = nonDE_unimodal_mu_range, nonDE_unimodal_theta_range = nonDE_unimodal_theta_range,
    nonDE_mixture_mu_range = nonDE_mixture_mu_range, nonDE_mixture_theta_range = nonDE_mixture_theta_range,
    nonDE_mixture_fc_range = nonDE_mixture_fc_range, nonDE_mixture_prob_range = nonDE_mixture_prob_range,
    DE_mu_range = DE_mu_range,
    DE_mixture_fc_range = DE_mixture_fc_range,
    DE_mixture_prob_range = DE_mixture_prob_range,
    sigma_sq_nonDE = sigma_sq_nonDE,
    sigma_sq_DE = sigma_sq_DE,
    sigma_sq = sigma_sq,
    de_unimodal_group = "group1",
    share_mode = share_mode,
    seed = seed
  )
  message(sprintf("seed %02d data generation complete", seed))
  flush.console()
  gene_info = sim_obj$gene_info
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  meta = meta[match(rownames(counts), meta$cell_id), , drop = FALSE]
  stopifnot(all(meta$cell_id == rownames(counts)))

  donor_idx = split(seq_len(nrow(counts)), meta$donor_id)
  donors = names(donor_idx)

  variance_mat = t(vapply(donors, function(d) {
    idx = donor_idx[[d]]
    apply(counts[idx, , drop = FALSE], 2, var)
  }, numeric(ncol(counts))))

  colnames(variance_mat) = colnames(counts)
  rownames(variance_mat) = donors

  group_df = unique(meta[, c("donor_id", "group")])
  group_df = group_df[match(rownames(variance_mat), group_df$donor_id), , drop = FALSE]

  variance_df = data.frame(
    donor_id = rep(rownames(variance_mat), times = ncol(variance_mat)),
    group = rep(group_df$group, times = ncol(variance_mat)),
    gene = rep(colnames(variance_mat), each = nrow(variance_mat)),
    donor_level_variance = as.vector(variance_mat),
    stringsAsFactors = FALSE
  )
  genes = unique(variance_df$gene)

  pvals = sapply(genes, function(g) {
    df_g = variance_df[variance_df$gene == g, , drop = FALSE]

    x1 = df_g$donor_level_variance[df_g$group == "group1"]
    x2 = df_g$donor_level_variance[df_g$group == "group2"]

    if (length(x1) == 0 || length(x2) == 0) return(NA_real_)
    if (length(unique(c(x1, x2))) <= 1) return(1)

    wilcox.test(x1, x2, exact = FALSE)$p.value
  })

  res_var_wilcox = data.frame(
    gene = genes,
    PValue = pvals,
    FDR = p.adjust(pvals, method = "BH"),
    stringsAsFactors = FALSE
  )

  rownames(res_var_wilcox) = res_var_wilcox$gene
  labels = gene_info$de_status[match(res_var_wilcox$gene, gene_info$gene)]

  metrics = get_FDR_TPR_metrics(
    pvals = res_var_wilcox$FDR,
    labels = labels
  )
  out = list(test = res_var_wilcox, metrics = metrics)
  return(out)
}

# ----- get vartest objects: preserve the original unnamed list structure -----
test_vartest = pblapply(seq_len(n_sim), function(s) run_one_vartest(seed = s))
fname = "NB_vartest_small_effect_result_5_26_2026.RDS"
out.file = file.path(out.dir, fname)
saveRDS(test_vartest, out.file)
