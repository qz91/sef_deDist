rm(list = ls())
library(pbapply)
library(moments)
library(DESeq2)
library(edgeR)

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
shift = "mean"
DE_mixture_prob_range = c(0.30, 0.70)
share_mode = TRUE # unimodal group shares the first mixture component mean

# Keep this reproducibility run's outputs with the submission code.
out.dir = file.path("")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
message("beginning strong mean shift pseudobulk analysis...")

# ----- run_one function wrapper based on seed using same params -----
run_one_pb = function(seed = 1) { 
  sim_obj = generate_NB_gene_data(
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    n_nonDE_unimodal = n_nonDE_unimodal, 
    n_nonDE_mixture = n_nonDE_mixture, 
    n_DE = n_DE, 
    nonDE_unimodal_mu_range = nonDE_unimodal_mu_range, nonDE_unimodal_theta_range = nonDE_unimodal_theta_range,
    nonDE_mixture_mu_range = nonDE_mixture_mu_range, nonDE_mixture_theta_range =nonDE_mixture_theta_range, 
    nonDE_mixture_fc_range = nonDE_mixture_fc_range, nonDE_mixture_prob_range = nonDE_mixture_prob_range,
    DE_mu_range = DE_mu_range,
    DE_theta_range = DE_theta_range,
    DE_mixture_fc_range = DE_mixture_fc_range,
    DE_mixture_prob_range = DE_mixture_prob_range,
    sigma_sq_nonDE = sigma_sq_nonDE,
    sigma_sq_DE = sigma_sq_DE,
    sigma_sq = sigma_sq,
    de_unimodal_group = de_unimodal_group,
    share_mode = share_mode,
    seed = seed
  )
  message(sprintf("seed %02d data generation complete", seed))
  flush.console()
  counts_pb = get_pseudobulk(sim_obj$count_matrix, sim_obj$metadata, prep_for_de = F) # aggregate cells together
  message(sprintf("seed %02d pseudobulking complete", seed))
  
  make_test_full = function(res) {
    if (!"gene" %in% colnames(res)) {
      res$gene = rownames(res)
    }
    
    if (!"p_raw" %in% colnames(res)) {
      if ("PValue" %in% colnames(res)) {
        res$p_raw = res$PValue
      } else if ("pvalue" %in% colnames(res)) {
        res$p_raw = res$pvalue
      }
    }
    
    if (!"FDR" %in% colnames(res) && "padj" %in% colnames(res)) {
      res$FDR = res$padj
    }
    
    if (!"p_raw" %in% colnames(res)) {
      stop("Result must contain raw p-values (PValue or pvalue).")
    }
    if (!"FDR" %in% colnames(res)) {
      stop("Result must contain FDR (or padj).")
    }
    
    out = merge(
      sim_obj$gene_info[, c("gene", "de_status")],
      res[, c("gene", "p_raw", "FDR")],
      by = "gene",
      all.x = TRUE,
      sort = FALSE
    )
    
    out$p_raw[is.na(out$p_raw)] = 1
    out$FDR[is.na(out$FDR)] = 1
    out
  }
  pb_methods = c("edgeR", "DESeq2", "wilcox")
  results_by_pb_method = setNames(lapply(pb_methods, function(method) {
    res = run_pseudobulk_methods(
      counts_pb,
      gene_info = sim_obj$gene_info,
      method = method
    )
    test = make_test_full(res)
    metrics = get_FDR_TPR_metrics(test$FDR, test$de_status)
    
    message(sprintf("seed %02d finished %s", seed, method))
    flush.console()
    
    list(test = test, metrics = metrics)
  }), pb_methods)
  
  out = list(
    seed = seed,
    results = results_by_pb_method
  )
  message(sprintf("seed %02d pseudobulk DE analysis complete", seed))
  return(out)
}

# ----- get pb objects -----
# should output testPB[[seed]]$results[[method]]$test / $metrics, and pb_avg_metrics
testPB = pblapply(seq_len(n_sim), function(s) run_one_pb(seed = s))
fname = "NB_strong_mean_shift_pb_metrics_5_22_2026.RDS"
out.file = file.path(out.dir, fname)
saveRDS(testPB, out.file)

# ----- get pb avgs -----
get_avg_pb_metrics = function(x, methods = NULL, metrics = c("TPR", "FDR")) {
  if (is.null(methods)) {
    methods <- names(x[[1]]$results)
  }
  
  out <- do.call(rbind, lapply(methods, function(method) {
    seed_df <- do.call(rbind, lapply(x, function(seed) {
      df <- data.frame(seed$results[[method]]$metrics)
      df[, metrics, drop = FALSE]
    }))
    
    avgs <- colMeans(seed_df, na.rm = TRUE)
    
    data.frame(
      method = method,
      t(avgs),
      row.names = NULL
    )
  }))
  
  rownames(out) <- NULL
  out
}

pb_avg_metrics = get_avg_pb_metrics(testPB)
fname2 = "NB_strong_mean_shift_avg_pb_metrics_5_22_2026.RDS"
out.file2 = file.path(out.dir, fname2)
saveRDS(pb_avg_metrics, out.file2)
print(pb_avg_metrics)
