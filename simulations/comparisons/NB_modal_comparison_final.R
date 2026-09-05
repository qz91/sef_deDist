# NB Modal comparison using pseudobulk methods
rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)
source("revisions_simulation_main_5_8_2026_final.R")
source("main_NB_4_30_2026_final.R")
source("lpc_NB_modality_analysis_final.R")

# ----- param instantiation -----
n1 = 100; n2 = 100
lower = 500
upper = 1000

n_null_unimodal = 1000
n_null_bimodal = 1000
n_modal = 200

target_mu_range = c(100, 100)
delta_range = c(30, 30)
target_marginal_variance = 1150

theta_mix_range = c(100, 100)
mixing_prob = 0.5
sigma_sq = 0.01

modal_unimodal_group = "group1"
return_donor_effects = TRUE
message("beginning modality shifted pseudobulk analysis...")
out.dir = file.path("")
# ----- run_one function wrapper based on seed using same params -----
run_one_pb = function(seed = 1) { 
  sim_obj = simulate_NB_modality_analysis(
    seed = seed,
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    n_null_unimodal = n_null_unimodal,
    n_null_bimodal = n_null_bimodal,
    n_modal = n_modal,
    target_mu_range = target_mu_range,
    delta_range = delta_range,
    target_marginal_variance = target_marginal_variance,
    theta_mix_range = theta_mix_range,
    mixing_prob = mixing_prob,
    sigma_sq = sigma_sq,
    modal_unimodal_group = modal_unimodal_group,
    return_donor_effects = return_donor_effects
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
  message(sprintf("seed %02d pseudobulk DE analysis for modality shift complete", seed))
  return(out)
}
## ----- get pb objects -----
testPB = pblapply(seq_len(50), function(s) run_one_pb(seed = s))
fname = paste0(
  "NB_modal_pb_metrics_6_16_2026.RDS"
)
out.file = file.path(out.dir, fname)
saveRDS(testPB, out.file)

## ----- get pb avgs -----
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
fname2 = paste0(
  "NB_modal_avg_pb_metrics_6_16_2026.RDS"
)
out.file2 = file.path(out.dir, fname2)
saveRDS(pb_avg_metrics, out.file2)

