# used to backtest the similarities
# Revised Real Data Analysis Scripts
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(patchwork)
library(grid)
library(Seurat)
library(pbmcapply)
library(parallel)
library(doSNOW)
library(foreach) 
library(doParallel)
library(Hmisc)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) 
library(BiocParallel)
library(future)
library(future.apply)
library(progressr)
# ----- pre-regression -----
lib_size_correction_sparse = function(mat) {
  # mat: sparse G x C matrix
  lcf = Matrix::colMeans(mat)
  mean_lcf = mean(lcf)
  # Scale each column without converting to dense
  mat = Matrix::t(Matrix::t(mat) / lcf) * mean_lcf
  return(mat)
}
filter_SLE = function(mat, donor_ids, min_cells = 300, min_prop = 0.03) {
  # mat: sparse G x C matrix (genes x cells) for a single cell type
  # donor_ids: character vector of donor assignments, length = ncol(mat)
  # min_cells: minimum cells per donor (default 10)
  # min_prop: minimum donor-level non-zero expression rate per gene (default 3%)
  
  stopifnot(length(donor_ids) == ncol(mat))
  donor_ids = as.character(donor_ids)
  
  # filter donors with >= min_cells
  cell_counts = table(donor_ids)
  passing_donors = names(cell_counts)[cell_counts >= min_cells]
  cell_keep = donor_ids %in% passing_donors
  
  mat = mat[, cell_keep, drop = FALSE]
  donor_ids = donor_ids[cell_keep]
  
  # next for each gene, compute non-zero rate per donor
  # keep genes that pass the threshold in ALL remaining donors
  prop_mat = sapply(passing_donors, function(donor) {
    cell_idx = which(donor_ids == donor)
    submat = mat[, cell_idx, drop = FALSE]
    Matrix::rowSums(submat != 0) / length(cell_idx)
  })
  rownames(prop_mat) = rownames(mat)
  
  genes_keep = rownames(prop_mat)[apply(prop_mat >= min_prop, 1, all)]
  mat = mat[genes_keep, , drop = FALSE]
  
  return(list(
    matrix = mat,
    donors = passing_donors,
    genes = genes_keep
  ))
}
filter_SLE_donors = function(mat, donor_ids, min_cells = 300) {
  stopifnot(length(donor_ids) == ncol(mat))
  donor_ids   = as.character(donor_ids)
  cell_counts = table(donor_ids)
  donors_keep = names(cell_counts)[cell_counts >= min_cells]
  cell_keep   = donor_ids %in% donors_keep
  list(
    matrix     = mat[, cell_keep, drop = FALSE],
    donor_ids  = donor_ids[cell_keep],
    donors     = donors_keep
  )
}

filter_SLE_genes = function(mat, donor_ids, min_prop = 0.03) {
  stopifnot(length(donor_ids) == ncol(mat))
  donor_ids = as.character(donor_ids)
  donors    = unique(donor_ids)
  
  prop_mat = sapply(donors, function(d) {
    cell_idx = which(donor_ids == d)
    submat   = mat[, cell_idx, drop = FALSE]
    Matrix::rowSums(submat != 0) / length(cell_idx)
  })
  rownames(prop_mat) = rownames(mat)
  
  rownames(prop_mat)[apply(prop_mat >= min_prop, 1, all)]
}

filter_extreme_coverage <- function(sObj, n_mad = 3, log_scale = F, verbose = T) {
  
  counts_per_cell = sObj$nCount_RNA  # get per-cell total UMI counts from metadata
  
  x = if (log_scale) log1p(counts_per_cell) else counts_per_cell # allow for filtered counts on transformations
  
  # compute median and MAD on the chosen scale
  med  = median(x)
  madv = mad(x)
  
  if (madv == 0) { # exception handling when MAD = 0 (e.g. very few or very uniform cells)
    warning("MAD is 0 — distribution too uniform or too few cells. Returning object unfiltered.")
    return(sObj)
  }
  
  keep = abs(x - med) <= n_mad * madv
  
  sObj_filt = subset(sObj, cells = colnames(sObj)[keep]) # subset seurat objects
  
  if (verbose) {
    n_before = ncol(sObj)
    n_after  = ncol(sObj_filt)
    n_low    = sum(x < med - n_mad * madv)
    n_high   = sum(x > med + n_mad * madv)
    scale_lbl = if (log_scale) "log1p scale" else "raw scale"
    
    # print summary of filtered objects
    cat("Extreme coverage filter (", n_mad, " MAD, ", scale_lbl, ")\n", sep = "")
    cat("  Median:          ", round(med, 3), "\n", sep = "")
    cat("  MAD:             ", round(madv, 3), "\n", sep = "")
    cat("  Cells before:    ", n_before, "\n", sep = "")
    cat("  Cells after:     ", n_after, "\n", sep = "")
    cat("  Removed (low):   ", n_low, "\n", sep = "")
    cat("  Removed (high):  ", n_high, "\n", sep = "")
  }
  
  return(sObj_filt)
}
# ----- SEF helper functions -----
get_sc_bin_counts = function(Y, bin_edges) {
  stopifnot(is.list(Y), is.numeric(bin_edges), length(bin_edges) >= 2)
  if (is.unsorted(bin_edges, strictly = TRUE)) stop("bin_edges must be strictly increasing")
  
  K = length(bin_edges) - 1
  counts_matrix = matrix(0, nrow = length(Y), ncol = K)
  
  for (i in seq_along(Y)) {
    idx = cut(as.numeric(Y[[i]]), breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
    counts_matrix[i, ] = tabulate(idx, nbins = K)
  }
  return(counts_matrix)
}
safe_solve = function(A, B = NULL, ridge = 1e-8) {
  A_reg = A + diag(ridge, nrow(A))
  if (is.null(B)) return(solve(A_reg))
  solve(A_reg, B)
}
solve_with_fallback = function(A, B = NULL, ridge = 1e-8) {
  out = tryCatch(
    {
      if (is.null(B)) {
        solve(A)
      } else {
        solve(A, B)
      }
    },
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  safe_solve(A, B = B, ridge = ridge)
}
estimate_cov_sj = function(Smat) {
  n = nrow(Smat)
  K = ncol(Smat)
  m_vec = rowSums(Smat)
  bar_pooled = colSums(Smat) / sum(m_vec)
  inflate = if (n > 1) n / (n - 1) else 0
  
  lapply(seq_len(n), function(j) {
    s_j = Smat[j, ]
    m_j = m_vec[j]
    bar_j = if (m_j > 0) s_j / m_j else rep(0, K)
    # within = diag(s_j, K) - tcrossprod(s_j) / max(m_j, 1)
    if (m_j > 1) {
      within = (m_j / (m_j - 1)) * (diag(s_j, K) - tcrossprod(s_j) / m_j)
    } else {
      within = matrix(0, K, K)
    }
    between = inflate * (m_j^2) * tcrossprod(bar_j - bar_pooled)
    within + between
  })
}
estimate_cov_sj_one = function(Smat, j, m_vec = NULL, bar_pooled = NULL, inflate = NULL) {
  if (is.null(m_vec)) {
    m_vec = rowSums(Smat)
  }
  if (is.null(bar_pooled)) {
    bar_pooled = colSums(Smat) / sum(m_vec)
  }
  if (is.null(inflate)) {
    n = nrow(Smat)
    inflate = if (n > 1) n / (n - 1) else 0
  }
  
  K = ncol(Smat)
  s_j = Smat[j, ]
  m_j = m_vec[j]
  bar_j = if (m_j > 0) s_j / m_j else rep(0, K)
  # within = diag(s_j, K) - tcrossprod(s_j) / max(m_j, 1)
  if (m_j > 1) {
    within = (m_j / (m_j - 1)) * (diag(s_j, K) - tcrossprod(s_j) / m_j)
  } else {
    within = matrix(0, K, K)
  }
  between = inflate * (m_j^2) * tcrossprod(bar_j - bar_pooled)
  within + between
}

# ----- SEF regression functions -----
sef_regression_new = function(exprVec, sObj_meta, p = 2, plot_flag = F, K = "FD", bandwidth = 0.5, eps = 1e-10){
  # function: conducts regression for a given gene; smoothing done with std Gaussian KDE
  # output: returns "sef_regression_object" (list) that contains all objects used and computed for regression; retained for inference
  # check for subsetting
  stopifnot(length(exprVec) == nrow(sObj_meta))
  
  y_agg = as.vector(exprVec) # should already be a vector
  n_total_cells = length(y_agg)
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  
  control_donors = unique(sObj_meta[sObj_meta$disease == "normal", "donor_id"]) #group 1
  y1 = lapply(control_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  disease_donors = unique(sObj_meta[sObj_meta$disease != "normal", "donor_id"]) #group 2
  y2 = lapply(disease_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  if (is.null(K) | K == "FD") {
    K = nclass.FD(y_agg)
  }
  K = max(10, K)
  
  est_grid = seq(l, u, length.out = K + 1)
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
  binwidth = diff(est_grid)[1]
  Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid) # bin matrix K columns for group 1
  Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid) # bin matrix K columns for group 2
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  # Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  # Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)
  
  kernel_mat = outer(
    est_midpoints,
    est_midpoints,
    function(a, b) dnorm(a - b, sd = bandwidth)
  )
  smooth_mat = kernel_mat / pooled_total
  carrier = as.vector(smooth_mat %*% pooled_counts)
  if (sum(carrier == 0) > 0) {
    warning("carrier contains ", sum(carrier == 0), " zero values. adding value eps = ", eps)
  }
  carrier = pmax(carrier, eps)
  x_center = mean(est_midpoints)
  x_scale = sd(est_midpoints)
  if (!is.finite(x_scale) || x_scale == 0) {
    x_scale = 1
  }
  z_mid = (est_midpoints - x_center) / x_scale
  
  X = cbind(
    Intercept = 1,
    vapply(seq_len(p), function(d) z_mid^d, numeric(K))
  )
  colnames(X) = c("Intercept", paste0("X_", seq_len(p)))
  
  if (plot_flag) {
    print(
      ggplot(data.frame(x = est_midpoints, y = carrier), aes(x = x, y = y)) +
        geom_area(fill = "dodgerblue", alpha = 0.25) +
        geom_line(color = "dodgerblue", linewidth = 1.2) +
        labs(
          x = "Transformed Expression",
          y = "Density",
          title = "Carrier Density (updated covariance version)"
        ) +
        theme_minimal()
    )
  }
  
  #setup for modeling
  cellSum1 = sum(rowSums(Smat1))
  cellSum2 = sum(rowSums(Smat2))
  carrier_scale1 = carrier * binwidth * cellSum1
  carrier_scale2 = carrier * binwidth * cellSum2
  
  df1 = data.frame(sum_cts = S1sum, carrier_scale = carrier_scale1, X)
  df2 = data.frame(sum_cts = S2sum, carrier_scale = carrier_scale2, X)
  
  formula = as.formula(
    paste0(
      "sum_cts ~ offset(log(carrier_scale)) + ",
      paste(colnames(X)[-1], collapse = " + ")
    )
  )
  ctrl = glm.control(maxit = 300, epsilon = 1e-8)
  fit1 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df1, control = ctrl),
    error = function(e) NULL
  )
  fit2 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df2, control = ctrl),
    error = function(e) NULL
  )
  
  if (is.null(fit1) || is.null(fit2)) {
    return("GLM failed to fit for one or more groups.")
  }
  
  beta_est1 = as.vector(fit1$coefficients)
  beta_est2 = as.vector(fit2$coefficients)
  beta_diff = beta_est1 - beta_est2 # differences in coefficients
  e1 = as.vector(exp(X %*% beta_est1))
  e2 = as.vector(exp(X %*% beta_est2))
  
  
  sef_df1 = as.vector(carrier * e1 ) # get estimated group-wise density 
  sef_df2 = as.vector(carrier * e2 ) # get estimated group-wise density 
  
  if (plot_flag == T) {
    df_curve = data.frame(
      est_midpoints = rep(est_midpoints, 2),
      sef_value = c(carrier * e1, carrier * e2),
      group = factor(rep(c("Control", "SLE"), each = K))
    )
    df_hist = data.frame(
      value = y_agg,
      group = factor(c(rep("Control", cellSum1), rep("SLE", cellSum2)))
    )
    print(
      ggplot() +
        geom_histogram(
          data = df_hist,
          aes(x = value, y = after_stat(density), fill = group),
          bins = 50,
          alpha = 0.25,
          position = "identity",
          color = NA
        ) +
        geom_ribbon(
          data = df_curve,
          aes(x = est_midpoints, ymin = 0, ymax = sef_value, fill = group),
          alpha = 0.3,
          color = NA
        ) +
        geom_line(
          data = df_curve,
          aes(x = est_midpoints, y = sef_value, color = group),
          linewidth = 1
        ) +
        scale_color_manual(values = c("Control" = "firebrick", "SLE" = "dodgerblue")) +
        scale_fill_manual(values = c("Control" = "firebrick", "SLE" = "dodgerblue")) +
        labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
        theme_minimal()
    )
  }
  return(list(y1 = y1, y2 = y2, sef_df1 = sef_df1, sef_df2 = sef_df2, df1 = df1, df2 = df2, smooth_mat = smooth_mat,
              carrier = carrier, e1 = e1, e2 = e2, est_midpoints = est_midpoints, est_grid = est_grid,
              X = X, K = K, p = p, n_total_cells = n_total_cells, binwidth = binwidth, bandwidth = bandwidth,
              Smat1 = Smat1, Smat2 = Smat2, Smat = Smat, cellSum1 = cellSum1, cellSum2 = cellSum2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, fit1 = fit1, fit2 = fit2))
}

sef_inference_new = function(sef_regression_obj = NULL, ridge = 1e-8)
{
  # function: conducts inference after using SEF regression
  # input: sef_regression object that contains design matrix, bin count, fitted values, regression coefficients for the p moments tested, etc.
  # output: test statistic, p-value, covariance, matrix, difference in beta estimates
  X              = sef_regression_obj$X
  K              = sef_regression_obj$K
  n_total_cells  = sef_regression_obj$n_total_cells
  e1             = sef_regression_obj$e1
  e2             = sef_regression_obj$e2
  carrier        = sef_regression_obj$carrier
  smooth_mat     = sef_regression_obj$smooth_mat
  cellSum1       = sef_regression_obj$cellSum1
  cellSum2       = sef_regression_obj$cellSum2
  est_midpoints  = sef_regression_obj$est_midpoints
  binwidth       = sef_regression_obj$binwidth
  bandwidth      = sef_regression_obj$bandwidth
  beta_est1      = sef_regression_obj$beta_est1
  beta_est2      = sef_regression_obj$beta_est2
  Smat1          = sef_regression_obj$Smat1
  Smat2          = sef_regression_obj$Smat2
  p              = sef_regression_obj$p
  
  # validate that all required fields are present and non-NULL
  required_fields = c("X", "K","n_total_cells", "e1", "e2", "carrier", "smooth_mat",
                      "cellSum1", "cellSum2", "est_midpoints", "binwidth", "bandwidth",
                      "beta_est1", "beta_est2", "Smat1", "Smat2", "p")
  missing_fields = required_fields[sapply(required_fields, function(nm) {
    is.null(sef_regression_obj[[nm]])
  })]
  if (length(missing_fields) > 0) {
    stop("sef_regression_obj is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }
  
  beta_diff = beta_est2 - beta_est1
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  
  # The unadjusted specialization of the PDF formulas uses |Y_K| = binwidth.
  lambda1 = carrier * e1 * binwidth
  lambda2 = carrier * e2 * binwidth
  G_1 = t(X) %*% (lambda1 * X)
  G_2 = t(X) %*% (lambda2 * X)
  G1_inv = solve_with_fallback(G_1, ridge = ridge)
  G2_inv = solve_with_fallback(G_2, ridge = ridge)
  
  cov_list1 = estimate_cov_sj(Smat1)
  cov_list2 = estimate_cov_sj(Smat2)
  
  Z_t1_inT1 = (1 / cellSum1) * diag(K) - (e1 * binwidth) * smooth_mat
  Z_t1_inT2 = -(e1 * binwidth) * smooth_mat
  Z_t2_inT1 = -(e2 * binwidth) * smooth_mat
  Z_t2_inT2 = (1 / cellSum2) * diag(K) - (e2 * binwidth) * smooth_mat
  
  XZ_t1_inT1 = t(X) %*% Z_t1_inT1
  XZ_t1_inT2 = t(X) %*% Z_t1_inT2
  XZ_t2_inT1 = t(X) %*% Z_t2_inT1
  XZ_t2_inT2 = t(X) %*% Z_t2_inT2
  
  GXZ_t1_inT1 = G1_inv %*% XZ_t1_inT1
  GXZ_t1_inT2 = G1_inv %*% XZ_t1_inT2
  GXZ_t2_inT1 = G2_inv %*% XZ_t2_inT1
  GXZ_t2_inT2 = G2_inv %*% XZ_t2_inT2
  
  Mid_t1 = matrix(0, ncol(X), ncol(X))
  Mid_t2 = matrix(0, ncol(X), ncol(X))
  Cov_bar = matrix(0, ncol(X), ncol(X))
  
  A_diff_T1 = GXZ_t1_inT1 - GXZ_t2_inT1
  for (j in seq_len(nrow(Smat1))) {
    cj = cov_list1[[j]]
    Mid_t1 = Mid_t1 + XZ_t1_inT1 %*% cj %*% t(XZ_t1_inT1)
    Mid_t2 = Mid_t2 + XZ_t2_inT1 %*% cj %*% t(XZ_t2_inT1)
    Cov_bar = Cov_bar + A_diff_T1 %*% cj %*% t(A_diff_T1)
  }
  
  A_diff_T2 = GXZ_t1_inT2 - GXZ_t2_inT2
  for (j in seq_len(nrow(Smat2))) {
    cj = cov_list2[[j]]
    Mid_t1 = Mid_t1 + XZ_t1_inT2 %*% cj %*% t(XZ_t1_inT2)
    Mid_t2 = Mid_t2 + XZ_t2_inT2 %*% cj %*% t(XZ_t2_inT2)
    Cov_bar = Cov_bar + A_diff_T2 %*% cj %*% t(A_diff_T2)
  }
  Cov1 = G1_inv %*% Mid_t1 %*% G1_inv
  Cov2 = G2_inv %*% Mid_t2 %*% G2_inv
  Cov1 = 0.5 * (Cov1 + t(Cov1))
  Cov2 = 0.5 * (Cov2 + t(Cov2))
  Cov_bar = 0.5 * (Cov_bar + t(Cov_bar))
  
  stat_mat = Cov_bar[-1, -1, drop = FALSE]
  chi_stat1 = as.numeric(beta_diff[-1] %*% solve_with_fallback(stat_mat, beta_diff[-1], ridge = ridge))
  pval_1 = pchisq(chi_stat1, df = p, lower.tail = FALSE)
  
  return(list(beta_diff = beta_diff, pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_bar = Cov_bar))
}


run_sef_procedure = function(exprVec, sObj_meta, p = 2, plot_flag = F, K = "FD", bandwidth = 0.5, eps = 1e-10, ridge = 1e-8){
  # function: put the regression and inference functions together; wrapper for entire inferential and density estimation process
  # input: expression array; cell-wise metadata, named list of ID arrays; moments to test; plot flag; bandwidth
  # output: regression statistics, etc.
  regression_object = sef_regression_new(exprVec = exprVec, sObj_meta = sObj_meta, p = p, K = "FD",
                                         bandwidth = bandwidth, eps = eps, plot_flag = plot_flag)
  test_object = sef_inference_new(sef_regression_obj = regression_object, ridge = ridge)
  to_remove = c("X", "Smat1", "Smat2", "Smat", "fit1", "fit2", "df1", "df2","e1","e2","smooth_mat","carrier","est_grid","y1","y2") # removed more things optionally
  regression_statistics = regression_object[!(names(regression_object) %in% to_remove)]
  return(c(regression_statistics, test_object))
}

# ----- plotting -----
group_model_comparison_plot = function(sef_object, gene) {
  # function: plot a group comparison
  # input: sef object; gene of interest
  # output: plot comparing two densities; controls is blue; SLE is red
  
  sef_1 = sef_object[[gene]]$sef_df1
  sef_2 = sef_object[[gene]]$sef_df2
  est_midpoints = sef_object[[gene]]$est_midpoints
  df_plot = data.frame(
    x = rep(est_midpoints, 2),
    y = c(sef_1, sef_2),
    group = rep(c("Control", "SLE"), each = length(est_midpoints))
  )
  
  p = ggplot(df_plot, aes(x = x, y = y, color = group, fill = group)) + #plot using values from sef_df1, sef_df2
    geom_line(size = 1.2) +
    geom_area(alpha = 0.2, position = "identity") +
    scale_color_manual(values = c("Control" = "dodgerblue", "SLE" = "firebrick")) +
    scale_fill_manual(values = c("Control" = "dodgerblue", "SLE" = "firebrick")) +
    labs(
      title = paste("Group-Level Model Density:", gene),
      x = "Expression",
      y = "Density",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal(base_size = 14)
  
  return(p)
}
