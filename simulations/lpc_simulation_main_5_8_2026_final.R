# revised additional simulation functions MAIN
library(moments)

# ----- helper functions for SEF -----
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

estimate_cov_sj_adj = function(Smat, weights) {
  n = nrow(Smat)
  K = ncol(Smat)
  m_vec = rowSums(Smat)
  weighted_sum = colSums(sweep(Smat, 1, weights, "*"))
  weighted_total = sum(m_vec * weights)
  bar_adj = weighted_sum / weighted_total
  inflate = if (n > 1) n / (n - 1) else 0
  
  lapply(seq_len(n), function(j) {
    s_j = Smat[j, ]
    m_j = m_vec[j]
    bar_j = s_j / m_j
    #within = diag(s_j, K) - tcrossprod(s_j) / m_j
    if (m_j > 1) {
      within = (m_j / (m_j - 1)) * (diag(s_j, K) - tcrossprod(s_j) / m_j)
    } else {
      within = matrix(0, K, K)
    }
    between = inflate * (m_j^2) * tcrossprod(bar_j - bar_adj)
    within + between
  })
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

continuous_discrete_kde = function(y, eval_points, bw = 0.05) {
  # function: continuous gaussian kernel density estimate evaluated at the integer values the data is in 
  # input: expression array, points to evaluate (similar to est_midpoints)
  # output: density estimate evaluated at parameter-passed points
  dx = diff(eval_points)[1]  # step size between grid points; assumes equally spaced grid
  raw_density = sapply(eval_points, function(x) {
    mean(dnorm((x - y) / bw)) / bw
  })
  normalized_density = raw_density / sum(raw_density * dx)
  return(normalized_density)
}

lib_size_correction_sparse = function(mat) {
  # mat: sparse G x C matrix
  lcf = Matrix::colMeans(mat)
  mean_lcf = mean(lcf)
  # Scale each column without converting to dense
  mat = Matrix::t(Matrix::t(mat) / lcf) * mean_lcf
  return(mat)
}


# ----- data generation -----
generate_unimodal = function(n, k, mu = 0, sigma_sq = 1, sigmaC_sq = 0.2){
  stopifnot(length(k) == n, all(k > 0), all(k %% 1 == 0))
  y = rnorm(n, mu, sqrt(sigma_sq)) # generate n1 total group 1-specific mean effects
  y_list = vector("list", n)
  for(i in seq_len(n)){
    eps_ij = rnorm(k[i], 0, sqrt(sigmaC_sq)) 
    yj_vec = rep(y[i], k[i]) + eps_ij
    y_list[[i]] = yj_vec
  }
  return(y_list)
}

generate_bimodal_mixture = function(n, k, pi = 0.5, mu = 0.9, sigma_sq = 0.19, sigmaC_sq =0.2){
  stopifnot(length(k) == n, all(k > 0), all(k %% 1 == 0))
  z = rbinom(n, 1, pi) # binomial var for mixture
  y =  rnorm(n,mean = ifelse(z == 1, mu, -mu), sd = sqrt(sigma_sq)) # generate n1 total group 1-specific mean effects
  y_list = vector("list", n)
  for(i in seq_len(n)){
    eps_ij = rnorm(k[i], 0, sqrt(sigmaC_sq)) 
    yj_vec = rep(y[i], k[i]) + eps_ij
    y_list[[i]] = yj_vec
  }
  return(y_list)
}

generate_NB_fdr = function(n, k, mu = 10, theta = 5, sigma_sq = 0.5) {
  stopifnot(length(k) == n, all(k > 0), all(k %% 1 == 0))
  
  # Log-normal individual means, parameterized so E[mu_i] = mu
  # If log(mu_i) ~ Normal(m, sigma_sq), then E[mu_i] = exp(m + sigma_sq/2)
  # So set m = log(mu) - sigma_sq/2 to preserve the marginal mean
  log_mean = log(mu) - sigma_sq / 2
  mu_i = rlnorm(n, meanlog = log_mean, sdlog = sqrt(sigma_sq))
  
  y_list = vector("list", n)
  for (i in seq_len(n)) {
    y_list[[i]] = rnbinom(k[i], size = theta, mu = mu_i[i])
  }
  return(y_list)
}

# ----- empirical FDR test data generation together -----
get_FDR_TPR_metrics = function(pvals, labels, alpha = 0.05) {
  # labels: "DE" or "nonDE" i.e. sim$gene_info$de_status; ground truth
  
  DE_rej = pvals <= alpha # which genes found to be differentially distributed 
  
  TP = sum(labels == "DE" & DE_rej)
  FP = sum(labels == "nonDE" & DE_rej)
  FN = sum(labels == "DE" & !DE_rej)
  
  TPR = ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  FDR = ifelse((TP + FP) > 0, FP / (TP + FP), 0)
  
  return(list(TPR = TPR, FDR = FDR, TP = TP, FP = FP, FN = FN))
}

# ----- run sef regression functions -----
# note: scaled and unscaled will yield identical p-values and density estimates
run_sef_regression_covfix = function(gene, counts_vec, meta, gene_info = NULL, ridge = 1e-8,
                                     K = NULL, bw = 2, p = 2, eps = 1e-10, plt_flag = FALSE) 
{
  #input: gene name, count vector associated to gene, gene_info not needed
  cell_names = names(counts_vec)
  meta_g = meta[meta$cell_id %in% cell_names, , drop = FALSE]
  meta_g = meta_g[match(cell_names, meta_g$cell_id), , drop = FALSE]
  
  g1 = meta_g$group == "group1"
  g2 = meta_g$group == "group2"
  
  y1_list = split(
    counts_vec[meta_g$cell_id[g1]],
    meta_g$donor_id[g1],
    drop = TRUE
  )
  y2_list = split(
    counts_vec[meta_g$cell_id[g2]],
    meta_g$donor_id[g2],
    drop = TRUE
  )
  
  y_agg = c(unlist(y1_list, use.names = FALSE), unlist(y2_list, use.names = FALSE))
  l = min(y_agg)
  u = max(y_agg)
  
  if (is.null(K)) {
    K = nclass.FD(y_agg)
  }
  K = max(10, K)
  
  est_grid = seq(l, u, length.out = K + 1)
  est_midpoints = (est_grid[-1] + est_grid[-length(est_grid)]) / 2
  binwidth = diff(est_grid)[1]
  
  Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)
  
  # Keep the carrier on the density scale and carry binwidth explicitly
  # through the Poisson mean and sandwich terms.
  kernel_mat = outer(
    est_midpoints,
    est_midpoints,
    function(a, b) dnorm(a - b, sd = bw)
  )
  smooth_mat = kernel_mat / pooled_total
  carrier = as.vector(smooth_mat %*% pooled_counts)
  if (sum(carrier == 0) > 0) {
    warning("carrier_est contains ", sum(carrier == 0), " zero values. adding value eps = ", eps)
  }
  carrier = pmax(carrier, eps)
  
  X = cbind(
    Intercept = 1,
    vapply(seq_len(p), function(d) est_midpoints^d, numeric(K))
  )
  colnames(X) = c("Intercept", paste0("X_", seq_len(p)))
  
  if (plt_flag) {
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
  
  fail_return = function(msg) {
    warning("Gene ", gene, ": ", msg)
    list(
      gene = gene,
      pval_1 = NA_real_,
      beta_est1 = NULL,
      beta_est2 = NULL,
      beta_diff = NULL,
      Cov_bar = NULL,
      Cov1 = NULL,
      Cov2 = NULL,
      K = K,
      kde_used = "gaussian-kernel-matrix-covfix",
      p = p,
      binwidth = binwidth,
      bw = bw
    )
  }
  
  fit1 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df1, control = ctrl),
    error = function(e) NULL
  )
  fit2 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df2, control = ctrl),
    error = function(e) NULL
  )
  
  if (is.null(fit1) || is.null(fit2)) {
    return(fail_return("GLM failed to fit for one or more groups."))
  }
  
  beta_est1 = as.vector(fit1$coefficients)
  beta_est2 = as.vector(fit2$coefficients)
  beta_diff = beta_est1 - beta_est2
  e1 = as.vector(exp(X %*% beta_est1))
  e2 = as.vector(exp(X %*% beta_est2))
  
  plt = NULL
  if (plt_flag) {
    df_curve = data.frame(
      est_midpoints = rep(est_midpoints, 2),
      sef_value = c(carrier * e1, carrier * e2),
      group = factor(rep(c("Group 1", "Group 2"), each = K))
    )
    df_hist = data.frame(
      value = y_agg,
      group = factor(c(rep("Group 1", cellSum1), rep("Group 2", cellSum2)))
    )
    plt = ggplot() +
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
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(
      plt
    )
  }
  
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
  
  list(
    gene = gene,
    pval_1 = pval_1,
    beta_est1 = beta_est1,
    beta_est2 = beta_est2,
    beta_diff = beta_diff,
    Cov_bar = Cov_bar,
    Cov1 = Cov1,
    Cov2 = Cov2,
    chi_stat1 = chi_stat1,
    K = K,
    kde_used = "gaussian-kernel-matrix-covfix",
    p = p,
    binwidth = binwidth,
    bw = bw,
    carrier = carrier,
    plt = plt
  )
}

run_sef_test_covfix = function(sim_obj, K = NULL, bw = 0.5, p = 2, verbose = TRUE, ridge = 1e-8) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)
  
  pVec = setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    if (verbose) {
      message(i, " of ", length(genes), " (covfix)")
    }
    out = run_sef_regression_covfix(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      K = K,
      ridge = ridge,
      bw = bw,
      p = p
    )
    pVec[i] = out$pval_1
  }
  
  cbind(gene_info, pVec = unname(pVec[genes]))
}

run_sef_regression_covfix_scaled_basis = function(gene, counts_vec, meta, gene_info = NULL, ridge = 1e-8,
                                                  K = NULL, bw = 2, p = 2, eps = 1e-10, plt_flag = FALSE) {
  # centered method helps enable covariance
  cell_names = names(counts_vec)
  meta_g = meta[meta$cell_id %in% cell_names, , drop = FALSE]
  meta_g = meta_g[match(cell_names, meta_g$cell_id), , drop = FALSE]
  
  g1 = meta_g$group == "group1"
  g2 = meta_g$group == "group2"
  
  y1_list = split(
    counts_vec[meta_g$cell_id[g1]],
    meta_g$donor_id[g1],
    drop = TRUE
  )
  y2_list = split(
    counts_vec[meta_g$cell_id[g2]],
    meta_g$donor_id[g2],
    drop = TRUE
  )
  
  y_agg = c(unlist(y1_list, use.names = FALSE), unlist(y2_list, use.names = FALSE))
  l = min(y_agg)
  u = max(y_agg)
  
  if (is.null(K)) {
    K = nclass.FD(y_agg)
  }
  K = max(10, K)
  
  est_grid = seq(l, u, length.out = K + 1)
  est_midpoints = (est_grid[-1] + est_grid[-length(est_grid)]) / 2
  binwidth = diff(est_grid)[1]
  
  Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)
  
  kernel_mat = outer(
    est_midpoints,
    est_midpoints,
    function(a, b) dnorm(a - b, sd = bw)
  )
  smooth_mat = kernel_mat / pooled_total
  carrier = as.vector(smooth_mat %*% pooled_counts)
  if (sum(carrier == 0) > 0) {
    warning("carrier_est contains ", sum(carrier == 0), " zero values. adding value eps = ", eps)
  }
  carrier = pmax(carrier, eps)
  
  # keep plotting/evaluation on the original midpoint scale
  x_mid = est_midpoints
  
  # scale only the regression basis
  x_center = mean(x_mid)
  x_scale = sd(x_mid)
  if (!is.finite(x_scale) || x_scale == 0) {
    x_scale = 1
  }
  z_mid = (x_mid - x_center) / x_scale
  
  X = cbind(
    Intercept = 1,
    vapply(seq_len(p), function(d) z_mid^d, numeric(K))
  )
  colnames(X) = c("Intercept", paste0("X_", seq_len(p)))
  
  if (plt_flag) {
    print(
      ggplot(data.frame(x = est_midpoints, y = carrier), aes(x = x, y = y)) +
        geom_area(fill = "dodgerblue", alpha = 0.25) +
        geom_line(color = "dodgerblue", linewidth = 1.2) +
        labs(
          x = "Transformed Expression",
          y = "Density",
          title = "Carrier Density"
        ) +
        theme_minimal()
    )
  }
  
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
  
  fail_return = function(msg) {
    warning("Gene ", gene, ": ", msg)
    list(
      gene = gene,
      pval_1 = NA_real_,
      beta_est1 = NULL,
      beta_est2 = NULL,
      beta_diff = NULL,
      Cov_bar = NULL,
      Cov1 = NULL,
      Cov2 = NULL,
      K = K,
      kde_used = "gaussian-kernel-matrix-covfix-scaled-basis",
      p = p,
      binwidth = binwidth,
      bw = bw
    )
  }
  
  fit1 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df1, control = ctrl),
    error = function(e) NULL
  )
  fit2 = tryCatch(
    glm(formula, family = poisson(link = "log"), data = df2, control = ctrl),
    error = function(e) NULL
  )
  
  if (is.null(fit1) || is.null(fit2)) {
    return(fail_return("GLM failed to fit for one or more groups."))
  }
  
  beta_est1 = as.vector(fit1$coefficients)
  beta_est2 = as.vector(fit2$coefficients)
  beta_diff = beta_est1 - beta_est2
  
  e1 = as.vector(exp(X %*% beta_est1))
  e2 = as.vector(exp(X %*% beta_est2))
  
  sef_df1 = carrier * e1
  sef_df2 = carrier * e2
  
  plt = NULL
  if (plt_flag) {
    df_curve = data.frame(
      est_midpoints = rep(est_midpoints, 2),
      sef_value = c(sef_df1, sef_df2),
      group = factor(rep(c("Group 1", "Group 2"), each = K))
    )
    df_hist = data.frame(
      value = y_agg,
      group = factor(c(rep("Group 1", cellSum1), rep("Group 2", cellSum2)))
    )
    
    plt = ggplot() +
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
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    
    print(plt)
  }
  
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
  chi_stat1 = as.numeric(
    beta_diff[-1] %*% solve_with_fallback(stat_mat, beta_diff[-1], ridge = ridge)
  )
  pval_1 = pchisq(chi_stat1, df = p, lower.tail = FALSE)
  
  list(
    gene = gene,
    pval_1 = pval_1,
    beta_est1 = beta_est1,
    beta_est2 = beta_est2,
    beta_diff = beta_diff,
    Cov_bar = Cov_bar,
    Cov1 = Cov1,
    Cov2 = Cov2,
    chi_stat1 = chi_stat1,
    K = K,
    kde_used = "gaussian-kernel-matrix-covfix-scaled-basis",
    p = p,
    binwidth = binwidth,
    bw = bw,
    carrier = carrier,
    est_midpoints = est_midpoints,
    z_midpoints = z_mid,
    basis_center = x_center,
    basis_scale = x_scale,
    sef_df1 = sef_df1,
    sef_df2 = sef_df2,
    plt = plt
  )
}
run_sef_test_covfix_scaled_basis = function(sim_obj, K = NULL, bw = 0.5, p = 2, verbose = TRUE, ridge = 1e-8) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)
  
  pVec = setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    if (verbose) {
      message(i, " of ", length(genes), " (covfix, scaled basis)")
    }
    out = run_sef_regression_covfix_scaled_basis(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      K = K,
      ridge = ridge,
      bw = bw,
      p = p
    )
    pVec[i] = out$pval_1
  }
  
  cbind(gene_info, pVec = unname(pVec[genes]))
}

# ----- comparison methods -----
get_pseudobulk = function(counts, meta, prep_for_de = FALSE) {
  # input checking
  if (is.null(rownames(counts))) {
    stop("counts must have rownames (cell IDs).")
  }
  if (!all(c("cell_id", "donor_id") %in% colnames(meta))) {
    stop("meta must contain columns: cell_id and donor_id.")
  }
  if (anyDuplicated(meta$cell_id)) {
    stop("meta$cell_id contains duplicates.")
  }
  
  common_cells = intersect(rownames(counts), meta$cell_id)
  if (length(common_cells) == 0) {
    stop("No overlapping cell IDs between counts and meta.")
  }
  
  
  counts2 = counts[common_cells, , drop = FALSE]
  meta2 = meta[match(common_cells, meta$cell_id), , drop = FALSE]
  
  keep = !is.na(meta2$donor_id)
  counts2 = counts2[keep, , drop = FALSE]
  donors = meta2$donor_id[keep]
  
  if (nrow(counts2) == 0) {
    stop("No cells left after removing NA donor_id.")
  }
  
  counts_pseudobulk = rowsum(counts2, group = donors, reorder = FALSE)
  
  if (prep_for_de) {
    modal_cols = grep("^modal", colnames(counts_pseudobulk), value = TRUE)
    
    if (length(modal_cols) > 0) {
      nonint = vapply(modal_cols, function(cl) { # loop through columns of the modality-shift DEGs
        x = counts_pseudobulk[, cl]
        any(abs(x - round(x)) > 1e-8, na.rm = TRUE)
      }, logical(1))
      
      to_round = modal_cols[nonint]
      
      if (length(to_round) > 0) {
        counts_pseudobulk[, to_round] = round(counts_pseudobulk[, to_round, drop = FALSE])
        
        global_min = min(counts_pseudobulk[, to_round, drop = FALSE], na.rm = TRUE)
        if (is.finite(global_min) && global_min < 0) {
          offset = as.integer(-global_min)
          counts_pseudobulk[, to_round] = counts_pseudobulk[, to_round, drop = FALSE] + offset
        }
      }
    }
  }
  
  return(t(counts_pseudobulk)) #genes as rows, samples as columns
}

run_pseudobulk_methods = function(counts_pb, gene_info, method = "edgeR"){
  genes = gene_info$gene
  DE_labels = gene_info$de_status
  donors = unique(colnames(counts_pb))
  prefix = substr(donors, 1, 2)
  stopifnot(all(prefix %in% c("G1", "G2")))
  groups_pb = ifelse(prefix == "G1", 1, 2)
  groups_pb = factor(groups_pb)
  
  if(method == "edgeR"){
    # reference source for proper comparison: https://github.com/xihuimeijing/DEGs_Analysis_FDR/blob/main/scripts/DEGs.R
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    y = edgeR::DGEList(counts = counts_pb, group = groups_pb)
    keep_gene = edgeR::filterByExpr(y)
    y = y[keep_gene, , keep.lib.sizes = FALSE]
    y = edgeR::calcNormFactors(y, method = "TMM")
    
    design = model.matrix(~ groups_pb)
    y = edgeR::estimateDisp(y, design)
    fit = edgeR::glmQLFit(y, design)
    qlf = edgeR::glmQLFTest(fit, coef = 2)
    
    res = edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    return(res)
    
  } else if(method == "DESeq2"){
    # reference source: https://github.com/xihuimeijing/DEGs_Analysis_FDR/blob/main/scripts/DEGs.R
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Please install DESeq2.")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    
    keep_gene = rowSums(counts_pb) > 0
    counts_pb = counts_pb[keep_gene, , drop = FALSE]
    
    coldata = data.frame(groups_pb = groups_pb)
    rownames(coldata) = colnames(counts_pb)
    
    dds = DESeq2::DESeqDataSetFromMatrix(
      countData = counts_pb,
      colData = coldata,
      design = ~ groups_pb
    )
    dds = DESeq2::DESeq(dds, quiet = TRUE)
    res = as.data.frame(DESeq2::results(dds))
    colnames(res)[colnames(res) == "padj"] = "FDR"
    return(res)
    
  } else if(method == "wilcox"){
    # Tutorial followed: https://rpubs.com/LiYumei/806213 to use edgeR style filtering and TMM normalization
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR for TMM normalization")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    y = edgeR::DGEList(counts = counts_pb, group = groups_pb)
    keep_gene = edgeR::filterByExpr(y)
    y = y[keep_gene, , keep.lib.sizes = FALSE]
    y = edgeR::calcNormFactors(y, method = "TMM")
    count_norm = edgeR::cpm(y)
    
    g1 = which(groups_pb == levels(groups_pb)[1])
    g2 = which(groups_pb == levels(groups_pb)[2])
    stopifnot(length(g1) > 0, length(g2) > 0)
    
    pvals = sapply(seq_len(nrow(count_norm)), function(i) {
      data = data.frame(
        gene = as.numeric(count_norm[i, ]),
        conditions = groups_pb
      )
      if (length(unique(data$gene)) <= 1) return(1)
      suppressWarnings(wilcox.test(gene ~ conditions, data = data, exact = FALSE)$p.value)
    })
    mean_g1 = rowMeans(count_norm[, g1, drop = FALSE])
    mean_g2 = rowMeans(count_norm[, g2, drop = FALSE])
    res = data.frame(
      gene = rownames(count_norm),
      logFC = log2(mean_g2/mean_g1),
      PValue = pvals,
      FDR = p.adjust(pvals, method = "fdr"),
      stringsAsFactors = FALSE
    )
    rownames(res) = res$gene
    return(res)
  }
  else{
    stop(sprintf("Error: invalid method '%s'. Valid methods are: edgeR, DESeq2, wilcoxon.", method))
  }
}

run_pseudobulk_methods_fix_size_factors = function(counts_pb, gene_info, method = "edgeR"){
  genes = gene_info$gene
  DE_labels = gene_info$de_status
  donors = unique(colnames(counts_pb))
  prefix = substr(donors, 1, 2)
  stopifnot(all(prefix %in% c("G1", "G2")))
  groups_pb = ifelse(prefix == "G1", 1, 2)
  groups_pb = factor(groups_pb)
  
  if(method == "edgeR"){
    # reference source for proper comparison: https://github.com/xihuimeijing/DEGs_Analysis_FDR/blob/main/scripts/DEGs.R
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    y = edgeR::DGEList(counts = counts_pb, group = groups_pb)
    keep_gene = edgeR::filterByExpr(y)
    y = y[keep_gene, , keep.lib.sizes = FALSE]
    y = edgeR::calcNormFactors(y, method = "none")   # sets norm.factors = 1
    
    design = model.matrix(~ groups_pb)
    y = edgeR::estimateDisp(y, design)
    fit = edgeR::glmQLFit(y, design)
    qlf = edgeR::glmQLFTest(fit, coef = 2)
    
    res = edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    return(res)
    
  } else if(method == "DESeq2"){
    # reference source: https://github.com/xihuimeijing/DEGs_Analysis_FDR/blob/main/scripts/DEGs.R
    if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Please install DESeq2.")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    
    keep_gene = rowSums(counts_pb) > 0
    counts_pb = counts_pb[keep_gene, , drop = FALSE]
    
    coldata = data.frame(groups_pb = groups_pb)
    rownames(coldata) = colnames(counts_pb)
    
    dds = DESeq2::DESeqDataSetFromMatrix(
      countData = counts_pb,
      colData = coldata,
      design = ~ groups_pb
    )
    libsize = colSums(counts_pb)
    sf = libsize / exp(mean(log(libsize)))   # edgeR method="none" analog
    DESeq2::sizeFactors(dds) = sf
    
    dds = DESeq2::DESeq(dds, quiet = TRUE)
    res = as.data.frame(DESeq2::results(dds))
    colnames(res)[colnames(res) == "padj"] = "FDR"
    return(res)
    
  } else if(method == "wilcox"){
    # Tutorial followed: https://rpubs.com/LiYumei/806213 to use edgeR style filtering and TMM normalization
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install edgeR for TMM normalization")
    stopifnot(identical(rownames(counts_pb), as.character(gene_info$gene)))
    y = edgeR::DGEList(counts = counts_pb, group = groups_pb)
    keep_gene = edgeR::filterByExpr(y)
    y = y[keep_gene, , keep.lib.sizes = FALSE]
    y = edgeR::calcNormFactors(y, method = "none")
    count_norm = edgeR::cpm(y)
    
    g1 = which(groups_pb == levels(groups_pb)[1])
    g2 = which(groups_pb == levels(groups_pb)[2])
    stopifnot(length(g1) > 0, length(g2) > 0)
    
    pvals = sapply(seq_len(nrow(count_norm)), function(i) {
      data = data.frame(
        gene = as.numeric(count_norm[i, ]),
        conditions = groups_pb
      )
      if (length(unique(data$gene)) <= 1) return(1)
      suppressWarnings(wilcox.test(gene ~ conditions, data = data, exact = FALSE)$p.value)
    })
    mean_g1 = rowMeans(count_norm[, g1, drop = FALSE])
    mean_g2 = rowMeans(count_norm[, g2, drop = FALSE])
    res = data.frame(
      gene = rownames(count_norm),
      logFC = log2(mean_g2/mean_g1),
      PValue = pvals,
      FDR = p.adjust(pvals, method = "fdr"),
      stringsAsFactors = FALSE
    )
    rownames(res) = res$gene
    return(res)
  }
  else{
    stop(sprintf("Error: invalid method '%s'. Valid methods are: edgeR, DESeq2, wilcoxon.", method))
  }
}


# ----- plotting -----
gg_qqplot = function(ps, ci = 0.95) {
  # function: plotting qq plot with 95% confidence interval band using input of raw p-values
  # input: p-value array; confidence interval 
  # output: QQ plot
  n  = length(ps)
  df = data.frame(
    observed = -log10(sort(ps)), 
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe = expression(paste("Expected -log"[10], plain(P)))
  log10Po = expression(paste("Observed -log"[10], plain(P)))
  max_lim = max(df$expected, df$observed, df$clower, df$cupper, na.rm = TRUE) 
  
  plt = ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1, fill = "steelblue"
    ) +
    geom_point(aes(expected, observed), shape = 21, size = 3, color = "black", fill = "steelblue" ) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    theme_bw(base_size = 22) +
    xlab(log10Pe) + ylab(log10Po) 
  return(plt)
} 
