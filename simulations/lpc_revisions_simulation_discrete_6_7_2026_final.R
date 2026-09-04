# Discrete singleton SEF variants.
#
# These functions depend on helpers defined in revisions_simulation_main_5_8_2026.R:
# get_sc_bin_counts(), estimate_cov_sj(), and solve_with_fallback().

run_sef_regression_discrete_covfix = function(gene, counts_vec, meta, gene_info = NULL,
                                              ridge = 1e-8, K = NULL, bw = NULL,
                                              p = 2, eps = 1e-10, plt_flag = FALSE)
{
  # bw is kept only for drop-in compatibility with the KDE-based functions.
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
  y_agg = as.numeric(y_agg)

  Kmax_obs = ceiling(max(y_agg, na.rm = TRUE))
  if (is.null(K)) {
    Kmax = Kmax_obs
  } else {
    Kmax = max(as.integer(K), Kmax_obs)
  }
  Kmax = max(0L, Kmax)

  est_midpoints = seq.int(0L, Kmax)
  est_grid = seq(-0.5, Kmax + 0.5, by = 1)
  K = length(est_midpoints)
  binwidth = 1

  Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)

  # Singleton PMF carrier: one count at k affects only the mass at k.
  smooth_mat = diag(K) / pooled_total
  carrier = as.vector(smooth_mat %*% pooled_counts)
  # Zero-mass singleton bins are expected when some integer counts are unobserved.
  carrier = pmax(carrier, eps)

  X = cbind(
    Intercept = 1,
    vapply(seq_len(p), function(d) est_midpoints^d, numeric(K))
  )
  colnames(X) = c("Intercept", paste0("X_", seq_len(p)))

  if (plt_flag) {
    print(
      ggplot(data.frame(x = est_midpoints, y = carrier), aes(x = x, y = y)) +
        geom_col(fill = "dodgerblue", alpha = 0.35, width = 0.9) +
        geom_point(color = "dodgerblue", size = 1.8) +
        labs(
          x = "Count",
          y = "Mass",
          title = "Singleton Carrier PMF"
        ) +
        theme_minimal()
    )
  }

  cellSum1 = sum(rowSums(Smat1))
  cellSum2 = sum(rowSums(Smat2))
  carrier_scale1 = carrier * cellSum1
  carrier_scale2 = carrier * cellSum2

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
      Kmax = Kmax,
      kde_used = "discrete-singleton-covfix",
      carrier_type = "singleton-pmf",
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
        binwidth = 1,
        boundary = -0.5,
        alpha = 0.20,
        position = "identity",
        color = NA
      ) +
      geom_col(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, fill = group),
        width = 0.9,
        alpha = 0.30,
        position = "identity",
        color = NA
      ) +
      geom_line(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, color = group),
        linewidth = 0.8
      ) +
      geom_point(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, color = group),
        size = 1.2
      ) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Count", y = "Mass", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }

  lambda1 = carrier * e1
  lambda2 = carrier * e2
  G_1 = t(X) %*% (lambda1 * X)
  G_2 = t(X) %*% (lambda2 * X)
  G1_inv = solve_with_fallback(G_1, ridge = ridge)
  G2_inv = solve_with_fallback(G_2, ridge = ridge)

  cov_list1 = estimate_cov_sj(Smat1)
  cov_list2 = estimate_cov_sj(Smat2)

  Z_t1_inT1 = (1 / cellSum1) * diag(K) - e1 * smooth_mat
  Z_t1_inT2 = -e1 * smooth_mat
  Z_t2_inT1 = -e2 * smooth_mat
  Z_t2_inT2 = (1 / cellSum2) * diag(K) - e2 * smooth_mat

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
    Kmax = Kmax,
    kde_used = "discrete-singleton-covfix",
    carrier_type = "singleton-pmf",
    p = p,
    binwidth = binwidth,
    bw = bw,
    carrier = carrier,
    est_midpoints = est_midpoints,
    sef_df1 = sef_df1,
    sef_df2 = sef_df2,
    plt = plt
  )
}

run_sef_test_discrete_covfix = function(sim_obj, K = NULL, bw = NULL, p = 2,
                                        verbose = TRUE, ridge = 1e-8, eps = 1e-10) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)

  pVec = setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    if (verbose) {
      message(i, " of ", length(genes), " (discrete singleton covfix)")
    }
    out = run_sef_regression_discrete_covfix(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      K = K,
      ridge = ridge,
      bw = bw,
      p = p,
      eps = eps
    )
    pVec[i] = out$pval_1
  }

  cbind(gene_info, pVec = unname(pVec[genes]))
}

run_sef_regression_discrete_covfix_scaled_basis = function(gene, counts_vec, meta,
                                                           gene_info = NULL,
                                                           ridge = 1e-8, K = NULL,
                                                           bw = NULL, p = 2,
                                                           eps = 1e-10,
                                                           plt_flag = FALSE) {
  # bw is kept only for drop-in compatibility with the KDE-based functions.
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
  y_agg = as.numeric(y_agg)

  Kmax_obs = ceiling(max(y_agg, na.rm = TRUE))
  if (is.null(K)) {
    Kmax = Kmax_obs
  } else {
    Kmax = max(as.integer(K), Kmax_obs)
  }
  Kmax = max(0L, Kmax)

  est_midpoints = seq.int(0L, Kmax)
  est_grid = seq(-0.5, Kmax + 0.5, by = 1)
  K = length(est_midpoints)
  binwidth = 1

  Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)

  # Singleton PMF carrier: one count at k affects only the mass at k.
  smooth_mat = diag(K) / pooled_total
  carrier = as.vector(smooth_mat %*% pooled_counts)
  # Zero-mass singleton bins are expected when some integer counts are unobserved.
  carrier = pmax(carrier, eps)

  x_mid = est_midpoints
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
        geom_col(fill = "dodgerblue", alpha = 0.35, width = 0.9) +
        geom_point(color = "dodgerblue", size = 1.8) +
        labs(
          x = "Count",
          y = "Mass",
          title = "Singleton Carrier PMF"
        ) +
        theme_minimal()
    )
  }

  cellSum1 = sum(rowSums(Smat1))
  cellSum2 = sum(rowSums(Smat2))
  carrier_scale1 = carrier * cellSum1
  carrier_scale2 = carrier * cellSum2

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
      Kmax = Kmax,
      kde_used = "discrete-singleton-covfix-scaled-basis",
      carrier_type = "singleton-pmf",
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
        binwidth = 1,
        boundary = -0.5,
        alpha = 0.20,
        position = "identity",
        color = NA
      ) +
      geom_col(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, fill = group),
        width = 0.9,
        alpha = 0.30,
        position = "identity",
        color = NA
      ) +
      geom_line(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, color = group),
        linewidth = 0.8
      ) +
      geom_point(
        data = df_curve,
        aes(x = est_midpoints, y = sef_value, color = group),
        size = 1.2
      ) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Count", y = "Mass", color = "Group", fill = "Group") +
      theme_minimal()

    print(plt)
  }

  lambda1 = carrier * e1
  lambda2 = carrier * e2
  G_1 = t(X) %*% (lambda1 * X)
  G_2 = t(X) %*% (lambda2 * X)
  G1_inv = solve_with_fallback(G_1, ridge = ridge)
  G2_inv = solve_with_fallback(G_2, ridge = ridge)

  cov_list1 = estimate_cov_sj(Smat1)
  cov_list2 = estimate_cov_sj(Smat2)

  Z_t1_inT1 = (1 / cellSum1) * diag(K) - e1 * smooth_mat
  Z_t1_inT2 = -e1 * smooth_mat
  Z_t2_inT1 = -e2 * smooth_mat
  Z_t2_inT2 = (1 / cellSum2) * diag(K) - e2 * smooth_mat

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
    Kmax = Kmax,
    kde_used = "discrete-singleton-covfix-scaled-basis",
    carrier_type = "singleton-pmf",
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

run_sef_test_discrete_covfix_scaled_basis = function(sim_obj, K = NULL, bw = NULL,
                                                     p = 2, verbose = TRUE,
                                                     ridge = 1e-8, eps = 1e-10) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)

  pVec = setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    if (verbose) {
      message(i, " of ", length(genes), " (discrete singleton covfix, scaled basis)")
    }
    out = run_sef_regression_discrete_covfix_scaled_basis(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      K = K,
      ridge = ridge,
      bw = bw,
      p = p,
      eps = eps
    )
    pVec[i] = out$pval_1
  }

  cbind(gene_info, pVec = unname(pVec[genes]))
}
