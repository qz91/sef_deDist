# LPC version to use parallel computing

`%||%` = function(x, y) if (is.null(x)) y else x

# P-value combination helpers adapted from hye_MoM/main_cov_shift_1d/main_func_new.R.
combine_pvalues_bonferroni = function(p) {
  p = p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  min(1, length(p) * min(p))
}

combine_pvalues_simes = function(p) {
  p = sort(p[!is.na(p)])
  if (!length(p)) return(NA_real_)
  min(1, min(length(p) * p / seq_along(p)))
}

combine_pvalues_fisher = function(p) {
  p = p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  if (any(p <= 0)) return(0)
  stat = -2 * sum(log(p))
  stats::pchisq(stat, df = 2 * length(p), lower.tail = FALSE)
}

combine_pvalues_acat = function(p, weights = NULL, eps = 1e-15) {
  keep = !is.na(p)
  p = p[keep]
  if (!length(p)) return(NA_real_)

  if (is.null(weights)) {
    weights = rep(1 / length(p), length(p))
  } else {
    weights = as.numeric(weights)[keep]
    if (length(weights) != length(p)) {
      stop("weights must have the same length as p.")
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("ACAT weights must be finite and non-negative.")
    }
    if (sum(weights) <= 0) {
      stop("ACAT weights must have positive total weight.")
    }
    weights = weights / sum(weights)
  }

  if (any(p <= 0)) return(0)
  if (all(p >= 1)) return(1)

  p_clip = pmin(pmax(p, eps), 1 - eps)
  stat = sum(weights * tan((0.5 - p_clip) * pi))
  p_acat = 0.5 - atan(stat) / pi
  pmin(pmax(p_acat, 0), 1)
}

extract_moment_pvalue_matrix = function(
    pvals,
    p_cols = NULL,
    family = NULL,
    gene_col = "gene"
) {
  if (is.data.frame(pvals)) {
    if (is.null(p_cols)) {
      if (is.null(family)) {
        stop("Provide p_cols, or provide family to infer columns like '<family>_X_1'.")
      }
      p_cols = grep(paste0("^", family, "_X_[0-9]+$"), names(pvals), value = TRUE)
    }
    if (!length(p_cols)) {
      stop("No moment p-value columns were found.")
    }
    missing_cols = setdiff(p_cols, names(pvals))
    if (length(missing_cols)) {
      stop("p_cols not found: ", paste(missing_cols, collapse = ", "))
    }
    pmat = as.matrix(pvals[, p_cols, drop = FALSE])
    genes = if (gene_col %in% names(pvals)) {
      as.character(pvals[[gene_col]])
    } else if (!is.null(rownames(pvals))) {
      rownames(pvals)
    } else {
      as.character(seq_len(nrow(pvals)))
    }
  } else {
    pmat = as.matrix(pvals)
    if (!is.numeric(pmat)) {
      stop("pvals must be numeric, or a data frame with numeric p_cols.")
    }
    p_cols = colnames(pmat) %||% paste0("p_", seq_len(ncol(pmat)))
    genes = rownames(pmat) %||% as.character(seq_len(nrow(pmat)))
  }

  storage.mode(pmat) = "numeric"
  bad = is.finite(pmat) & (pmat < 0 | pmat > 1)
  if (any(bad)) {
    stop("All finite p-values must lie in [0, 1].")
  }
  list(pmat = pmat, genes = genes, p_cols = p_cols)
}

combine_moment_pvalues = function(
    pvals,
    p_cols = NULL,
    family = NULL,
    gene_col = "gene",
    methods = c("bonferroni", "simes", "acat", "fisher"),
    acat_weights = NULL,
    adjust_method = NULL,
    adjust_for = NULL,
    na_rm = FALSE
) {
  methods = match.arg(
    methods,
    choices = c("bonferroni", "simes", "acat", "fisher"),
    several.ok = TRUE
  )
  parsed = extract_moment_pvalue_matrix(
    pvals = pvals,
    p_cols = p_cols,
    family = family,
    gene_col = gene_col
  )
  pmat = parsed$pmat
  genes = parsed$genes
  p_cols = parsed$p_cols

  if (!is.null(acat_weights) && length(acat_weights) != ncol(pmat)) {
    stop("acat_weights must have length equal to the number of moment p-value columns.")
  }

  combine_one = function(z) {
    if (!isTRUE(na_rm) && anyNA(z)) {
      out = c(n_moments = sum(!is.na(z)), min_p = NA_real_)
      for (method in methods) out[paste0("p_", method)] = NA_real_
      return(out)
    }
    keep = !is.na(z)
    z_use = z[keep]
    out = c(
      n_moments = length(z_use),
      min_p = if (length(z_use)) min(z_use) else NA_real_
    )
    if (!length(z_use)) {
      for (method in methods) out[paste0("p_", method)] = NA_real_
      return(out)
    }
    for (method in methods) {
      out[paste0("p_", method)] = switch(
        method,
        bonferroni = combine_pvalues_bonferroni(z_use),
        simes = combine_pvalues_simes(z_use),
        acat = combine_pvalues_acat(
          z,
          weights = if (is.null(acat_weights)) NULL else acat_weights
        ),
        fisher = combine_pvalues_fisher(z_use)
      )
    }
    out
  }

  combined = t(apply(pmat, 1, combine_one))
  out = data.frame(
    gene = genes,
    n_moments = as.integer(combined[, "n_moments"]),
    min_p = as.numeric(combined[, "min_p"]),
    stringsAsFactors = FALSE
  )
  for (method in methods) {
    col = paste0("p_", method)
    out[[col]] = as.numeric(combined[, col])
  }

  if (!is.null(adjust_method)) {
    adjust_for = adjust_for %||% paste0("p_", methods[1])
    if (!(adjust_for %in% names(out))) {
      stop("adjust_for must be one of the combined p-value columns.")
    }
    out$p_adjust = p.adjust(out[[adjust_for]], method = adjust_method)
  }
  attr(out, "p_cols") = p_cols
  attr(out, "methods") = methods
  out
}

combine_moment_pvalues_bonferroni = function(
    pvals,
    p_cols = NULL,
    family = NULL,
    gene_col = "gene",
    adjust_method = NULL,
    na_rm = FALSE
) {
  out = combine_moment_pvalues(
    pvals = pvals,
    p_cols = p_cols,
    family = family,
    gene_col = gene_col,
    methods = "bonferroni",
    adjust_method = adjust_method,
    adjust_for = "p_bonferroni",
    na_rm = na_rm
  )
  names(out)[names(out) == "p_bonferroni"] = "p_bonf"
  out
}

resolve_moment_indices = function(p_total = NULL, moment_idx = NULL) {
  if (is.null(moment_idx)) {
    if (is.null(p_total)) {
      stop("Provide p_total or moment_idx.")
    }
    moment_idx = seq_len(as.integer(p_total))
  } else {
    moment_idx = as.integer(moment_idx)
  }
  stopifnot(
    "moment_idx must be positive integers" = length(moment_idx) > 0L && all(moment_idx >= 1L),
    "moment_idx must not contain duplicates" = !anyDuplicated(moment_idx)
  )
  moment_idx
}

moment_pvalue_colnames = function(moment_idx) {
  paste0("X_", as.integer(moment_idx))
}

make_pvec_df = function(gene_info, pVec) {
  genes = as.character(gene_info$gene)
  pVec = setNames(as.numeric(pVec), names(pVec))
  cbind(gene_info, pVec = unname(pVec[genes]))
}

combine_momentwise_pmat = function(
    gene_info,
    pmat,
    moment_idx = NULL,
    methods = c("bonferroni", "simes", "acat", "fisher"),
    acat_weights = NULL,
    na_rm = FALSE
) {
  genes = as.character(gene_info$gene)
  if (nrow(pmat) != length(genes)) {
    stop("pmat must have one row per gene in gene_info.")
  }
  rownames(pmat) = genes

  if (is.null(moment_idx)) {
    moment_idx = as.integer(sub("^X_", "", colnames(pmat)))
  }
  moment_idx = as.integer(moment_idx)
  p_cols = moment_pvalue_colnames(moment_idx)
  missing_cols = setdiff(p_cols, colnames(pmat))
  if (length(missing_cols)) {
    stop("Missing moment-wise p-value columns: ", paste(missing_cols, collapse = ", "))
  }

  combined = combine_moment_pvalues(
    pvals = pmat[, p_cols, drop = FALSE],
    methods = methods,
    acat_weights = acat_weights,
    na_rm = na_rm
  )

  out = lapply(methods, function(method) {
    p_col = paste0("p_", method)
    pVec = setNames(combined[[p_col]], combined$gene)
    make_pvec_df(gene_info, pVec[genes])
  })
  names(out) = methods
  attr(out, "combined_pvalues") = combined
  out
}

run_sef_regression_covfix_combination = function(gene, counts_vec, meta, gene_info = NULL,
                                                 K = NULL, bw = 2, p = 2,
                                                 p_total = NULL, p_specific = NULL,
                                                 carrier_type = c("continuous", "discrete"),
                                                 discrete_grid = c("full", "quantile"),
                                                 discrete_K = NULL,
                                                 discrete_quantile = 0.999,
                                                 scaled_basis = FALSE,
                                                 eps = 1e-10,
                                                 eps_prob = 1e-8,
                                                 ridge = NULL,
                                                 plt_flag = FALSE) {
  carrier_type = match.arg(carrier_type)
  discrete_grid = match.arg(discrete_grid)
  if (is.null(ridge)) {
    ridge = eps_prob
  }

  if (is.null(p_total)) {
    if (is.null(p_specific)) {
      p_total = p
    } else {
      p_total = max(p, max(as.integer(p_specific)))
    }
  }
  if (is.null(p_specific)) {
    moment_idx = seq_len(p_total)
  } else {
    moment_idx = as.integer(p_specific)
    stopifnot(
      "p_specific must be a positive integer vector" = all(moment_idx >= 1L),
      "max(p_specific) must be <= p_total" = max(moment_idx) <= p_total
    )
  }
  p_eff = length(moment_idx)

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

  if (carrier_type == "continuous") {
    l = min(y_agg)
    u = max(y_agg)
    if (is.null(K)) {
      K = nclass.FD(y_agg)
    }
    K = max(10, K)
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints = (est_grid[-1] + est_grid[-length(est_grid)]) / 2
    binwidth = diff(est_grid)[1]
    discrete_has_tail = FALSE
    discrete_cutoff = NA_integer_
  } else {
    y_int = as.integer(round(y_agg))
    if (any(y_int < 0, na.rm = TRUE)) {
      stop("Discrete carrier requires non-negative integer counts.")
    }
    max_y = max(y_int, na.rm = TRUE)
    if (discrete_grid == "full") {
      Kmax_obs = ceiling(max_y)
      if (is.null(K)) {
        Kmax = Kmax_obs
      } else {
        Kmax = max(as.integer(K), Kmax_obs)
      }
      Kmax = max(0L, Kmax)
      est_midpoints = as.numeric(seq.int(0L, Kmax))
      est_grid = seq(-0.5, Kmax + 0.5, by = 1)
      discrete_has_tail = FALSE
      discrete_cutoff = Kmax
    } else {
      K_arg = discrete_K
      if (is.null(K_arg) && !is.null(K)) {
        K_arg = K
      }
      if (is.null(K_arg)) {
        K_arg = as.integer(ceiling(stats::quantile(
          y_int,
          discrete_quantile,
          names = FALSE,
          na.rm = TRUE
        )))
      }
      K_arg = max(as.integer(K_arg), 1L)
      K_arg = min(K_arg, max_y)
      if (max_y <= K_arg) {
        est_grid = seq(-0.5, K_arg + 0.5, by = 1)
        est_midpoints = as.numeric(seq.int(0L, K_arg))
        discrete_has_tail = FALSE
      } else {
        est_grid = c(seq(-0.5, K_arg + 0.5, by = 1), max_y + 0.5)
        tail_mean = mean(y_int[y_int > K_arg])
        est_midpoints = c(as.numeric(seq.int(0L, K_arg)), tail_mean)
        discrete_has_tail = TRUE
      }
      discrete_cutoff = K_arg
    }
    K = length(est_midpoints)
    binwidth = 1
  }
  K_bins = length(est_midpoints)

  Smat1 = get_sc_bin_counts(Y = y1_list, bin_edges = est_grid)
  Smat2 = get_sc_bin_counts(Y = y2_list, bin_edges = est_grid)
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  pooled_counts = S1sum + S2sum
  pooled_total = sum(pooled_counts)

  if (carrier_type == "continuous") {
    kernel_mat = outer(
      est_midpoints,
      est_midpoints,
      function(a, b) dnorm(a - b, sd = bw)
    )
    smooth_mat = kernel_mat / pooled_total
  } else {
    smooth_mat = diag(K_bins) / pooled_total
  }
  carrier = as.vector(smooth_mat %*% pooled_counts)
  if (sum(carrier == 0) > 0) {
    warning("carrier_est contains ", sum(carrier == 0), " zero values. adding value eps = ", eps)
  }
  carrier = pmax(carrier, eps)

  x_mid = est_midpoints
  if (isTRUE(scaled_basis)) {
    basis_center = mean(x_mid)
    basis_scale = sd(x_mid)
    if (!is.finite(basis_scale) || basis_scale == 0) {
      basis_scale = 1
    }
    basis_midpoints = (x_mid - basis_center) / basis_scale
    basis_type = "scaled"
  } else {
    basis_center = 0
    basis_scale = 1
    basis_midpoints = x_mid
    basis_type = "raw"
  }

  X = cbind(
    Intercept = 1,
    vapply(moment_idx, function(d) basis_midpoints^d, numeric(K_bins))
  )
  colnames(X) = c("Intercept", paste0("X_", moment_idx))

  if (plt_flag) {
    carrier_data = data.frame(x = est_midpoints, y = carrier)
    if (carrier_type == "continuous") {
      print(
        ggplot(carrier_data, aes(x = x, y = y)) +
          geom_area(fill = "dodgerblue", alpha = 0.25) +
          geom_line(color = "dodgerblue", linewidth = 1.2) +
          labs(x = "Expression", y = "Density", title = "Carrier Density") +
          theme_minimal()
      )
    } else {
      print(
        ggplot(carrier_data, aes(x = x, y = y)) +
          geom_col(fill = "dodgerblue", alpha = 0.35, width = 0.9) +
          geom_point(color = "dodgerblue", size = 1.8) +
          labs(x = "Count", y = "Mass", title = "Singleton Carrier PMF") +
          theme_minimal()
      )
    }
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
    nm = paste0("X_", moment_idx)
    na_vec = setNames(rep(NA_real_, p_eff), nm)
    list(
      gene = gene,
      pval_1 = NA_real_,
      pval_omnibus = NA_real_,
      pval_marginal_z = na_vec,
      pval_marginal_chi = na_vec,
      pval_decor_z = na_vec,
      pval_decor_chi = na_vec,
      z_marginal = na_vec,
      z_decorrelated = na_vec,
      beta_est1 = NULL,
      beta_est2 = NULL,
      beta_diff = NULL,
      Cov_bar = NULL,
      Cov1 = NULL,
      Cov2 = NULL,
      K = K,
      kde_used = if (carrier_type == "continuous") {
        "gaussian-kernel-matrix-covfix-specific"
      } else {
        "discrete-singleton-covfix-specific"
      },
      carrier_type = carrier_type,
      discrete_grid = if (carrier_type == "discrete") discrete_grid else NA_character_,
      p = p,
      p_total = p_total,
      p_specific = p_specific,
      moment_idx = moment_idx,
      basis_type = basis_type,
      binwidth = binwidth,
      bw = bw,
      ridge = ridge,
      eps_prob = eps_prob
    )
  }

  fit1 = tryCatch(
    suppressWarnings(glm(formula, family = poisson(link = "log"), data = df1, control = ctrl)),
    error = function(e) NULL
  )
  fit2 = tryCatch(
    suppressWarnings(glm(formula, family = poisson(link = "log"), data = df2, control = ctrl)),
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
      group = factor(rep(c("Group 1", "Group 2"), each = K_bins))
    )
    df_hist = data.frame(
      value = y_agg,
      group = factor(c(rep("Group 1", cellSum1), rep("Group 2", cellSum2)))
    )
    plt = ggplot() +
      geom_histogram(
        data = df_hist,
        aes(x = value, y = after_stat(density), fill = group),
        bins = if (carrier_type == "continuous") 50 else NULL,
        binwidth = if (carrier_type == "discrete") 1 else NULL,
        boundary = if (carrier_type == "discrete") -0.5 else NULL,
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
      labs(x = if (carrier_type == "continuous") "Expression" else "Count",
           y = if (carrier_type == "continuous") "Density" else "Mass",
           color = "Group",
           fill = "Group") +
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

  Z_t1_inT1 = (1 / cellSum1) * diag(K_bins) - (e1 * binwidth) * smooth_mat
  Z_t1_inT2 = -(e1 * binwidth) * smooth_mat
  Z_t2_inT1 = -(e2 * binwidth) * smooth_mat
  Z_t2_inT2 = (1 / cellSum2) * diag(K_bins) - (e2 * binwidth) * smooth_mat

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
  beta_diff_no_int = beta_diff[-1]
  nm = paste0("X_", moment_idx)
  names(beta_diff_no_int) = nm

  chi_stat1 = as.numeric(
    beta_diff_no_int %*% solve_with_fallback(stat_mat, beta_diff_no_int, ridge = ridge)
  )
  pval_omnibus = pchisq(chi_stat1, df = p_eff, lower.tail = FALSE)

  diag_var = pmax(diag(stat_mat), 0)
  se_marginal = sqrt(diag_var)
  z_marginal = beta_diff_no_int / se_marginal
  pval_marginal_z = 2 * pnorm(-abs(z_marginal))
  pval_marginal_chi = pchisq(z_marginal^2, df = 1, lower.tail = FALSE)

  inv_sqrt = inv_sqrt_psd(stat_mat, ridge = ridge)
  z_decorrelated = as.vector(inv_sqrt %*% beta_diff_no_int)
  names(z_decorrelated) = nm
  pval_decor_z = 2 * pnorm(-abs(z_decorrelated))
  pval_decor_chi = pchisq(z_decorrelated^2, df = 1, lower.tail = FALSE)

  list(
    gene = gene,
    pval_1 = pval_omnibus,
    pval_omnibus = pval_omnibus,
    pval_marginal_z = pval_marginal_z,
    pval_marginal_chi = pval_marginal_chi,
    pval_decor_z = pval_decor_z,
    pval_decor_chi = pval_decor_chi,
    z_marginal = z_marginal,
    z_decorrelated = z_decorrelated,
    chi_stat1 = chi_stat1,
    beta_est1 = beta_est1,
    beta_est2 = beta_est2,
    beta_diff = beta_diff,
    Cov_bar = Cov_bar,
    Cov1 = Cov1,
    Cov2 = Cov2,
    K = K,
    K_bins = K_bins,
    kde_used = if (carrier_type == "continuous") {
      "gaussian-kernel-matrix-covfix-specific"
    } else {
      "discrete-singleton-covfix-specific"
    },
    carrier_type = carrier_type,
    discrete_grid = if (carrier_type == "discrete") discrete_grid else NA_character_,
    discrete_cutoff = discrete_cutoff,
    discrete_has_tail = discrete_has_tail,
    p = p,
    p_total = p_total,
    p_specific = p_specific,
    moment_idx = moment_idx,
    basis_type = basis_type,
    basis_center = basis_center,
    basis_scale = basis_scale,
    binwidth = binwidth,
    bw = bw,
    eps = eps,
    eps_prob = eps_prob,
    ridge = ridge,
    carrier = carrier,
    est_midpoints = est_midpoints,
    basis_midpoints = basis_midpoints,
    sef_df1 = sef_df1,
    sef_df2 = sef_df2,
    plt = plt
  )
}

run_sef_test_covfix_combination = function(sim_obj, K = NULL, bw = 0.5, p = 2,
                                           p_total = NULL, p_specific = NULL,
                                           carrier_type = c("continuous", "discrete"),
                                           discrete_grid = c("full", "quantile"),
                                           discrete_K = NULL,
                                           discrete_quantile = 0.999,
                                           scaled_basis = FALSE,
                                           eps = 1e-10,
                                           eps_prob = 1e-8,
                                           ridge = NULL,
                                           verbose = TRUE) {
  carrier_type = match.arg(carrier_type)
  discrete_grid = match.arg(discrete_grid)

  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)

  pVec = setNames(rep(NA_real_, length(genes)), genes)
  for (i in seq_along(genes)) {
    if (verbose) {
      message(i, " of ", length(genes), " (covfix combination)")
    }
    out = run_sef_regression_covfix_combination(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      K = K,
      bw = bw,
      p = p,
      p_total = p_total,
      p_specific = p_specific,
      carrier_type = carrier_type,
      discrete_grid = discrete_grid,
      discrete_K = discrete_K,
      discrete_quantile = discrete_quantile,
      scaled_basis = scaled_basis,
      eps = eps,
      eps_prob = eps_prob,
      ridge = ridge
    )
    pVec[i] = out$pval_1
  }

  cbind(gene_info, pVec = unname(pVec[genes]))
}

run_sef_test_covfix_momentwise = function(sim_obj, K = NULL, bw = 0.5,
                                          p_total = 4,
                                          moment_idx = NULL,
                                          carrier_type = c("continuous", "discrete"),
                                          discrete_grid = c("full", "quantile"),
                                          discrete_K = NULL,
                                          discrete_quantile = 0.999,
                                          scaled_basis = FALSE,
                                          eps = 1e-10,
                                          eps_prob = 1e-8,
                                          ridge = NULL,
                                          verbose = TRUE) {
  carrier_type = match.arg(carrier_type)
  discrete_grid = match.arg(discrete_grid)
  moment_idx = resolve_moment_indices(p_total = p_total, moment_idx = moment_idx)
  p_total = max(moment_idx)

  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)
  p_cols = moment_pvalue_colnames(moment_idx)
  pmat = matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(moment_idx),
    dimnames = list(genes, p_cols)
  )
  moment_fits = vector("list", length(moment_idx))
  names(moment_fits) = p_cols

  for (j in seq_along(moment_idx)) {
    d = moment_idx[j]
    if (verbose) {
      message("SEF moment-wise test for moment ", d)
    }
    fit = run_sef_test_covfix_combination(
      sim_obj = sim_obj,
      K = K,
      bw = bw,
      p = p_total,
      p_total = p_total,
      p_specific = d,
      carrier_type = carrier_type,
      discrete_grid = discrete_grid,
      discrete_K = discrete_K,
      discrete_quantile = discrete_quantile,
      scaled_basis = scaled_basis,
      eps = eps,
      eps_prob = eps_prob,
      ridge = ridge,
      verbose = FALSE
    )
    fit_p = setNames(as.numeric(fit$pVec), as.character(fit$gene))
    pmat[genes, j] = fit_p[genes]
    moment_fits[[j]] = fit
  }

  out = cbind(gene_info, as.data.frame(pmat, stringsAsFactors = FALSE))
  attr(out, "moment_fits") = moment_fits
  attr(out, "moment_idx") = moment_idx
  attr(out, "p_total") = p_total
  out
}

run_sef_test_covfix_pcombo = function(sim_obj, K = NULL, bw = 0.5,
                                      p_total = 4,
                                      moment_idx = NULL,
                                      methods = c("bonferroni", "simes", "acat", "fisher"),
                                      acat_weights = NULL,
                                      na_rm = FALSE,
                                      carrier_type = c("continuous", "discrete"),
                                      discrete_grid = c("full", "quantile"),
                                      discrete_K = NULL,
                                      discrete_quantile = 0.999,
                                      scaled_basis = FALSE,
                                      eps = 1e-10,
                                      eps_prob = 1e-8,
                                      ridge = NULL,
                                      verbose = TRUE) {
  momentwise = run_sef_test_covfix_momentwise(
    sim_obj = sim_obj,
    K = K,
    bw = bw,
    p_total = p_total,
    moment_idx = moment_idx,
    carrier_type = carrier_type,
    discrete_grid = discrete_grid,
    discrete_K = discrete_K,
    discrete_quantile = discrete_quantile,
    scaled_basis = scaled_basis,
    eps = eps,
    eps_prob = eps_prob,
    ridge = ridge,
    verbose = verbose
  )
  moment_idx = attr(momentwise, "moment_idx")
  p_cols = moment_pvalue_colnames(moment_idx)
  pmat = as.matrix(momentwise[, p_cols, drop = FALSE])

  combined = combine_momentwise_pmat(
    gene_info = sim_obj$gene_info,
    pmat = pmat,
    moment_idx = moment_idx,
    methods = methods,
    acat_weights = acat_weights,
    na_rm = na_rm
  )
  attr(combined, "momentwise") = momentwise
  combined
}

run_hotelling_mom_regression_covfix_combination = function(gene, counts_vec, meta,
                                                           gene_info = NULL,
                                                           p = 4,
                                                           p_total = NULL,
                                                           p_specific = NULL,
                                                           eps_prob = 1e-8,
                                                           ridge = NULL) {
  if (is.null(ridge)) {
    ridge = eps_prob
  }
  if (is.null(p_total)) {
    if (is.null(p_specific)) {
      p_total = p
    } else {
      p_total = max(p, max(as.integer(p_specific)))
    }
  }
  if (is.null(p_specific)) {
    moment_idx = seq_len(p_total)
  } else {
    moment_idx = as.integer(p_specific)
    stopifnot(
      "p_specific must be a positive integer vector" = all(moment_idx >= 1L),
      "max(p_specific) must be <= p_total" = max(moment_idx) <= p_total
    )
  }
  p_eff = length(moment_idx)
  nm = paste0("M_", moment_idx)

  as_cov_matrix_local = function(x, nm) {
    out = stats::cov(as.matrix(x))
    if (is.null(dim(out))) {
      out = matrix(out, nrow = 1L, ncol = 1L)
    }
    dimnames(out) = list(nm, nm)
    out
  }

  get_donor_moment_average_matrix_local = function(y_list, moment_idx) {
    out = do.call(
      rbind,
      lapply(y_list, function(y) colMeans(get_moment_matrix(y, moment_idx)))
    )
    colnames(out) = paste0("M_", moment_idx)
    rownames(out) = names(y_list)
    out
  }

  fail_return = function(msg) {
    warning("Gene ", gene, ": ", msg)
    na_vec = setNames(rep(NA_real_, p_eff), nm)
    list(
      gene = gene,
      pval_1 = NA_real_,
      pval_omnibus = NA_real_,
      hotelling_t2 = NA_real_,
      f_stat = NA_real_,
      df1 = p_eff,
      df2 = NA_real_,
      component_t = na_vec,
      component_t2 = na_vec,
      component_pval = na_vec,
      component_df = NA_real_,
      diff_mean = na_vec,
      mean1 = na_vec,
      mean2 = na_vec,
      S1 = matrix(NA_real_, p_eff, p_eff, dimnames = list(nm, nm)),
      S2 = matrix(NA_real_, p_eff, p_eff, dimnames = list(nm, nm)),
      S_pooled = matrix(NA_real_, p_eff, p_eff, dimnames = list(nm, nm)),
      V_diff = matrix(NA_real_, p_eff, p_eff, dimnames = list(nm, nm)),
      p = p,
      p_total = p_total,
      p_specific = p_specific,
      moment_idx = moment_idx,
      ridge = ridge,
      eps_prob = eps_prob,
      reference_distribution = "hotelling_f_unweighted"
    )
  }

  cell_names = names(counts_vec)
  meta_g = meta[meta$cell_id %in% cell_names, , drop = FALSE]
  meta_g = meta_g[match(cell_names, meta_g$cell_id), , drop = FALSE]

  g1 = meta_g$group == "group1"
  g2 = meta_g$group == "group2"
  y1_list = split(counts_vec[meta_g$cell_id[g1]], meta_g$donor_id[g1], drop = TRUE)
  y2_list = split(counts_vec[meta_g$cell_id[g2]], meta_g$donor_id[g2], drop = TRUE)

  n1 = length(y1_list)
  n2 = length(y2_list)
  if (n1 < 2L || n2 < 2L) {
    return(fail_return("Hotelling test needs at least two donors per group."))
  }
  if (n1 + n2 <= p_eff + 1L) {
    return(fail_return("Hotelling test residual degrees of freedom are non-positive."))
  }

  avg1 = get_donor_moment_average_matrix_local(y1_list, moment_idx)
  avg2 = get_donor_moment_average_matrix_local(y2_list, moment_idx)
  mean1 = colMeans(avg1)
  mean2 = colMeans(avg2)
  diff_mean = mean1 - mean2
  names(mean1) = nm
  names(mean2) = nm
  names(diff_mean) = nm

  S1 = as_cov_matrix_local(avg1, nm)
  S2 = as_cov_matrix_local(avg2, nm)
  S_pooled = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  S_pooled = 0.5 * (S_pooled + t(S_pooled))
  dimnames(S_pooled) = list(nm, nm)

  hotelling_t2 = as.numeric(
    (n1 * n2) / (n1 + n2) *
      diff_mean %*% solve_with_fallback(S_pooled, diff_mean, ridge = ridge)
  )
  df1 = p_eff
  df2 = n1 + n2 - p_eff - 1L
  f_stat = (df2 / (df1 * (n1 + n2 - 2))) * hotelling_t2
  pval_omnibus = stats::pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)

  component_df = n1 + n2 - 2L
  component_se = sqrt(diag(S_pooled) * (1 / n1 + 1 / n2))
  component_t = diff_mean / component_se
  component_t2 = component_t^2
  component_pval = stats::pf(component_t2, df1 = 1, df2 = component_df, lower.tail = FALSE)
  names(component_t) = nm
  names(component_t2) = nm
  names(component_pval) = nm

  V_diff = S1 / n1 + S2 / n2
  V_diff = 0.5 * (V_diff + t(V_diff))
  dimnames(V_diff) = list(nm, nm)

  list(
    gene = gene,
    pval_1 = pval_omnibus,
    pval_omnibus = pval_omnibus,
    hotelling_t2 = hotelling_t2,
    f_stat = f_stat,
    df1 = df1,
    df2 = df2,
    component_t = component_t,
    component_t2 = component_t2,
    component_pval = component_pval,
    component_df = component_df,
    diff_mean = diff_mean,
    mean1 = mean1,
    mean2 = mean2,
    donor_avg_group1 = avg1,
    donor_avg_group2 = avg2,
    S1 = S1,
    S2 = S2,
    S_pooled = S_pooled,
    V_diff = V_diff,
    p = p,
    p_total = p_total,
    p_specific = p_specific,
    moment_idx = moment_idx,
    n1 = n1,
    n2 = n2,
    ridge = ridge,
    eps_prob = eps_prob,
    reference_distribution = "hotelling_f_unweighted"
  )
}

run_hotelling_mom_test_covfix_combination = function(sim_obj,
                                                     p = 4,
                                                     p_total = NULL,
                                                     p_specific = NULL,
                                                     eps_prob = 1e-8,
                                                     ridge = NULL,
                                                     verbose = TRUE,
                                                     parallel = FALSE,
                                                     n_cores = 1L) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)

  fit_one_gene = function(i) {
    if (verbose) {
      message(i, " of ", length(genes), " (Hotelling covfix combination)")
    }
    out = run_hotelling_mom_regression_covfix_combination(
      gene = genes[i],
      counts_vec = counts[, genes[i]],
      meta = meta,
      gene_info = gene_info,
      p = p,
      p_total = p_total,
      p_specific = p_specific,
      eps_prob = eps_prob,
      ridge = ridge
    )
    list(pval = out$pval_omnibus, detail = out)
  }

  n_cores_int = suppressWarnings(as.integer(n_cores))
  if (is.na(n_cores_int) || n_cores_int < 1L) {
    n_cores_int = 1L
  }
  use_parallel = isTRUE(parallel) && length(genes) > 1L && n_cores_int > 1L
  if (use_parallel) {
    out_list = parallel::mclapply(seq_along(genes), fit_one_gene, mc.cores = n_cores_int)
  } else {
    out_list = lapply(seq_along(genes), fit_one_gene)
  }

  pVec = setNames(vapply(out_list, function(z) z$pval, numeric(1)), genes)
  detail_list = setNames(lapply(out_list, function(z) z$detail), genes)
  out_df = cbind(gene_info, pVec = unname(pVec[genes]))
  attr(out_df, "details") = detail_list
  out_df
}

run_hotelling_mom_test_covfix_componentwise = function(sim_obj,
                                                       p = 4,
                                                       p_total = NULL,
                                                       moment_idx = NULL,
                                                       eps_prob = 1e-8,
                                                       ridge = NULL,
                                                       verbose = TRUE,
                                                       parallel = FALSE,
                                                       n_cores = 1L) {
  if (is.null(p_total)) {
    p_total = if (is.null(moment_idx)) p else max(as.integer(moment_idx))
  }
  moment_idx = resolve_moment_indices(p_total = p_total, moment_idx = moment_idx)
  p_total = max(moment_idx)

  test = run_hotelling_mom_test_covfix_combination(
    sim_obj = sim_obj,
    p = p_total,
    p_total = p_total,
    p_specific = moment_idx,
    eps_prob = eps_prob,
    ridge = ridge,
    verbose = verbose,
    parallel = parallel,
    n_cores = n_cores
  )
  details = attr(test, "details")
  if (is.null(details)) {
    stop("Hotelling component p-values require detail output.")
  }

  gene_info = sim_obj$gene_info
  genes = as.character(gene_info$gene)
  p_cols = moment_pvalue_colnames(moment_idx)
  pmat = matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(moment_idx),
    dimnames = list(genes, p_cols)
  )

  for (gene in genes) {
    detail = details[[gene]]
    if (is.null(detail) || is.null(detail$component_pval)) next
    component_names = paste0("M_", moment_idx)
    component_p = setNames(as.numeric(detail$component_pval), names(detail$component_pval))
    pmat[gene, ] = component_p[component_names]
  }

  out = cbind(gene_info, as.data.frame(pmat, stringsAsFactors = FALSE))
  attr(out, "omnibus_test") = test
  attr(out, "details") = details
  attr(out, "moment_idx") = moment_idx
  attr(out, "p_total") = p_total
  out
}

run_hotelling_mom_test_covfix_pcombo = function(sim_obj,
                                                p = 4,
                                                p_total = NULL,
                                                moment_idx = NULL,
                                                methods = c("bonferroni", "simes", "acat", "fisher"),
                                                acat_weights = NULL,
                                                na_rm = FALSE,
                                                eps_prob = 1e-8,
                                                ridge = NULL,
                                                verbose = TRUE,
                                                parallel = FALSE,
                                                n_cores = 1L) {
  componentwise = run_hotelling_mom_test_covfix_componentwise(
    sim_obj = sim_obj,
    p = p,
    p_total = p_total,
    moment_idx = moment_idx,
    eps_prob = eps_prob,
    ridge = ridge,
    verbose = verbose,
    parallel = parallel,
    n_cores = n_cores
  )
  moment_idx = attr(componentwise, "moment_idx")
  p_cols = moment_pvalue_colnames(moment_idx)
  pmat = as.matrix(componentwise[, p_cols, drop = FALSE])

  combined = combine_momentwise_pmat(
    gene_info = sim_obj$gene_info,
    pmat = pmat,
    moment_idx = moment_idx,
    methods = methods,
    acat_weights = acat_weights,
    na_rm = na_rm
  )
  attr(combined, "componentwise") = componentwise
  attr(combined, "omnibus_test") = attr(componentwise, "omnibus_test")
  combined
}
