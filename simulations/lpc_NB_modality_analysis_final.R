# Hierarchical negative-binomial modality-only simulation.
#
# Goal:
#   Simulate donor/cell count data where modality differs between groups
#   (unimodal NB vs. bimodal NB mixture) but the model-level marginal mean
#   and marginal variance are identical after donor effects are included.
#
# Donor hierarchy:
#   S_i ~ LogNormal(-sigma_sq / 2, sigma_sq), so E[S_i] = 1.
#   Counts are then generated conditional on donor scale S_i.
#
# Moment matching:
#   Let tau_sq = Var(S_i) = exp(sigma_sq) - 1.
#   For a common target mean mu, define
#     A_mix = sum_c w_c * mu_c^2 / theta_c + sum_c w_c * (mu_c - mu)^2.
#   The bimodal mixture has
#     Var(Y | S) = S * mu + S^2 * A_mix.
#   A unimodal NB with theta_unimodal = mu^2 / A_mix has the same
#   conditional variance form:
#     Var(Y | S) = S * mu + S^2 * A_mix.
#   Therefore both groups have identical marginal moments:
#     E[Y] = mu
#     Var(Y) = mu + (1 + tau_sq) * A_mix + tau_sq * mu^2.

modal_nb_component_means = function(mu, delta, mixing_prob = 0.5) {
  stopifnot(length(mu) == 1, mu > 0)
  stopifnot(length(delta) == 1, delta > 0)
  stopifnot(length(mixing_prob) == 1, mixing_prob > 0, mixing_prob < 1)
  
  w1 = mixing_prob
  w2 = 1 - mixing_prob
  mu1 = mu - delta
  mu2 = (mu - w1 * mu1) / w2
  
  if (mu1 <= 0 || mu2 <= 0) {
    stop("Component means must be positive. Reduce delta or change mixing_prob.")
  }
  
  c(mu1, mu2)
}

modal_nb_conditional_A = function(mu, mu_components, theta_components, mixing_prob = 0.5) {
  weights = c(mixing_prob, 1 - mixing_prob)
  stopifnot(length(mu_components) == 2, all(mu_components > 0))
  stopifnot(length(theta_components) == 2, all(theta_components > 0))
  
  sum(weights * mu_components^2 / theta_components) +
    sum(weights * (mu_components - mu)^2)
}

modal_nb_equal_mv_params = function(mu = 20,
                                    delta = 8,
                                    theta_mix = 28.34304,
                                    mixing_prob = 0.5,
                                    sigma_sq = 0.05) {
  stopifnot(length(theta_mix) == 1, theta_mix > 0)
  stopifnot(length(sigma_sq) == 1, sigma_sq >= 0)
  
  mu_components = modal_nb_component_means(mu, delta, mixing_prob)
  theta_components = rep(theta_mix, 2)
  conditional_A = modal_nb_conditional_A(
    mu = mu,
    mu_components = mu_components,
    theta_components = theta_components,
    mixing_prob = mixing_prob
  )
  
  theta_unimodal = mu^2 / conditional_A
  tau_sq = expm1(sigma_sq)
  marginal_variance = mu + (1 + tau_sq) * conditional_A + tau_sq * mu^2
  
  list(
    mu = mu,
    delta = delta,
    mixing_prob = mixing_prob,
    mu_components = mu_components,
    theta_components = theta_components,
    theta_unimodal = theta_unimodal,
    conditional_A = conditional_A,
    sigma_sq = sigma_sq,
    donor_scale_variance = tau_sq,
    marginal_mean = mu,
    marginal_variance = marginal_variance
  )
}

modal_nb_equal_mv_params_from_variance = function(mu = 20,
                                                  delta = 8,
                                                  marginal_variance = 125,
                                                  mixing_prob = 0.5,
                                                  sigma_sq = 0.05) {
  stopifnot(length(mu) == 1, mu > 0)
  stopifnot(length(delta) == 1, delta > 0)
  stopifnot(length(marginal_variance) == 1, marginal_variance > mu)
  stopifnot(length(sigma_sq) == 1, sigma_sq >= 0)
  
  weights = c(mixing_prob, 1 - mixing_prob)
  mu_components = modal_nb_component_means(mu, delta, mixing_prob)
  tau_sq = expm1(sigma_sq)
  
  target_A = (marginal_variance - mu - tau_sq * mu^2) / (1 + tau_sq)
  between_components = sum(weights * (mu_components - mu)^2)
  weighted_second = sum(weights * mu_components^2)
  
  if (target_A <= 0) {
    stop("Requested marginal variance is too small for this mu and sigma_sq.")
  }
  if (target_A <= between_components) {
    stop(
      "Requested marginal variance is too small for this delta and donor-effect variance. ",
      "Reduce delta, reduce sigma_sq, or increase marginal_variance."
    )
  }
  
  theta_mix = weighted_second / (target_A - between_components)
  
  modal_nb_equal_mv_params(
    mu = mu,
    delta = delta,
    theta_mix = theta_mix,
    mixing_prob = mixing_prob,
    sigma_sq = sigma_sq
  )
}

modal_nb_draw_donor_scale = function(n, sigma_sq) {
  stopifnot(length(n) == 1, n >= 0, n %% 1 == 0)
  stopifnot(length(sigma_sq) == 1, sigma_sq >= 0)
  
  if (sigma_sq == 0) {
    return(rep(1, n))
  }
  
  rlnorm(n, meanlog = -sigma_sq / 2, sdlog = sqrt(sigma_sq))
}

modal_nb_draw_unimodal = function(k, mu, theta, donor_scale) {
  n = length(k)
  stopifnot(all(k > 0), all(k %% 1 == 0))
  stopifnot(length(mu) == 1, mu > 0)
  stopifnot(length(theta) == 1, theta > 0)
  stopifnot(length(donor_scale) == n, all(donor_scale > 0))
  
  y_list = vector("list", n)
  for (i in seq_len(n)) {
    y_list[[i]] = rnbinom(k[i], size = theta, mu = donor_scale[i] * mu)
  }
  
  y_list
}

modal_nb_draw_mixture = function(k,
                                 mu_components,
                                 theta_components,
                                 mixing_prob,
                                 donor_scale) {
  n = length(k)
  weights = c(mixing_prob, 1 - mixing_prob)
  stopifnot(all(k > 0), all(k %% 1 == 0))
  stopifnot(length(mu_components) == 2, all(mu_components > 0))
  stopifnot(length(theta_components) == 2, all(theta_components > 0))
  stopifnot(length(mixing_prob) == 1, mixing_prob > 0, mixing_prob < 1)
  stopifnot(length(donor_scale) == n, all(donor_scale > 0))
  
  y_list = vector("list", n)
  for (i in seq_len(n)) {
    component_id = sample.int(2, size = k[i], replace = TRUE, prob = weights)
    y_i = integer(k[i])
    
    for (component in 1:2) {
      idx = component_id == component
      if (any(idx)) {
        y_i[idx] = rnbinom(
          sum(idx),
          size = theta_components[component],
          mu = donor_scale[i] * mu_components[component]
        )
      }
    }
    
    y_list[[i]] = y_i
  }
  
  y_list
}

modal_nb_flatten = function(y_list) {
  as.integer(unlist(y_list, use.names = FALSE))
}

modal_nb_draw_param_table = function(n,
                                     target_mu_range = c(20, 20),
                                     delta_range = c(8, 8),
                                     target_marginal_variance = 125,
                                     theta_mix_range = c(100, 100),
                                     mixing_prob = 0.5,
                                     sigma_sq = 0.05) {
  stopifnot(length(n) == 1, n >= 0, n %% 1 == 0)
  
  if (n == 0) {
    return(data.frame())
  }
  
  mu = runif(n, min = target_mu_range[1], max = target_mu_range[2])
  delta_hi = pmin(delta_range[2], 0.95 * mu)
  delta_lo = pmin(delta_range[1], delta_hi)
  delta = runif(n, min = delta_lo, max = delta_hi)
  
  if (is.null(target_marginal_variance)) {
    theta_mix = runif(n, min = theta_mix_range[1], max = theta_mix_range[2])
    params = lapply(seq_len(n), function(g) {
      modal_nb_equal_mv_params(
        mu = mu[g],
        delta = delta[g],
        theta_mix = theta_mix[g],
        mixing_prob = mixing_prob,
        sigma_sq = sigma_sq
      )
    })
  } else {
    if (length(target_marginal_variance) == 1) {
      marginal_variance = rep(target_marginal_variance, n)
    } else if (length(target_marginal_variance) == 2) {
      marginal_variance = runif(
        n,
        min = target_marginal_variance[1],
        max = target_marginal_variance[2]
      )
    } else if (length(target_marginal_variance) == n) {
      marginal_variance = target_marginal_variance
    } else {
      stop("target_marginal_variance must be NULL, length 1, length 2, or length n.")
    }
    
    params = lapply(seq_len(n), function(g) {
      modal_nb_equal_mv_params_from_variance(
        mu = mu[g],
        delta = delta[g],
        marginal_variance = marginal_variance[g],
        mixing_prob = mixing_prob,
        sigma_sq = sigma_sq
      )
    })
  }
  
  data.frame(
    target_mean = vapply(params, `[[`, numeric(1), "mu"),
    target_marginal_variance = vapply(params, `[[`, numeric(1), "marginal_variance"),
    delta = vapply(params, `[[`, numeric(1), "delta"),
    mixing_prob = vapply(params, `[[`, numeric(1), "mixing_prob"),
    mu_component_1 = vapply(params, function(x) x$mu_components[1], numeric(1)),
    mu_component_2 = vapply(params, function(x) x$mu_components[2], numeric(1)),
    theta_component_1 = vapply(params, function(x) x$theta_components[1], numeric(1)),
    theta_component_2 = vapply(params, function(x) x$theta_components[2], numeric(1)),
    theta_unimodal = vapply(params, `[[`, numeric(1), "theta_unimodal"),
    conditional_A = vapply(params, `[[`, numeric(1), "conditional_A"),
    sigma_sq = vapply(params, `[[`, numeric(1), "sigma_sq"),
    donor_scale_variance = vapply(params, `[[`, numeric(1), "donor_scale_variance"),
    model_marginal_mean = vapply(params, `[[`, numeric(1), "marginal_mean"),
    model_marginal_variance = vapply(params, `[[`, numeric(1), "marginal_variance"),
    stringsAsFactors = FALSE
  )
}

simulate_NB_modality_analysis = function(
    seed = NULL,
    n1 = 20,
    n2 = 20,
    lower = 500,
    upper = 1000,
    k1 = NULL,
    k2 = NULL,
    n_null_unimodal = 50,
    n_null_bimodal = 50,
    n_modal = 50,
    target_mu_range = c(20, 20),
    delta_range = c(8, 8),
    target_marginal_variance = 125,
    theta_mix_range = c(100, 100),
    mixing_prob = 0.5,
    sigma_sq = 0.05,
    modal_unimodal_group = "group1",
    return_donor_effects = TRUE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  stopifnot(length(n1) == 1, n1 > 0, n1 %% 1 == 0)
  stopifnot(length(n2) == 1, n2 > 0, n2 %% 1 == 0)
  stopifnot(length(lower) == 1, length(upper) == 1, lower > 0, upper >= lower)
  stopifnot(modal_unimodal_group %in% c("group1", "group2"))
  
  if (is.null(k1)) {
    k1 = sample(lower:upper, n1, replace = TRUE)
  }
  if (is.null(k2)) {
    k2 = sample(lower:upper, n2, replace = TRUE)
  }
  stopifnot(length(k1) == n1, all(k1 > 0), all(k1 %% 1 == 0))
  stopifnot(length(k2) == n2, all(k2 > 0), all(k2 %% 1 == 0))
  
  total_genes = n_null_unimodal + n_null_bimodal + n_modal
  if (total_genes == 0) {
    stop("At least one gene must be simulated.")
  }
  
  donor_ids_1 = sprintf("G1_D%03d", seq_len(n1))
  donor_ids_2 = sprintf("G2_D%03d", seq_len(n2))
  total_cells = sum(k1) + sum(k2)
  
  metadata = data.frame(
    cell_id = sprintf("Cell_%07d", seq_len(total_cells)),
    donor_id = c(rep(donor_ids_1, times = k1), rep(donor_ids_2, times = k2)),
    donor_index_within_group = c(rep(seq_len(n1), times = k1), rep(seq_len(n2), times = k2)),
    group = c(rep("group1", sum(k1)), rep("group2", sum(k2))),
    group_indicator = c(rep(1L, sum(k1)), rep(2L, sum(k2))),
    stringsAsFactors = FALSE
  )
  
  gene_ids = c(
    sprintf("null_unimodal_%03d", seq_len(n_null_unimodal)),
    sprintf("null_bimodal_%03d", seq_len(n_null_bimodal)),
    sprintf("modal_equal_mv_%03d", seq_len(n_modal))
  )
  gene_class = c(
    rep("null_unimodal", n_null_unimodal),
    rep("null_bimodal", n_null_bimodal),
    rep("modality_only_equal_mean_variance", n_modal)
  )
  
  group1_modality = c(
    rep("unimodal", n_null_unimodal),
    rep("bimodal", n_null_bimodal),
    rep(ifelse(modal_unimodal_group == "group1", "unimodal", "bimodal"), n_modal)
  )
  group2_modality = c(
    rep("unimodal", n_null_unimodal),
    rep("bimodal", n_null_bimodal),
    rep(ifelse(modal_unimodal_group == "group2", "unimodal", "bimodal"), n_modal)
  )
  
  param_table = modal_nb_draw_param_table(
    n = total_genes,
    target_mu_range = target_mu_range,
    delta_range = delta_range,
    target_marginal_variance = target_marginal_variance,
    theta_mix_range = theta_mix_range,
    mixing_prob = mixing_prob,
    sigma_sq = sigma_sq
  )
  
  count_matrix = matrix(NA_integer_, nrow = total_cells, ncol = total_genes)
  colnames(count_matrix) = gene_ids
  rownames(count_matrix) = metadata$cell_id
  
  donor_effect_rows = vector("list", if (return_donor_effects) 2 * total_genes else 0)
  donor_effect_counter = 0L
  
  simulate_group = function(g, k, modality, donor_scale) {
    if (modality == "unimodal") {
      return(modal_nb_draw_unimodal(
        k = k,
        mu = param_table$target_mean[g],
        theta = param_table$theta_unimodal[g],
        donor_scale = donor_scale
      ))
    }
    
    modal_nb_draw_mixture(
      k = k,
      mu_components = c(param_table$mu_component_1[g], param_table$mu_component_2[g]),
      theta_components = c(param_table$theta_component_1[g], param_table$theta_component_2[g]),
      mixing_prob = param_table$mixing_prob[g],
      donor_scale = donor_scale
    )
  }
  
  record_donor_effects = function(g, group, donor_ids, donor_scale, modality) {
    data.frame(
      gene = gene_ids[g],
      gene_class = gene_class[g],
      group = group,
      donor_id = donor_ids,
      modality = modality,
      donor_scale = donor_scale,
      stringsAsFactors = FALSE
    )
  }
  
  for (g in seq_len(total_genes)) {
    donor_scale_1 = modal_nb_draw_donor_scale(n1, sigma_sq)
    donor_scale_2 = modal_nb_draw_donor_scale(n2, sigma_sq)
    
    y1 = simulate_group(g, k1, group1_modality[g], donor_scale_1)
    y2 = simulate_group(g, k2, group2_modality[g], donor_scale_2)
    count_matrix[, g] = c(modal_nb_flatten(y1), modal_nb_flatten(y2))
    
    if (return_donor_effects) {
      donor_effect_counter = donor_effect_counter + 1L
      donor_effect_rows[[donor_effect_counter]] = record_donor_effects(
        g = g,
        group = "group1",
        donor_ids = donor_ids_1,
        donor_scale = donor_scale_1,
        modality = group1_modality[g]
      )
      
      donor_effect_counter = donor_effect_counter + 1L
      donor_effect_rows[[donor_effect_counter]] = record_donor_effects(
        g = g,
        group = "group2",
        donor_ids = donor_ids_2,
        donor_scale = donor_scale_2,
        modality = group2_modality[g]
      )
    }
  }
  
  gene_info = cbind(
    data.frame(
      gene = gene_ids,
      de_status = ifelse(gene_class == "modality_only_equal_mean_variance", "DE", "nonDE"),
      effect_type = ifelse(
        gene_class == "modality_only_equal_mean_variance",
        "modality_only_equal_mean_variance",
        "none"
      ),
      gene_class = gene_class,
      group1_modality = group1_modality,
      group2_modality = group2_modality,
      stringsAsFactors = FALSE
    ),
    param_table
  )
  
  out = list(
    count_matrix = count_matrix,
    metadata = metadata,
    gene_info = gene_info,
    k1 = k1,
    k2 = k2,
    settings = list(
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
      moment_matching = "model-level equal marginal mean and variance after donor effects"
    )
  )
  
  if (return_donor_effects) {
    out$donor_effects = do.call(rbind, donor_effect_rows)
  }
  
  out
}

