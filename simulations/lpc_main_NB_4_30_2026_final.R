# member functions to help generate NB data used in simulation revisions: FDR test

generate_NB_fdr = function(n, k, mu = 10, theta = 5, sigma_sq = 0.5) {
  # simulates 1 gene for n donors that follow unimodal NB distribution
  # start from a global mean mu
  # for each subject/donor i, draw a subject-specific mean: mu_i ~ lognormal(log(mu) - sigma_sq/2, sigma_sq)
  # generate that subject’s cell counts from Y_ij | mu_i ~ NB(mean = mu_i, size = theta)
  stopifnot(length(k) == n, all(k > 0), all(k %% 1 == 0))
  stopifnot(length(mu) == 1, mu > 0, length(theta) == 1, theta > 0, sigma_sq >= 0)
  
  # E[mu_i] is preserved at mu under the log-normal subject effect.
  log_mean = log(mu) - sigma_sq / 2
  mu_i = rlnorm(n, meanlog = log_mean, sdlog = sqrt(sigma_sq))
  
  y_list = vector("list", n)
  for (i in seq_len(n)) {
    y_list[[i]] = rnbinom(k[i], size = theta, mu = mu_i[i])
  }
  
  y_list
}


generate_NB_mixture_fdr = function(n,
                                   k,
                                   mu = c(8, 16),
                                   theta = c(5, 5),
                                   mixing_prob = 0.5,
                                   sigma_sq = 0.5) {
  # simulates 1 gene for n donors that follow mixture of NB distribution
  stopifnot(length(k) == n, all(k > 0), all(k %% 1 == 0))
  stopifnot(length(mu) == 2, all(mu > 0))
  stopifnot(length(theta) %in% c(1, 2), all(theta > 0))
  stopifnot(length(mixing_prob) == 1, mixing_prob > 0, mixing_prob < 1)
  stopifnot(sigma_sq >= 0)
  
  if (length(theta) == 1) {
    theta = rep(theta, 2)
  }
  
  subject_scale = rlnorm(n, meanlog = -sigma_sq / 2, sdlog = sqrt(sigma_sq))
  
  y_list = vector("list", n)
  for (i in seq_len(n)) {
    mu_i = mu * subject_scale[i]
    component_id = rbinom(k[i], size = 1, prob = mixing_prob) + 1L
    y_i = integer(k[i])
    
    idx1 = component_id == 1L
    idx2 = !idx1
    
    if (any(idx1)) {
      y_i[idx1] = rnbinom(sum(idx1), size = theta[1], mu = mu_i[1])
    }
    if (any(idx2)) {
      y_i[idx2] = rnbinom(sum(idx2), size = theta[2], mu = mu_i[2])
    }
    
    y_list[[i]] = y_i
  }
  
  y_list
}


.draw_truncated_lognormal = function(n,
                                     lower,
                                     upper,
                                     central_prob = c(0.05, 0.95)) {
  # helper function to simulate positive parameters: mean (mu), dispersion (theta), fold change from LN distribution
  # lower and upper provide bounds so it is truncated to prevent overly extreme values
  stopifnot(length(n) == 1, n >= 0, n %% 1 == 0)
  stopifnot(length(lower) == 1, length(upper) == 1, lower > 0, upper > lower)
  stopifnot(length(central_prob) == 2, all(central_prob > 0), all(central_prob < 1))
  stopifnot(central_prob[1] < central_prob[2])
  
  if (n == 0) {
    return(numeric(0))
  }
  
  z = qnorm(central_prob)
  meanlog = (log(lower) + log(upper)) / 2
  sdlog = (log(upper) - log(lower)) / (z[2] - z[1])
  
  lower_cdf = plnorm(lower, meanlog = meanlog, sdlog = sdlog)
  upper_cdf = plnorm(upper, meanlog = meanlog, sdlog = sdlog)
  qlnorm(runif(n, min = lower_cdf, max = upper_cdf), meanlog = meanlog, sdlog = sdlog)
}


.flatten_subject_draws = function(y_list) {
  # helper that converts list of per-donor count vectors into a long numeric vector 
  as.numeric(unlist(y_list, use.names = FALSE))
}



generate_NB_gene_data = function(
    n1, n2,
    lower = 500, upper = 1000,
    k1 = NULL, k2 = NULL,
    n_nonDE_unimodal = 600,
    n_nonDE_mixture = 600,
    n_DE = 200,
    sigma_sq_nonDE = 0.1,
    sigma_sq_DE = sigma_sq_nonDE,
    sigma_sq = NULL,
    nonDE_unimodal_mu_range = c(5, 20),
    nonDE_unimodal_theta_range = c(3, 10),
    nonDE_mixture_mu_range = c(5, 15),
    nonDE_mixture_theta_range = c(3, 10),
    nonDE_mixture_fc_range = c(2, 4),
    nonDE_mixture_prob_range = c(0.30, 0.70),
    DE_mu_range = c(5, 15),
    DE_theta_range = c(3, 10),
    DE_mixture_fc_range = c(2, 4),
    DE_mixture_prob_range = c(0.30, 0.70),
    de_unimodal_group = "group1",
    share_mode = TRUE,
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.null(sigma_sq)) {
    sigma_sq_nonDE = sigma_sq
    sigma_sq_DE = sigma_sq
  }
  
  stopifnot(length(lower) == 1, length(upper) == 1, lower > 0, upper >= lower)
  stopifnot(de_unimodal_group %in% c("group1", "group2"))
  stopifnot(sigma_sq_nonDE >= 0, sigma_sq_DE >= 0)
  stopifnot(all(nonDE_unimodal_mu_range > 0), all(nonDE_unimodal_theta_range > 0))
  stopifnot(all(nonDE_mixture_mu_range > 0), all(nonDE_mixture_theta_range > 0))
  stopifnot(all(DE_mu_range > 0), all(DE_theta_range > 0))
  
  if (is.null(k1)) {
    k1 = sample(lower:upper, n1, replace = TRUE)
  }
  if (is.null(k2)) {
    k2 = sample(lower:upper, n2, replace = TRUE)
  }
  
  stopifnot(length(k1) == n1, length(k2) == n2)
  stopifnot(all(k1 > 0), all(k2 > 0), all(k1 %% 1 == 0), all(k2 %% 1 == 0))
  
  donor_ids_1 = sprintf("G1_D%03d", seq_len(n1)) # donor names for group 1
  donor_ids_2 = sprintf("G2_D%03d", seq_len(n2)) # donor names for group 2
  total_cells = sum(k1) + sum(k2) # total cells
  
  metadata = data.frame( # make metadata
    cell_id = sprintf("Cell_%07d", seq_len(total_cells)), # cell IDs
    donor_id = c(rep(donor_ids_1, times = k1), rep(donor_ids_2, times = k2)), # donor assignment
    donor_index_within_group = c(rep(seq_len(n1), times = k1), rep(seq_len(n2), times = k2)), # indices for donor
    group = c(rep("group1", sum(k1)), rep("group2", sum(k2))), # group assignment
    group_indicator = c(rep(1L, sum(k1)), rep(2L, sum(k2))), # indicator for group assignment
    stringsAsFactors = FALSE
  )
  
  total_genes = n_nonDE_unimodal + n_nonDE_mixture + n_DE
  gene_ids = c( # set gene names
    sprintf("nonDE_unimodal_%03d", seq_len(n_nonDE_unimodal)),
    sprintf("nonDE_mixture_%03d", seq_len(n_nonDE_mixture)),
    sprintf("DE_%03d", seq_len(n_DE))
  )
  
  simulate_block = function(n_genes, draw_fun, prefix) {
    # helper function to help build count matrix based on gene type (nonDE unimodal, nonDE mixture, DE mixture)
    if (n_genes <= 0) {
      return(matrix(numeric(0), nrow = total_cells, ncol = 0))
    }
    
    out = vapply(seq_len(n_genes), draw_fun, numeric(total_cells))
    colnames(out) = sprintf("%s_%03d", prefix, seq_len(n_genes))
    out
  }
  
  mu_nonDE_unimodal = .draw_truncated_lognormal( 
    # get mean parameters for nonDE unimodal
    n_nonDE_unimodal,
    lower = nonDE_unimodal_mu_range[1],
    upper = nonDE_unimodal_mu_range[2]
  )
  theta_nonDE_unimodal = .draw_truncated_lognormal(
    # get dispersion parameters for nonDE unimodal
    n_nonDE_unimodal,
    lower = nonDE_unimodal_theta_range[1],
    upper = nonDE_unimodal_theta_range[2]
  )
  
  base_mu_nonDE_mixture = .draw_truncated_lognormal(
    # get dispersion parameters for nonDE mixture
    n_nonDE_mixture,
    lower = nonDE_mixture_mu_range[1],
    upper = nonDE_mixture_mu_range[2]
  )
  theta_nonDE_mixture = .draw_truncated_lognormal(
    # get dispersion parameters for nonDE mixture
    n_nonDE_mixture,
    lower = nonDE_mixture_theta_range[1],
    upper = nonDE_mixture_theta_range[2]
  )
  fc_nonDE_mixture = .draw_truncated_lognormal(
    # get FC value for nonDE mixture
    n_nonDE_mixture,
    lower = nonDE_mixture_fc_range[1],
    upper = nonDE_mixture_fc_range[2]
  )
  mix_prob_nonDE_mixture = runif(
    # bernoulli probability for nonDE mixture
    n_nonDE_mixture,
    min = nonDE_mixture_prob_range[1],
    max = nonDE_mixture_prob_range[2]
  )
  
  base_mu_DE = .draw_truncated_lognormal(
    # get mean parameters for DE mixture
    n_DE,
    lower = DE_mu_range[1],
    upper = DE_mu_range[2]
  )
  theta_DE = .draw_truncated_lognormal(
    # get dispersion parameters for DE mixture
    n_DE,
    lower = DE_theta_range[1],
    upper = DE_theta_range[2]
  )
  fc_DE = .draw_truncated_lognormal(
    # get FC for DE mixture
    n_DE,
    lower = DE_mixture_fc_range[1],
    upper = DE_mixture_fc_range[2]
  )
  mix_prob_DE = runif(
    # get mixture probs
    n_DE,
    min = DE_mixture_prob_range[1],
    max = DE_mixture_prob_range[2]
  )
  
  mu1_nonDE_mixture = base_mu_nonDE_mixture
  mu2_nonDE_mixture = base_mu_nonDE_mixture * fc_nonDE_mixture
  mu1_DE = base_mu_DE
  mu2_DE = base_mu_DE * fc_DE
  mu_unimodal_DE = if (share_mode) mu1_DE else (mu1_DE + mu2_DE) / 2
  
  nonDE_unimodal_mat = simulate_block(
    n_nonDE_unimodal,
    function(g) {
      c(
        .flatten_subject_draws(generate_NB_fdr(n1, k1, mu_nonDE_unimodal[g], theta_nonDE_unimodal[g], sigma_sq_nonDE)),
        .flatten_subject_draws(generate_NB_fdr(n2, k2, mu_nonDE_unimodal[g], theta_nonDE_unimodal[g], sigma_sq_nonDE))
      )
    },
    "nonDE_unimodal"
  )
  
  nonDE_mixture_mat = simulate_block(
    n_nonDE_mixture,
    function(g) {
      mu_vec = c(mu1_nonDE_mixture[g], mu2_nonDE_mixture[g])
      theta_vec = c(theta_nonDE_mixture[g], theta_nonDE_mixture[g])
      c(
        .flatten_subject_draws(generate_NB_mixture_fdr(n1, k1, mu_vec, theta_vec, mix_prob_nonDE_mixture[g], sigma_sq_nonDE)),
        .flatten_subject_draws(generate_NB_mixture_fdr(n2, k2, mu_vec, theta_vec, mix_prob_nonDE_mixture[g], sigma_sq_nonDE))
      )
    },
    "nonDE_mixture"
  )
  
  DE_mat = simulate_block(
    n_DE,
    function(g) {
      mixture_mu = c(mu1_DE[g], mu2_DE[g])
      theta_vec = c(theta_DE[g], theta_DE[g])
      
      if (de_unimodal_group == "group1") {
        c(
          .flatten_subject_draws(generate_NB_fdr(n1, k1, mu_unimodal_DE[g], theta_DE[g], sigma_sq_DE)),
          .flatten_subject_draws(generate_NB_mixture_fdr(n2, k2, mixture_mu, theta_vec, mix_prob_DE[g], sigma_sq_DE))
        )
      } else {
        c(
          .flatten_subject_draws(generate_NB_mixture_fdr(n1, k1, mixture_mu, theta_vec, mix_prob_DE[g], sigma_sq_DE)),
          .flatten_subject_draws(generate_NB_fdr(n2, k2, mu_unimodal_DE[g], theta_DE[g], sigma_sq_DE))
        )
      }
    },
    "DE"
  )
  
  count_matrix = cbind(nonDE_unimodal_mat, nonDE_mixture_mat, DE_mat)
  rownames(count_matrix) = metadata$cell_id
  
  gene_info = data.frame(
    gene = gene_ids,
    modality = c(
      rep("unimodal", n_nonDE_unimodal),
      rep("mixture", n_nonDE_mixture),
      rep("mixture", n_DE)
    ),
    de_status = c(
      rep("nonDE", n_nonDE_unimodal + n_nonDE_mixture),
      rep("DE", n_DE)
    ),
    gene_class = c(
      rep("nonDE_unimodal", n_nonDE_unimodal),
      rep("nonDE_mixture", n_nonDE_mixture),
      rep("DE_unimodal_vs_mixture", n_DE)
    ),
    mu_unimodal = c(mu_nonDE_unimodal, rep(NA_real_, n_nonDE_mixture), mu_unimodal_DE),
    mu_component_1 = c(rep(NA_real_, n_nonDE_unimodal), mu1_nonDE_mixture, mu1_DE),
    mu_component_2 = c(rep(NA_real_, n_nonDE_unimodal), mu2_nonDE_mixture, mu2_DE),
    theta = c(theta_nonDE_unimodal, rep(NA_real_, n_nonDE_mixture), theta_DE),
    theta_component_1 = c(rep(NA_real_, n_nonDE_unimodal), theta_nonDE_mixture, theta_DE),
    theta_component_2 = c(rep(NA_real_, n_nonDE_unimodal), theta_nonDE_mixture, theta_DE),
    mixing_prob = c(rep(NA_real_, n_nonDE_unimodal), mix_prob_nonDE_mixture, mix_prob_DE),
    fold_change = c(rep(NA_real_, n_nonDE_unimodal), fc_nonDE_mixture, fc_DE),
    de_unimodal_group = c(
      rep(NA_character_, n_nonDE_unimodal + n_nonDE_mixture),
      rep(de_unimodal_group, n_DE)
    ),
    shared_mode = c(
      rep(NA, n_nonDE_unimodal + n_nonDE_mixture),
      rep(share_mode, n_DE)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
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
      k1 = k1,
      k2 = k2,
      sigma_sq_nonDE = sigma_sq_nonDE,
      sigma_sq_DE = sigma_sq_DE,
      param_sampling = "truncated_lognormal"
    )
  )
}

