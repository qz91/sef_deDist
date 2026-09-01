library(dplyr) 
library(ggplot2)
library(tidyverse)
library(mvtnorm)
library(MASS) 
library(patchwork)
library(grid)
library(parallel)
library(truncnorm)
library(pbmcapply)
library(parallel)
library(doSNOW)
library(foreach)
library(doParallel) 
library(Hmisc)
library(Matrix)

########## 1. DATA GENERATION FUNCTIONS ##########
simu_poigamma_vary_cellcts = function(seed = 100, lower = 300, upper = 1000, alpha1, alpha2, beta1, beta2, n1, n2){
  # function: simulate Poisson-Gamma mixture model data with alpha and beta are the shape and size parameters
  # input: seed; lower bound; upper bound; alpha group 1/2; beta group 1/2, number of samples in group 1/2
  # output: list of 2 lists (length n1 and n2) with each element being an array 
  Y1_list = vector("list", n1) 
  Y2_list = vector("list", n2)
  set.seed(seed = seed)
  for (i in 1:n1) { #GROUP 1
    n_cells_i = sample(lower:upper, 1)
    lambda1i = rgamma(1, shape = alpha1, rate = beta1)
    Y1_list[[i]] = rpois(n_cells_i, lambda1i)
  }
  for (i in 1:n2) { #GROUP 2 
    n_cells_i = sample(lower:upper, 1)
    lambda2i = rgamma(1, shape = alpha2, rate = beta2)
    Y2_list[[i]] <- rpois(n_cells_i, lambda2i)
  }  
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}
simu_zinb_vary_cellcts = function(n1, n2, mu, mu2, theta, sigma_sq = 0.5, pi = 0.5, pi2 = 0.5,
                                   b_pi = 1, b_theta = 1, de_type = 'dispersion', seed = 100, lower = 300, upper = 1000) {
  # function: simulate Zero-inflated Negative Binomial model data, offering dispersion or mean
  # input: seed; lower bound; upper bound; alpha group 1/2; beta group 1/2, number of samples in group 1/2
  # output: list of 2 lists (length n1 and n2) with each element being an array 
  set.seed(seed)
  # Output lists
  Y1_list = vector("list", n1)
  Y2_list = vector("list", n2)
  
  # Individual-specific means
  mu_i_group1 = rtruncnorm(n1, a = 0, mean = mu, sd = sqrt(sigma_sq))
  
  if (de_type == 'dispersion') {
    mu2 = mu
    theta2 = theta * b_theta
    pi2 = pi
    mu_i_group2 = rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
    
  } else if (de_type == "mean") {
    theta2 = theta
    pi2 = pi
    mu_i_group2 = rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
  } else {
    stop('Invalid de_type specified. Choose either "mean", "dispersion".')
  }
  
  # Default ZINB simulation for mean/dispersion de
  for (i in seq_len(n1)) { # group 1
    n_cells_i = sample(lower:upper, 1)
    Y1_list[[i]] = sapply(seq_len(n_cells_i), function(j) {
      if (runif(1) < pi) 0 else rnbinom(1, size = theta, mu = mu_i_group1[i])
    })
  }
  
  for (i in seq_len(n2)) { # group 2
    n_cells_i = sample(lower:upper, 1)
    Y2_list[[i]] = sapply(seq_len(n_cells_i), function(j) {
      if (runif(1) < pi2) 0 else rnbinom(1, size = theta2, mu = mu_i_group2[i])
    })
  }
  
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}

########## 2. SEF HELPER FUNCTIONS ##########
#  ---- DISCRETE GAUSSIAN KERNEL FUNCTION ----
continuous_discrete_kde = function(y, eval_points, bw = 0.25) {
  # function: continuous gaussian kernel density estimate evaluated at the integer values the data is in 
  # input: expression array, points to evaluate (similar to est_midpoints)
  # output: density estimate evaluated at parameter-passed points
  dx <- diff(eval_points)[1]  # step size between grid points; assumes equally spaced grid
  raw_density <- sapply(eval_points, function(x) {
    mean(dnorm((x - y) / bw)) / bw
  })
  normalized_density <- raw_density / sum(raw_density * dx)
  return(normalized_density)
}
#  ---- BIN COUNT HELPER FOR DGK ----
get_countdata_bin_counts = function(Y, grid_points){
  # function: tabulates expression values for each individual into bins with estimated grid points
  # input: list of arrays; grid points to use for bins
  # output: n x K matrix of bins, K is # of grid points
  
  n = length(Y)
  K = length(grid_points)
  counts_matrix = matrix(0, nrow = n, ncol = K) # init n x K matrix to store bin counts
  
  for (i in seq_along(Y)) {
    tab = table(factor(Y[[i]], levels = grid_points)) #factor() to ensure that values potentially not in a bin are tabulated still before counting
    counts_matrix[i, ] = as.numeric(tab) #revert back to integers
  }
  
  return(counts_matrix)
}

#  ---- POWER ANALYSIS SPECIFIC  ----

gamma_params = function(mu, V) {
  # function: given mean and variance, provides alpha and beta values to match the params
  # input: mean; variance
  # output: alpha; beta
  if(V <= mu) stop("Variance must be greater than the mean for overdispersion.")
  alpha = mu^2/(V - mu)
  beta  = mu/(V - mu)
  return(list(alpha = alpha, beta = beta))
}

compute_moments = function(y, p) {
  return(sapply(1:p, function(k) y^k))
} 

mom_cov_diff_cell_cts = function(Y, p) {
  # function: compute covariance matrix for Moments (for MoM estimator) with different cell counts
  # input: list of arrays of expression values; moments to test
  # output: covariance matrix (p x p)
  n = length(Y)  # number of individuals (length of the list)
  total_samples = 0
  overall_sum = rep(0, p)  # sum of overall mean
  within_var_sum = matrix(0, nrow = p, ncol = p)  # within-individual variance contributions
  between_var_sum = matrix(0, nrow = p, ncol = p)  # between-individual variance contributions 
  mean_list = matrix(0, nrow = n, ncol = p)  # matrix of p moments for each individual
  m_vector = numeric(n)  # individual-wise cell counts
  
  # within-individual variance
  for (i in 1:n) {  # iterate through individual
    samples_i = Y[[i]] # expression array for individual i
    m_i = length(samples_i) 
    m_vector[i] = m_i
    mean_samples_i = mean(samples_i)
    
    moments_i = t(sapply(samples_i, compute_moments, p = p))  # ncell_i x p moment 
    
    mean_t_i = colMeans(moments_i)  # compute mean of p moments for individual
    mean_list[i, ] = mean_t_i
    
    # update overall sum for computing overall mean
    overall_sum = overall_sum + m_i * mean_t_i
    total_samples = total_samples + m_i
    
    # compute within-individual variance
    diffs = t(moments_i) - mean_t_i  # (T_i - Tbar_i)
    within_var_sum = within_var_sum + diffs %*% t(diffs)  # sum (T_i - Tbar_i)^2
  }
  
  overall_mean <- overall_sum / total_samples # compute overall mean of moments
  
  # between-individual variance
  for (i in 1:n) {
    m_i = m_vector[i]
    mean_t_i = mean_list[i, ]
    diff <- mean_t_i - overall_mean
    between_var_sum <- between_var_sum + m_i^2 * (diff %*% t(diff))  # m_i^2, not m_i
  }
  
  Cov = within_var_sum / total_samples^2 + between_var_sum / total_samples^2 # covariance calculation
  return(Cov)
}

compute_empirical_power = function(pval_matrix, alpha = 0.05){
  # function: compute empirical power using the #of rejections we make divided by the total number of simulations
  # input: pvalue matrix
  # output: array of type i error values
  pmat = as.matrix(pval_matrix)
  return(colMeans(pmat < alpha))
}

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
########## 3. SEF FUNCTIONS ##########
#  ---- POISSON GAMMA  ----
simulatePGRealistic = function(repID = NULL, de_type = "mean", p = 2,  K = NULL,  lower = 300, upper = 1000, n1 = 100, n2 = 100, plot_flag = F, idx, alpha1, alpha2, beta1, beta2){
  # function: run Poisson-Gamma model
  # input: parameters of interest
  # output: 
  if (de_type == "mean") {
    effect_size = alpha2/beta2 - alpha1/beta1
  } else if(de_type == "variance"){
    var1 = alpha1/beta1 + alpha1/(beta1^2)
    var2 = alpha2/beta2 + alpha2/(beta2^2)
    effect_size = var2 - var1 #variance test
  } else{
    warning("Shift type not found. Please use valid shift type.")
    return(NULL)
  }
  seed = idx + 62
  set.seed(seed)
  simData = simu_poigamma_vary_cellcts(seed =seed, lower = lower, upper = upper, alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2,n1 = n1, n2 = n2)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  print(paste0("total cells: ", n_total_cells))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  
  if (is.null(K)) {
    est_grid = seq(l,u)
    est_midpoints = est_grid
    K = length(est_midpoints)
    print(paste0("Number of midpoints K not specified. Using total number of boundary points in pooled data as new K: ", K))
    
    carrier_est = continuous_discrete_kde(y = y_agg, eval_points = est_midpoints, bw = 0.05)
    
  } else{
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
  }
  
  hist_df <- data.frame(y_agg = y_agg)
  breaks <- seq(l, u, length.out = K + 2)
  line_df <- data.frame(x = est_midpoints, y = carrier_est)
  carrier_plot = NULL
  
  if(plot_flag == T){
    carrier_plot <- ggplot(hist_df, aes(x = y_agg)) +
      geom_histogram(aes(y = ..density..),
                     breaks = breaks,
                     fill = "skyblue1", alpha = 0.75,
                     color = "white") +
      geom_line(data = line_df, aes(x = x, y = y),
                color = "dodgerblue", linewidth = 1) +
      labs(
        title = "",
        x = "Expression",
        y = "Density"
      ) +
      theme_minimal()
    print(carrier_plot)
  }
  
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  Smat1 = get_countdata_bin_counts(Y = Y1, grid_points = est_grid)
  Smat2 = get_countdata_bin_counts(Y = Y2, grid_points = est_grid)
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  
  print("running SEF regression model")
  #Construct design matrix
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)

  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  #setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) 
  scale_factor2 = cellSum2/sum(carrier_est)
  carrier_scale1 = carrier_est*scale_factor1 
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  df1 <- df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 <- df2[df2$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  # fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  
  # handle issue for overfitting
  convergence_issue <- (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for one fewer 
    print(paste0("new test with p = ", p))
    varb = paste0("X_", 1:p)
    X = rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X = cbind(X, tk^dim)
    }
    formula = as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + "))) # rebuild formula
    df1_reduced = df1[, -ncol(df1), drop = FALSE] # drop last columns
    df2_reduced = df2[, -ncol(df2), drop = FALSE]
    
    sef_gp1 = tryCatch(glm(formula, family = poisson(link = "log"), data = df1_reduced),
                        error = function(e) {
                          message("Retry for sef_gp1 failed: ", conditionMessage(e))
                          return(NULL)
                        })
    sef_gp2 = tryCatch(glm(formula, family = poisson(link = "log"), data = df2_reduced),
                        error = function(e) {
                          message("Retry for sef_gp2 failed: ", conditionMessage(e))
                          return(NULL)
                        })
  }
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  
  # plot differences
  df_plot = data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  plt = NULL
  if(plot_flag == T){
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  # marginal distribution 
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1) 
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M = dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  Z11 = t(X) %*% ( diag(K) -  cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # Z_{t1, j} for j in T1
  Z12 = t(X) %*% ( -cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # For j in T2
  Z21 = t(X) %*% ( -cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2)  ) # For j in T1
  Z22 = t(X) %*% ( diag(K) -  cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2) ) # For j in T2
  L1 = solve(G_1, Z11) - solve(G_2, Z21) # For S1 data 
  L2 = solve(G_1, Z12) - solve(G_2, Z22)# For S2 data  
  m_vec1 = rowSums(Smat1) 
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1/m_vec1
  Smat2_ave = Smat2/m_vec2
  Smat1_centered = t(Smat1_ave) - colSums(Smat1)/sum(m_vec1)  
  Smat2_centered = t(Smat2_ave) - colSums(Smat2)/sum(m_vec2)
  tmp1 = diag(S1sum) - 2* t(Smat1)%*%Smat1_ave + t(Smat1_ave)%*%Smat1_ave + Smat1_centered%*%( t(Smat1_centered) * m_vec1^2 )
  tmp2 = diag(S2sum) - 2* t(Smat2)%*%Smat2_ave + t(Smat2_ave)%*%Smat2_ave + Smat2_centered%*%( t(Smat2_centered) * m_vec2^2 )
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  chi_stat1 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1], beta_diff[-1]) ) #marginal distribution test
  pval_1 = 1 - pchisq( chi_stat1, df = p)
  print(pval_1)
  
  return(list(de_type = de_type, effect_size = effect_size, y_agg = y_agg, est_midpoints = est_midpoints, carrier_est = carrier_est, binwidth = binwidth, K = K, y1 = y1, y2 = y2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, beta_diff = beta_diff,
              sef_df1 = sef_df1, sef_df2 = sef_df2, pval_1 = pval_1, Cov_1 = Cov_bar,
              carrier_plt = carrier_plot, comparison_plt = plt, breaks = breaks))
}
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL  ----
simulateZINB = function(n1, n2, de_type = "dispersion", idx, mu, mu2, sigma_sq = 0.5, theta, b_theta = 1, b_pi = 1, pi = 0.5, pi2 = 0.5, lower = 300, upper = 1000, K = NULL, p = 2, repID = NULL, plot_flag = F){
  #save values: effect size, binwidth, betas, covariances
  if (de_type == "dispersion") {
    theta2 = b_theta*theta
    effect_size = theta2 - theta
  } else if(de_type == "mean"){
    effect_size = mu2 - mu
  } else{
    stop(paste0("Please use valid de_type."))
  }
  seed = idx + 62
  simData = simu_zinb_vary_cellcts(n1 = n1, n2 = n2, lower = lower, upper = upper, mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, pi = pi, pi2 = pi2, 
                                   b_pi = b_pi, b_theta = b_theta, de_type = de_type, seed = seed)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  
  n_total_cells = length(y_agg)
  print(paste0("total cells: ", n_total_cells))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  
  if (is.null(K)) {
    est_grid = seq(l,u)
    est_midpoints = est_grid
    K = length(est_midpoints)
    print(paste0("Number of midpoints K not specified. Using total number of boundary points in pooled data as new K: ", K))
    carrier_est = continuous_discrete_kde(y = y_agg, eval_points = est_midpoints, bw = 0.05)
  } else{
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
  }
  
  hist_df = data.frame(y_agg = y_agg)
  breaks = seq(l, u, length.out = K + 2)
  line_df = data.frame(x = est_midpoints, y = carrier_est)
  
  carrier_plot = NULL
  if(plot_flag == T){
    carrier_plot = ggplot(hist_df, aes(x = y_agg)) +
      geom_histogram(aes(y = ..density..),
                     breaks = breaks,
                     fill = "skyblue1", alpha = 0.75,
                     color = "white") +
      geom_line(data = line_df, aes(x = x, y = y),
                color = "dodgerblue", linewidth = 1) +
      labs(
        title = "",
        x = "Expression",
        y = "Density"
      ) +
      theme_minimal()
    print(carrier_plot)
  }
  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  Smat1 = get_countdata_bin_counts(Y = Y1, grid_points = est_grid)
  Smat2 = get_countdata_bin_counts(Y = Y2, grid_points = est_grid)
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  
  print("running SEF regression model")
  # construct design matrix
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  if (unique_pts < 75) {
    binwidth = median(diff(est_midpoints))
  }else{
    binwidth = diff(est_midpoints)[1] # is identical to (u - l)/K in large enough settings
  }
  
  # setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) 
  scale_factor2 = cellSum2/sum(carrier_est)
  carrier_scale1 = carrier_est*scale_factor1 
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  df1 = df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 = df2[df2$carrier_scale > 0, , drop = FALSE]
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  # fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  
  convergence_issue <- (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for one fewer 
    print(paste0("new test with p = ", p))
    varb = paste0("X_", 1:p)
    X = rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X = cbind(X, tk^dim)
    }
    formula = as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + "))) # rebuild formula
    df1_reduced = df1[, -ncol(df1), drop = FALSE] 
    df2_reduced = df2[, -ncol(df2), drop = FALSE]
    
    sef_gp1 = tryCatch(glm(formula, family = poisson(link = "log"), data = df1_reduced),
                        error = function(e) {
                          message("Retry for sef_gp1 failed: ", conditionMessage(e))
                          return(NULL)
                        })
    sef_gp2 = tryCatch(glm(formula, family = poisson(link = "log"), data = df2_reduced),
                        error = function(e) {
                          message("Retry for sef_gp2 failed: ", conditionMessage(e))
                          return(NULL)
                        })
  }
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  
  # plot differences
  df_plot = data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  plt = NULL
  if(plot_flag == T){
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  # marginal distribution
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1) 
  G_1
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2)
  G_2
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M = dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  Z11 = t(X) %*% ( diag(K) -  cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # Z_{t1, j} for j in T1
  Z12 = t(X) %*% ( -cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # For j in T2
  Z21 = t(X) %*% ( -cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2)  ) # For j in T1
  Z22 = t(X) %*% ( diag(K) -  cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2) ) # For j in T2
  L1 = solve(G_1, Z11) - solve(G_2, Z21) # For S1 data
  L2 = solve(G_1, Z12) - solve(G_2, Z22)# For S2 data
  m_vec1 = rowSums(Smat1)
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1/m_vec1
  Smat2_ave = Smat2/m_vec2
  Smat1_centered = t(Smat1_ave) - colSums(Smat1)/sum(m_vec1)
  Smat2_centered = t(Smat2_ave) - colSums(Smat2)/sum(m_vec2)
  tmp1 = diag(S1sum) - 2* t(Smat1)%*%Smat1_ave + t(Smat1_ave)%*%Smat1_ave + Smat1_centered%*%( t(Smat1_centered) * m_vec1^2 )
  tmp2 = diag(S2sum) - 2* t(Smat2)%*%Smat2_ave + t(Smat2_ave)%*%Smat2_ave + Smat2_centered%*%( t(Smat2_centered) * m_vec2^2 )
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  
  chi_stat1 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1], beta_diff[-1]) ) #marginal distribution test
  pval_1 = 1 - pchisq( chi_stat1, df = p)
  print(pval_1)
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, n1 = n1, n2 = n2, binwidth = binwidth, X = X, y_agg = y_agg, y1 = y1, y2 = y2,
              carrier_est = carrier_est, est_midpoints = est_midpoints, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2,
              sef_df1 = sef_df1, sef_df2 = sef_df2, pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_1 = Cov_bar, beta_diff = beta_diff,
              carrier_plt = carrier_plot, comparison_plt = plt))
}
########## 4. COMPETING METHOD FUNCTIONS ##########
#  ---- METHOD OF MOMENTS  ----
simulateMoMRealistic = function(distribution = c("pg", "zinb"), de_type = "mean", p = 2, repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
  # function: run method of moments estimator
  # input: distribution type; shift type (mean, variance, dispersion); moments tested; replicate ID; sample sizes; index; lower and upper cell count bounds
  # output: test statistic and p-value; distribution; true effect size; group-wise covariance matrices; time elapsed in computation
  seed = idx + 62
  set.seed(seed)

  effect_size  = NULL
  extra_args <- list(...)
  
  if (distribution == "pg") {
    # simu_poigamma uses alpha1, alpha2, beta1, beta2
    simData = simu_poigamma_vary_cellcts(seed = seed, lower = 300, upper = 1000, n1 = n1, n2 = n2,...)
  } else if (distribution == "zinb"){
    # simu_zinb uses mu, mu2, theta, b_theta for dispersion
    effect_size = switch(de_type,
                          dispersion = {
                            theta2 <- b_theta * theta
                            theta2 - theta
                          },
                          mean = {
                            warning("Mean shift inherently includes variance shift")
                            mu2 - mu
                          },
                          {
                            warning("Invalid de_type for zinb; returning NULL")
                            NULL
                          })
    simData = simu_zinb_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper,  de_type = de_type, ...)
  } else{
    print("Invalid distribution type")
    return(NULL)
  }
  
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) 
  y2 = as.vector(unlist(Y2))
  
  start.time.mom = Sys.time()
  Cov1 = mom_cov_diff_cell_cts(Y1, p) # group 1 MoM covariance estimator
  Cov2 = mom_cov_diff_cell_cts(Y2, p) # group 2 MoM covariance estimatorv
  T1 = y1
  for (dim in 2:p){
    T1 = cbind(T1, y1^dim)
  }
  T2 = y2
  for (dim in 2:p){
    T2 = cbind(T2, y2^dim)
  }
  T1_bar = colMeans(T1)
  T2_bar = colMeans(T2)
  
  chisq_mom = as.numeric(t((T1_bar - T2_bar)) %*% diag((1 / diag(Cov1 + Cov2))) %*% (T1_bar - T2_bar)) #chi sq test statistic
  pv_mom = 1 - pchisq(chisq_mom, df = p) # we use p degrees of freedom
  end.time.mom = Sys.time()
  time.taken.mom = end.time.mom - start.time.mom
  return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, pval = pv_mom, chisq_mom = chisq_mom, Cov1 = Cov1, Cov2 = Cov2, p = p, time_taken = time.taken.mom))
}
#  ---- PSEUDOBULK  ----
simulatePseudobulkRealistic = function(distribution = c("pg", "zinb"), test_type = c("ttest","ks","ftest"), de_type = "mean", repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
  # function: run pseudobulk average estimator; we opt pseudobulk average due to nature of testing 1 gene at a time, so aggregation with library size is not suitable
  # input: distribution type; shift type (mean, variance, dispersion); replicate ID; sample sizes; index; lower and upper cell count bounds
  # output: test statistic and p-value; distribution; true effect size; pseudobulk data for each group
  seed = idx + 62
  set.seed(seed)

  effect_size  = NULL
  extra_args <- list(...)
  for (argname in names(extra_args)) {
    assign(argname, extra_args[[argname]])
  }
  
  if (distribution == "pg") {
    # simu_poigamma uses alpha1, alpha2, beta1, beta2
    simData = simu_poigamma_vary_cellcts(seed = seed, lower = lower, upper = upper, n1 = n1, n2 = n2,...)
  } else if (distribution == "zinb"){
    # simu_zinb uses mu, mu2, theta, b_theta for dispersion
    effect_size = switch(de_type,
                          dispersion = {
                            theta2 <- b_theta * theta
                            theta2 - theta
                          },
                          mean = {
                            warning("Mean shift inherently includes variance shift")
                            mu2 - mu
                          },
                          {
                            warning("Invalid de_type for zinb; returning NULL")
                            NULL
                          })
    simData = simu_zinb_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper,  de_type = de_type, ...)
  } else{
    print("Invalid distribution type")
    return(NULL)
  }
  
  Y1 = simData$Y1 
  Y2 = simData$Y2
  pbY1 = vapply(Y1, mean, numeric(1)) # pb values
  pbY2 = vapply(Y2, mean, numeric(1)) # pb values
  pb.df1 = data.frame(pbY1, rep(0,length(n1)))
  colnames(pb.df1) = c("pseudobulk", "group")
  pb.df2 = data.frame(pbY2, rep(1,length(n2)))
  colnames(pb.df2) = c("pseudobulk", "group")
  pb.total = rbind(pb.df1, pb.df2)
  
  if (test_type == "ttest") { # standard t-test
    pb.test = t.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    pv_pb = pb.test$p.value
    mean1_pb = pb.test$estimate[1]
    mean2_pb = pb.test$estimate[2]
    test_stat_pb = pb.test$statistic
    sd_pb = pb.test$stderr
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, mean1_pb = mean1_pb, mean2_pb = mean2_pb, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, sd = sd_pb, test_type = test_type))
    
  } else if(test_type == "ks"){ # standard komolgorov-smirnov test 
    pb.test = ks.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    pv_pb = pb.test$p.value
    test_stat_pb = pb.test$statistic
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
    
  } else if(test_type == "ftest"){ # f-test
    pb.test = var.test(pb.df1$pseudobulk, pb.df2$pseudobulk, ratio = 1, alternative = "two.sided")
    test_stat_pb = pb.test$statistic
    pv_pb = pb.test$p.value
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
    
  } else{
    warning("Please use valid test")
    return(NULL)
  }
}

