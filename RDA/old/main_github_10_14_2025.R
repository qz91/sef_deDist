library(dplyr) 
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(MASS) 
library(patchwork)
library(grid)
library(parallel)
library(truncnorm)
library(Seurat)
library(pbmcapply)
library(parallel)
library(doSNOW)
library(foreach) 
library(doParallel)
library(Ake)
library(Hmisc)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) 


#  ---- BIN COUNT HELPER FUNCTION ----
get_sc_bin_counts = function(Y, bin_edges, K){
  # function: tabulates expression values for each individual into bins with estimated grid points (length K)
  # input: list of arrays; grid points to use for bins
  # output: n x K matrix of bins
  n = length(Y)
  counts_matrix = matrix(0, nrow = n, ncol = K) # init n x K matrix to store bin counts
  
  for (i in seq_along(Y)) {
    # cut() classifies Y[[i]] elements into bins
    bin_indices = cut(Y[[i]], breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
    counts_matrix[i, ] = tabulate(bin_indices, nbins = K) # count occurrences in each bin
  }
  return(counts_matrix)
}

#  ---- PROCESSING HELPER FUNCTIONS  ----
subset_exprMat_by_donors = function(exprMat, metadata, donor_vector) {
  # function: subsets expression array based on cell IDs
  # input: vector/array of expression counts; cell-based metadata dataframe; list of IDs
  # output: subsetted array of expression counts that retain only the cells from the list of provided donors
  # assumption: assumes metadata and expression array have aligned indices
  selected_cells = rownames(metadata)[metadata$donor_id %in% donor_vector]
  exprMat_subset = exprMat[, selected_cells, drop = FALSE]
  return(exprMat_subset)
}

filter_genes_by_disease_donor_counts = function(donors_to_use_per_gene, covariate_df, min_donors = 50) {
  # function: filtration function to retain genes that have at least min_donors threshold
  # input: gene-specific named list of donors that satisfy filtration criteria; disease status-corresponding to sample IDs; minimum # of satisfying samples per group
  # function retains elements (list of donors for a gene i) where there are at least min_donors in each group
  if (!all(c("donor_id", "disease_status") %in% colnames(covariate_df))) {
    stop("covariate_df must contain 'donor_id' and 'disease_status' cols")
  }
  disease_levels = sort(unique(covariate_df$disease_status))  
  
  if (length(disease_levels) != 2) {
    stop("disease_status requires binary condition")
  }
  
  filtered_donors <- Filter(function(donors) {
    # quick function to extract counts of SLE and control counts for a gene
    matched_df = covariate_df[covariate_df$donor_id %in% donors, ] #keep donors based on the donor IDs in the ith element of donors_to_use_per_gene
    matched_df$disease_status = factor(matched_df$disease_status, levels = disease_levels) # convert to factor
    donor_counts = table(matched_df$disease_status) #make table of # SLE, # controls
    all(donor_counts >= min_donors)
  }, donors_to_use_per_gene)
  
  return(filtered_donors)
}

filtration_method = function(sObj, gene_vector, exprMat, n_cores = parallel::detectCores() - 1) {
  # function: preprocessing function that identifies, for each gene listed, donors with at least 100 cells and 20% non-zero expression
  # inputs: seurat object, list of genes, count expression matrix (gene by cell), cores to use
  # output: returns a list of arrays, each element of the encompassing list corresponds to a gene and the donors that satisfy filtration criteria (at least 100 cells, no less than 20% of counts are non-zero)
  donors = unique(sObj$donor_id)
  
  stopifnot(ncol(exprMat) == ncol(sObj)) # sanity check for alignment
  
  donor_ids = as.character(sObj$donor_id)
  
  # run in parallel over genes; must loop sequentially through all donors thru each gene
  results = pbmclapply(gene_vector, function(gene) {
    if (!(gene %in% rownames(exprMat))) return(NULL) #sanity check to make sure gene name is actually in the expression matrix
    
    gene_expr = exprMat[gene, ] #subset the expression matrix into 
    passing_donors = c() #vector of donors that satisfy the criteria
    
    for (donor in donors) {
      cell_idx = which(donor_ids == donor)
      if (length(cell_idx) < 100) next #must have at least 100 cells
      
      donor_expr = gene_expr[cell_idx]
      nonzero_prop = sum(donor_expr != 0) / length(donor_expr)
      
      if (nonzero_prop >= 0.2) { #must have at least 20% non-zero in expression
        passing_donors = c(passing_donors, donor)
      }
    }
    return(passing_donors)
  }, mc.cores = n_cores)
  
  names(results) = gene_vector #should be a list with names for each gene
  return(results) 
}

#  ---- RUN MODEL FUNCTIONS  ----
sef_regression = function(exprMat, sObj_meta, donor_list, p = 2, plot_flag = F, h = NULL, K = NULL){
  # function: conducts regression for a given gene; smoothing done with either a gaussian kernel or a gamma kernel
  # inputs: gene-specific expression array; cell-wise metadata; list of donors satisfying filtration criteria for the same gene; moments to test for; plot flag; bandwidth
  # output: returns "sef_regression_object" (list) that contains all objects used and computed for regression; retained for inference
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(exprMat = exprMat,metadata = sObj_meta, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation

  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ #bandwidth selection based on unique points
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
      print(K)
    } else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1) # grid to discretize expression array into counts
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2 # use the midpoint of the bin
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u) # gamma kernel
    dev.off()
    carrier_est = carrier$est.fn # estimated carrier density
    hist_plot_breaks = carrier$hist$breaks # used for plotting
  } else{
    K = 100
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    print(paste0("greater than or equal to 100 unique points: using equidistant bins as points for Gaussian KDE; carrier density estimate (", unique_pts," unique points), setting K bin count to: ", K))
    print(K)
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  if (plot_flag == T) {
    hist(y_agg,
         probability = TRUE,       
         breaks = hist_plot_breaks,
         col = "lightgray",
         border = "white",
         main = "Histogram + Carrier Density",
         xlab = "Transformed Expression",
         ylab = "Density")
    lines(est_midpoints, carrier_est, col = "dodgerblue", lwd = 2)
  }
  
  print("computing group 1 bin matrices")
  #counts in each bin for both groups
  control_donors =unique(sObj_meta[sObj_meta$disease == "normal", "donor_id"]) #group 1
  y1 = lapply(control_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  print("computing group 2 bin matrices")
  disease_donors = unique(sObj_meta[sObj_meta$disease != "normal", "donor_id"]) #group 2
  y2 = lapply(disease_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K) # bin matrix K columns for group 1
  Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) # bin matrix K columns for group 2
  
  binwidth = diff(est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  print("running SEF regression model")
  
  tk = est_midpoints # make design matrix for SEF regression
  X = rep(1, length(est_midpoints)) 
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  
  S1sum = colSums(Smat1) # y for group 1
  S2sum = colSums(Smat2) # y for group 2
  
  #setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) # scale factor group 1
  scale_factor2 = cellSum2/sum(carrier_est) # scale factor group 2
  carrier_scale1 = carrier_est*scale_factor1 # group 1 scale factor onto carrier density
  carrier_scale2 = carrier_est*scale_factor2 # group 2 scale factor onto carrier density
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  df1 = df1[df1$carrier_scale > 0, , drop = FALSE] # avoid issue of log(0)
  df2 = df2[df2$carrier_scale > 0, , drop = FALSE] # avoid issue of log(0)
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  
  # if experience convergence issue, we simply test on p-1 moments
  convergence_issue = (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # test p-1 moments
    varb = paste0("X_", 1:p)
    X = rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X = cbind(X, tk^dim)
    }
    # new formula
    formula = as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + ")))
    df1_reduced = df1[, -ncol(df1), drop = FALSE] # drop pth moment columns
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
  beta_est1 = as.vector(sef_gp1$coefficients) # get group 1 coefficients
  beta_est2 = as.vector(sef_gp2$coefficients) # get group 2 coefficients
  beta_diff = beta_est1 - beta_est2 # differences in coefficients
  
  
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) ) # get estimated group-wise density 
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) ) # get estimated group-wise density 
  
  df_plot = data.frame( # used for plotting
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  
  if (plot_flag == T) {
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  return(list(y1 = y1, y2 = y2, sef_df1 = sef_df1, sef_df2 = sef_df2, df1 = df1, df2 = df2, 
              carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid, hist_plot_breaks = hist_plot_breaks,
              X =X, K = K, p = p, n_total_cells = n_total_cells, binwidth = binwidth,
              Smat1 = Smat1, Smat2 = Smat2, Smat = Smat, cellSum1 = cellSum1, cellSum2 = cellSum2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2))
}
sef_inference = function(X = NULL, K = NULL, sef_df1 = NULL, sef_df2 = NULL, cellSum1 = NULL, cellSum2 = NULL, est_midpoints = NULL, binwidth = NULL, beta_est1 = NULL, beta_est2 = NULL, 
                          Smat1 = NULL, Smat2 = NULL, p = NULL, sef_regression_obj = NULL)
{
  # function: conducts inference after using SEF regression
  # input: sef_regression object that contains design matrix, bin count, fitted values, regression coefficients for the p moments tested, etc.
  # output: test statistic, p-value, covariance, matrix, difference in beta estimates
  if (!is.null(sef_regression_obj)) {   # extract components from sef_regression object
    X = sef_regression_obj$X
    K = sef_regression_obj$K
    n_total_cells = sef_regression_obj$n_total_cells
    sef_df1 = sef_regression_obj$sef_df1
    sef_df2 = sef_regression_obj$sef_df2
    cellSum1 = sef_regression_obj$cellSum1
    cellSum2 = sef_regression_obj$cellSum2
    est_midpoints = sef_regression_obj$est_midpoints
    binwidth = sef_regression_obj$binwidth
    beta_est1 = sef_regression_obj$beta_est1
    beta_est2 = sef_regression_obj$beta_est2
    Smat1 = sef_regression_obj$Smat1
    Smat2 = sef_regression_obj$Smat2
    p = sef_regression_obj$p
  }
  
  if (any(sapply(list(X, K, sef_df1, sef_df2, cellSum1, cellSum2, n_total_cells,
                      est_midpoints, binwidth, beta_est1, beta_est2, 
                      Smat1, Smat2, p), is.null))) {
    stop("Some required parameters are missing. Either provide them individually or via sef_regression_obj.")
  }
  beta_diff = beta_est2 - beta_est1
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  
  # marginal distribution test statistic computation
  G_1 = t(X) %*% (sef_df1 * X) * cellSum1 / sum(sef_df1) # G from Thm 4.1
  G_2 = t(X) %*% (sef_df2 * X) * cellSum2 / sum(sef_df2)
  
  pairwise_diff = outer(est_midpoints, est_midpoints, "-")
  n_total_cells = cellSum1 + cellSum2
  M = dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  
  Z11 = t(X) %*% (diag(K) - cellSum1 * as.vector(exp(X %*% beta_est1)) * M / sum(sef_df1)) # Z from Thm 4.1
  Z12 = t(X) %*% (-cellSum1 * as.vector(exp(X %*% beta_est1)) * M / sum(sef_df1))
  Z21 = t(X) %*% (-cellSum2 * as.vector(exp(X %*% beta_est2)) * M / sum(sef_df2))
  Z22 = t(X) %*% (diag(K) - cellSum2 * as.vector(exp(X %*% beta_est2)) * M / sum(sef_df2))
  
  L1 = solve(G_1, Z11) - solve(G_2, Z21)
  L2 = solve(G_1, Z12) - solve(G_2, Z22)
  
  m_vec1 = rowSums(Smat1)
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1 / m_vec1
  Smat2_ave = Smat2 / m_vec2
  
  Smat1_centered = t(Smat1_ave) - colSums(Smat1) / sum(m_vec1)
  Smat2_centered = t(Smat2_ave) - colSums(Smat2) / sum(m_vec2)
  
  tmp1 = diag(S1sum) - 2 * t(Smat1) %*% Smat1_ave + t(Smat1_ave) %*% Smat1_ave + 
    Smat1_centered %*% (t(Smat1_centered) * m_vec1^2)
  tmp2 = diag(S2sum) - 2 * t(Smat2) %*% Smat2_ave + t(Smat2_ave) %*% Smat2_ave + 
    Smat2_centered %*% (t(Smat2_centered) * m_vec2^2)
  
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  
  chi_stat1 = as.numeric(beta_diff[-1] %*% solve(Cov_bar[-1, -1], beta_diff[-1]))
  pval_1 = 1 - pchisq(chi_stat1, df = p)
  print(paste0("marginal distribution pval: ", pval_1))
  
  return(list(beta_diff = beta_diff, pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_bar = Cov_bar))
}

stablerunSeuratCounts_Concise = function(exprMat, sObj_meta, donor_list, p = 2, plot_flag = T, h = NULL){
  # function: put the regression and inference functions together; wrapper for entire inferential and density estimation process
  # input: expression array; cell-wise metadata, named list of ID arrays; moments to test; plot flag; bandwidth
  # output: regression statistics, etc.
  regression_object = sef_regression(exprMat = exprMat, sObj_meta = sObj_meta, donor_list = donor_list, p = p, plot_flag = plot_flag, h = h) # , transformation = transformation)
  test_object = sef_inference(sef_regression_obj = regression_object)
  to_remove = c("X", "Smat1", "Smat2", "Smat", "sef_gp1", "sef_gp2") # large and unnecessary at the end, so we omit 
  regression_statistics = regression_object[!(names(regression_object) %in% to_remove)]
  return(c(regression_statistics, test_object))
}

#  ---- PERMUTATION FUNCTIONS (changed to match sef_regression)  ----
sef_permuted_regression = function(exprMat, sObj_meta, donor_list, covariate_df, p = 2, plot_flag = F, h = NULL, seed = 100){
  # function: SEF regression formatted for permutation test; note that inference is the same, uses the original sef_inference() function
  # input: same as original sef_regression() with addition of a seed parameter for reproducibility 
  # output: permuted regression still returns an sef_regression lsit object
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(exprMat = exprMat, metadata = sObj_meta, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  #apply poisson VST based on sqrt + 1/2 transformation
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation
  #Gaussian KDE evaluated at integers
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ 
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
      print(K)
    } else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    dev.off()
    carrier_est = carrier$est.fn
    hist_plot_breaks = carrier$hist$breaks
  } else{
    K = 100
    print(paste0("greater than or equal to 100 unique points: using equidistant bins as points for Gaussian KDE; carrier density estimate (", unique_pts," unique points), setting K bin count to: ", K))
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  if (plot_flag == T) {
    hist(y_agg,
         probability = TRUE,     
         breaks = hist_plot_breaks,
         col = "lightgray",
         border = "white",
         main = "Histogram + Carrier Density",
         xlab = "Transformed Expression",
         ylab = "Density")
    lines(est_midpoints, carrier_est, col = "dodgerblue", lwd = 2)
  }
  
  print(paste0("permuting samples with SLE and healthy controls with seed ", seed, "..."))
  set.seed(seed) # for reproducibility when permuting disease status across samples
  covariate_df = subset(covariate_df, donor_id %in% donor_list)
  permuted = covariate_df
  permuted$disease_status = sample(permuted$disease_status)

  
  print("done")
  
  print("computing group 1 bin matrices")
  # counts in each bin for both groups
  control_donors = subset(permuted, disease_status == 0)$donor_id #group 1
  y1 = lapply(control_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  print("computing group 2 bin matrices")
  disease_donors = subset(permuted, disease_status == 1)$donor_id #group 2
  y2 = lapply(disease_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K)
  Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) 
  
  binwidth = diff(est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  print("running SEF regression model")
  tk = est_midpoints
  X = rep(1, length(est_midpoints)) 
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  S1sum = colSums(Smat1) # y for group 1
  S2sum = colSums(Smat2) # y for group 2
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
  df1 = df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 = df2[df2$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  
  convergence_issue = (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for p-1 moments
    varb = paste0("X_", 1:p)
    X = rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X = cbind(X, tk^dim)
    }
    formula = as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + "))) #same formula as sef_regression to adjust for convergence issues
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
  
  df_plot = data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  
  if (plot_flag == T) {
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  return(list(control_donors = control_donors, disease_donors = disease_donors, y1 = y1, y2 = y2, 
              sef_df1 = sef_df1, sef_df2 = sef_df2, df1 = df1, df2 = df2, 
              carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid,
              X =X, K = K, p = p, n_total_cells = n_total_cells, binwidth = binwidth,
              Smat1 = Smat1, Smat2 = Smat2, Smat = Smat, cellSum1 = cellSum1, cellSum2 = cellSum2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2))
}

stablePermutationTest = function(exprMat, sObj_meta, donor_list, covariate_df, p = 2, plot_flag = T, h = NULL, seed = 100){
  regression_object = sef_permuted_regression(exprMat = exprMat, sObj_meta = sObj_meta, donor_list = donor_list, covariate_df = covariate_df, p = p, plot_flag = plot_flag, h = h, seed = seed) 
  test_object = sef_inference(sef_regression_obj = regression_object)
  to_remove = c("y1","y2","sef_gp1","sef_gp2", "df1","df2", "X", "Smat1", "Smat2", "Smat", "sef_gp1", "sef_gp2") # too large and unnecessary at the end
  regression_statistics = regression_object[!(names(regression_object) %in% to_remove)]
  
  return(c(regression_statistics, test_object))
}

####### FIGURE FUNCTIONS ####### 
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
plot_carrier_with_hist = function(exprMat, gene, sObj_meta, donors_to_use_list, h = NULL){
  # function: plots carrier density with histogram overlayed
  # input: expression count matrix; gene of interest; cell-wise metadata; list of donor IDs; bandwidth
  # output: plot of carrier density with histogram overlaying
  exprMat = exprMat[gene, , drop = F]
  donor_list = donors_to_use_list[[gene]]
  
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(metadata = sObj_meta, exprMat = exprMat, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  #apply poisson VST based on sqrt + 1/2 transformation
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation
  #Gaussian KDE evaluated at integers
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ #uneven bins
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
    }
    else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    dev.off()
    carrier_est = carrier$est.fn
    hist_plot_breaks = carrier$hist$breaks
  } else{
    K = 100
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    print(carrier$bw)
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  p = ggplot() +
    geom_histogram(aes(x = y_agg, y = ..density..), #histogram density mode
                   breaks = hist_plot_breaks,
                   fill = "skyblue1", alpha = 0.75,color = "white") +
    geom_line(aes(x = est_midpoints, y = carrier_est), #carrier density
              color = "dodgerblue", size = 1) +
    labs(
      title = paste0("Carrier Density for ", gene),
      x = "Expression",
      y = "Density"
    ) + theme_minimal()
  
  return(p)
}
group_model_comparison_plot = function(sef_object, gene) {
  # function: plot a group comparison
  # input: sef object; gene of interest
  # output: plot comparing two densities; controls is blue; SLE is red

  sef_1 = sef_object[[gene]]$sef_df1
  sef_2 = sef_object[[gene]]$sef_df2
  est_grid = sef_object[[gene]]$est_grid
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
get_individual_sef_counts_model = function(donor_id, gene, exprMat, sObj_meta, donor_list, p = 2, h = NULL, plot_flag = F){
  # function: plot overlaying histogram of data for specific donor, overlaying individual-specific model and group-wise model to which the donor belongs to
  # input: donor ID; gene of interest; expression matrix; cell-wise metadata; list of donor IDs; moments to test; bandwidth; plot flag
  # output: plot
  
  # run global model first, receive the carrier density as well
  res = stablerunSeuratCounts_Concise(exprMat = exprMat[gene, , drop = F], sObj_meta = sObj_meta, donor_list = donor_list[[gene]], p = p, h = h, plot_flag = plot_flag) #group-specific model execution
  print(paste0("now running individual-specific model for donor ", donor_id))
  
  #get midpoints, carrier density estimate, gridpoints, status
  ind_carrier_est = res$carrier_est
  ind_est_midpoints = res$est_midpoints
  ind_est_grid = res$est_grid
  ind_status = as.character(sObj_meta[sObj_meta$donor_id == donor_id,]$disease[1])
  y_agg = c(unlist(res$y1), unlist(res$y2))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  n_total_cells = res$n_total_cells
  K = res$K
  h = res$h
  p = res$p
  hist_plot_breaks = res$hist_plot_breaks
  
  if (ind_status == "normal") {
    grp_model = res$sef_df1 # Normal
  }else{
    grp_model = res$sef_df2 # SLE
  }
  
  ind_sObj_meta = sObj_meta[sObj_meta$donor_id == donor_id, , drop = FALSE]
  ind_cell_barcodes = rownames(ind_sObj_meta) # to subset individual-specific vector of expression counts
  y_ind = as.vector(exprMat[gene, ind_cell_barcodes, drop = F])
  y_ind = sqrt(y_ind + 1/2)
  n_indiv_cells = length(y_ind)
  
  tk_ind = ind_est_midpoints
  X_ind = rep(1, length(ind_est_midpoints))
  for (dim in 1:p){
    X_ind = cbind(X_ind, tk_ind^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X_ind) = c("Intercept", varb)
  
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  
  Smat_ind = get_sc_bin_counts(Y = list(y_ind), bin_edges = ind_est_grid, K = K)
  ind_SSum = colSums(Smat_ind)
  binwidth = diff(ind_est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  
  # setup for modeling individual-wise
  cellSum_ind = length(unlist(y_ind))
  scale_factor_ind = cellSum_ind/sum(ind_carrier_est)
  ind_carrier_scale = ind_carrier_est*scale_factor_ind
  df1 = as.data.frame(cbind( ind_SSum, ind_carrier_scale, X_ind))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  df1 = df1[df1$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+')))
  sef_gp_ind = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  
  beta_est_ind = as.vector(sef_gp_ind$coefficients)
  sef_df1 = as.vector( ind_carrier_est * exp(X_ind %*% beta_est_ind) )
  
  # plotting machinery; make dataframe for ggplot2
  p_df = data.frame(
    x = ind_est_midpoints, 
    y = sef_df1, # donor-specific estimates
    group = "Donor-Specific"
  )
  p_df2 = data.frame(
    x = ind_est_midpoints,
    y = grp_model, # group-specific estimates
    group = "Group Model"
  )
  p_df_combined = rbind(p_df, p_df2) #dfs used for ggplot()
  df_hist = data.frame(y_ind = y_ind)
  
  plt = ggplot() +
    geom_histogram(data = df_hist, aes(x = y_ind, y = ..density..), # histogram density mode
                   breaks = hist_plot_breaks,
                   fill = "skyblue1", alpha = 0.75, color = "white") +
    geom_line(data = p_df_combined, aes(x = x, y = y, color = group), size = 1.2) + #group versus donor specific density estimates
    scale_color_manual(values = c("Donor-Specific" = "dodgerblue",
                                  "Group Model" = "firebrick")) +
    theme_minimal() +
    labs(
      x = "Expression",
      y = "Density",
      title = paste0("Donor ", donor_id, " - Specific Density: ", gene),
      color = "Model"
    )
  return(plt)
  print("done")
  
}
run_go_enrichment_plot = function(gene_symbols, orgDb = org.Hs.eg.db, ontology = "BP", 
                                   cutoff = 0.6, universe = NULL, pvalueCutoff = 0.05,
                                   minGSSize = 10, maxGSSize = 500) {
  # function: run enrichment analysis with clusterProfiler; converts from HGNC to ENTREZID and back 
  # clusterProfiler hypergeometric test for over-representation analysis
  # input: list of gene symbols; organism used (human); ontology aspect; simplify() cutoff, etc.
  # output: clusterProfiler enrichment object
  require(clusterProfiler)
  require(enrichplot)
  require(org.Hs.eg.db) 
  gene_df = clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDb)
  if( is.null(universe) == F){
    entrez_universe = clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDb)$ENTREZID
  } else{
    entrez_universe = NULL
  }
  
  if (nrow(gene_df) == 0) {
    stop("Cannot reformat HGNC symbols to ENSG IDs. Re-check gene list.")
  }
  
  # GO enrichment hypergeometric test
  ego = enrichGO(gene = gene_df$ENTREZID, OrgDb = orgDb,
                  universe = entrez_universe, pvalueCutoff = pvalueCutoff,
                  minGSSize = minGSSize, maxGSSize = maxGSSize,
                  keyType = "ENTREZID", ont = ontology,
                  pAdjustMethod = "BH", readable = TRUE)
  # run simplify() to remove redundant enriched terms to identify general BPs occurring
  ego = clusterProfiler::simplify(ego, cutoff = cutoff, by = "p.adjust")
  return(ego)
}



