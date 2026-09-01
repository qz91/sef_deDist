# supplementary simulations
rm(list = ls())  
library(dplyr) 
library(truncnorm)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(MASS) 
library(patchwork)
library(grid)
library(parallel)
library(doSNOW)
library(foreach)
library(reshape)
source("./density_estimation/simulations/simulations_github/simulations_main_10_13_2025.r") 

########## Fig S.1 CARRIER DENSITY ##########
#  ---- POISSON GAMMA  ----
DataType = "pg"
n1 = 100; n2 = 100
p = 2
alpha1 = 20 # change to 5 for additional viz
alpha2 = 20 # change to 5 for additional viz
beta1 = 2; beta2 = 2
repID = 77  
test = simulatePGRealistic(
  repID = repID,
  idx = 1,
  alpha1 = alpha1, alpha2 = alpha2,
  beta1 = beta1, beta2 = beta2,
  n1 = n1, n2 = n2,
  de_type = "mean", p = p,
  plot_flag = T
)
GraphName = paste0("./carrier-density-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
print(test$carrier_plt)
dev.off() 
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL  ----
DataType = "zinb"
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_theta = 1
pi = 0.5 ; pi2 = 0.5
b_pi = 1
de_type = "dispersion"
repID = 77
p = 2
test = simulateZINB(
  repID = repID, idx = 1,
  mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
  theta = theta, b_theta = 1,
  pi = pi, pi2 = pi2,
  n1 = n1, n2 = n2,
  de_type = de_type, p = p,
  plot_flag = T
)
GraphName = paste0("./carrier-density-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
print(test$carrier_plt)
dev.off() 



########## Top Half Fig S.2, S.3 QQ Plots ##########
#  ---- POISSON GAMMA p = 3  ----
DataType = "pg"
maxIter = 300 
n1 = 100; n2 = 100
p = 3
alpha1 = 20; alpha2 = 20
beta1 = 2; beta2 = 2
repID = 77  
pvVec_margin = NULL 
pv_list = pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res = simulatePGRealistic(
      repID = repID,
      idx = i,
      alpha1 = alpha1, alpha2 = alpha2,
      beta1 = beta1, beta2 = beta2,
      n1 = n1, n2 = n2,
      de_type = "mean",
      p = p
    )
    pv = res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = 4, mc.set.seed = TRUE)

pvVec_margin = unlist(pv_list, use.names = FALSE)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./varying_cell_cts_qqplot-15-10-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Poisson Gamma (p = 3)")
print(realistic_qqplot)
dev.off() 

#  ---- POISSON GAMMA n1 = 100 n2 = 200  ----
DataType = "pg"
maxIter = 300 
n1 = 100; n2 = 200
p = 2
alpha1 = 20; alpha2 = 20
beta1 = 2; beta2 = 2
repID = 77  
pvVec_margin = NULL 
pv_list = pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res <- simulatePGRealistic(
      repID = repID,
      idx = i,
      alpha1 = alpha1,
      alpha2 = alpha2,
      beta1 = beta1,
      beta2 = beta2,
      n1 = n1, n2 = n2,
      de_type = "mean",
      p = p
    )
    pv = res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = 4, mc.set.seed = TRUE)
pvVec_margin = unlist(pv_list, use.names = FALSE)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./varying_cell_cts_qqplot-15-10", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Poisson Gamma Model")
print(realistic_qqplot)
dev.off() 

#  ---- ZINB p = 3  ----
DataType = "zinb"
maxIter = 300 
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_theta = 1
pi = 0.5 ; pi2 = 0.5
b_pi = 1
de_type = "dispersion"
repID = 77
p = 3
cores = 5
pvVec_margin = NULL 
pv_list = pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res = simulateZINB(
      repID = repID, idx = i,
      mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
      theta = theta, b_theta = 1,
      pi = pi, pi2 = pi2,
      n1 = n1, n2 = n2,
      de_type = de_type, p = p
    )
    pv = res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  
  })
}, mc.cores = cores, mc.set.seed = TRUE)
pvVec_margin = unlist(pv_list, use.names = FALSE)

gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB (p = 3)")
GraphName = paste0("./varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot =gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB (p = 3)")
print(realistic_qqplot)
dev.off() 


#  ---- ZINB n1=100, n2=200  ----
DataType = "zinb"
maxIter = 300 
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 200
sigma_sq = 5
b_theta = 1
pi = 0.5 ; pi2 = 0.5
b_pi = 1
de_type = "dispersion"
repID = 77
p = 2
cores = 5
pvVec_margin = NULL 
pv_list = pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res = simulateZINB(
      repID = repID, idx = i,
      mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
      theta = theta, b_theta = 1,
      pi = pi, pi2 = pi2,
      n1 = n1, n2 = n2,
      de_type = de_type, p = p
    )
    pv = res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL 
  })
}, mc.cores = cores, mc.set.seed = TRUE)
pvVec_margin = unlist(pv_list, use.names = FALSE)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot =gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB Model")
print(realistic_qqplot)
dev.off() 

########## Bottom Half Fig S.2, S.3 Power Analysis ##########
#  ---- POISSON GAMMA VARIANCE p = 3----
DataType = "pg"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 3
repID = 77 
mu1  = 10
V1   = 15
params1 = gamma_params(mu1, V1)
alpha1 = params1$alpha
beta1  = params1$beta
target_means = 10
target_variance = seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s = target_means^2 / (target_variance - target_means)
beta2s  = target_means / (target_variance - target_means)
result_df = data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 = alpha2s[idx]
  beta2 = beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }

  pvMat_margin = cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom  = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_variance - V1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))

df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL VARIANCE  p = 3 ----
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_thetas = seq(1, 1.70, by =  0.07)
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "dispersion"
repID = 77
K = 100
p = 3
DataType = "zinb"
maxIter = 300 
nCPUS = 8
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' )
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(b_thetas))
for (idx_b in 1:length(b_thetas)) {
  b_theta = b_thetas[idx_b]
  print(paste0("scale factor b_theta =", b_theta, " or (variance) = ", b_theta*theta))
  
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    tryCatch({
      simulateZINB(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                   b_pi = b_pi, pi = pi, pi2 = pi2,
                   repID, idx = i, n1 = n1, n2 = n2, de_type = de_type, K = NULL, p = p)
    }, error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(NULL)
    })
  }
  output1 = Filter(Negate(is.null), output)
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    simulateMoMRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                         b_pi = b_pi, pi = pi, pi2 = pi2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID, n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                           b_pi = b_pi, pi = pi, pi2 = pi2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                       b_pi = b_pi, pi = pi, pi2 = pi2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # sanity check in case of irregularity
  bad_test_cases[idx_b] = length(output) - length(output1)
  
  if (is.null(pvMat_margin)) {
    pvMat_margin = list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] = sapply(output1, function(x) x$pval_1)
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing

Allmethods = c("margin", "mom", "ttest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef <- sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw = df_pw[1:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size Factor (" ~ theta^{(2)}/theta^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)

GraphName = paste0("./newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 


#  ---- POISSON GAMMA VARIANCE n1 = 100; n2 = 200 ----
DataType = "pg"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 200
p = 2
repID = 77 
mu1 = 10
V1 = 15
params1 = gamma_params(mu1, V1)
alpha1 = params1$alpha
beta1 = params1$beta
target_means = 10
target_variance = seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s = target_means^2 / (target_variance - target_means)
beta2s = target_means / (target_variance - target_means)
result_df = data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 = alpha2s[idx]
  beta2 = beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  pvMat_margin = cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_variance - V1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- ZERO-INFLATED NEGATIVE BINOMIAL VARIANCE n1 = 100; n2 = 200----
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 200
sigma_sq = 5
b_thetas = seq(1, 1.70, by =  0.07)
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "dispersion"
repID = 77
K = 100
p = 2
DataType = "zinb"
maxIter = 300 
nCPUS = 8
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(b_thetas))
for (idx_b in 1:length(b_thetas)) {
  b_theta = b_thetas[idx_b]
  print(paste0("scale factor b_theta =", b_theta, " or (variance) = ", b_theta*theta))
  
  
  print("SEF:")
  output = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    tryCatch({
      simulateZINB(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                   b_pi = b_pi, pi = pi, pi2 = pi2,
                   repID, idx = i, n1 = n1, n2 = n2, de_type = de_type, K = NULL, p = p)
    }, error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(NULL)
    })
  }
  output1 = Filter(Negate(is.null), output) 
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    simulateMoMRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                         b_pi = b_pi, pi = pi, pi2 = pi2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID, n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                           b_pi = b_pi, pi = pi, pi2 = pi2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                       b_pi = b_pi, pi = pi, pi2 = pi2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # sanity check in case of irregularity
  bad_test_cases[idx_b] = length(output) - length(output1)
  
  if (is.null(pvMat_margin)) {
    pvMat_margin = list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] <- sapply(output1, function(x) x$pval_1)
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing

Allmethods = c("margin", "mom", "ttest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef <- sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw = df_pw[1:8,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size Factor (" ~ theta^{(2)}/theta^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

########## Fig S.4 Power Analysis ##########
#  ---- POISSON GAMMA F-TEST ----
DataType = "pg"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 2
repID = 77 
mu1 = 10
V1 = 15
params1 = gamma_params(mu1, V1)
alpha1 = params1$alpha
beta1 = params1$beta
target_means = 10
target_variance = seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s = target_means^2 / (target_variance - target_means)
beta2s = target_means / (target_variance - target_means)
result_df = data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ftest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 = alpha2s[idx]
  beta2 = beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    f_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ftest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(f_testPB = f_testPB, ksPB = ksPB)
  }
  
  pvMat_margin = cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom  = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ftest = cbind(pvMat_ftest, sapply(outputPB, function(x) x$f_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ftest", "ks")
df_pw = data.frame(effect_size = (target_variance - V1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ftest = colMeans(pvMat_ftest < q),
                   ks = colMeans(pvMat_ks < q))

df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ftest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)

#  ---- ZERO-INFLATED NEGATIVE BINOMIAL F-TEST ----
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_thetas = seq(1, 1.70, by =  0.07)
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "dispersion"
repID = 77
K = 100
p = 2
DataType = "zinb"
maxIter = 300 
nCPUS = 8
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ftest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' )
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(b_thetas))
for (idx_b in 1:length(b_thetas)) {
  b_theta = b_thetas[idx_b]
  print(paste0("scale factor b_theta =", b_theta, " or (variance) = ", b_theta*theta))
  
  print("SEF:")
  output = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    tryCatch({
      simulateZINB(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                   b_pi = b_pi, pi = pi, pi2 = pi2,
                   repID, idx = i, n1 = n1, n2 = n2, de_type = de_type, K = NULL, p = p)
    }, error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(NULL)
    })
  }
  output1 = Filter(Negate(is.null), output)
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    simulateMoMRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                         b_pi = b_pi, pi = pi, pi2 = pi2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID, n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    f_testPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                           b_pi = b_pi, pi = pi, pi2 = pi2,
                                           distribution = DataType,test_type = "ftest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                       b_pi = b_pi, pi = pi, pi2 = pi2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(f_testPB = f_testPB, ksPB = ksPB)
  }
  
  # sanity check in case of irregularity
  bad_test_cases[idx_b] = length(output) - length(output1)
  
  if (is.null(pvMat_margin)) {
    pvMat_margin = list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] = sapply(output1, function(x) x$pval_1)
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ftest = cbind(pvMat_ftest, sapply(outputPB, function(x) x$f_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing

Allmethods = c("margin", "mom", "ftest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef = sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ftest = colMeans(pvMat_ftest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw = df_pw[1:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ftest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size Factor (" ~ theta^{(2)}/theta^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)

GraphName = paste0("./newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

