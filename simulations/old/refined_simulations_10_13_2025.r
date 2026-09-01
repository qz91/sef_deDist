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

########## QQ Plots ##########
#  ---- POISSON GAMMA  ----
nCPUS = 10
DataType = "pg"
maxIter = 300 
n1 = 100; n2 = 100
p = 2
alpha2 = 20; alpha1 = 20
beta1 = 2; beta2 = 2
repID = 77  
pvVec_margin = NULL 
output <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    simulatePGRealistic(
      repID = repID,
      idx = i,
      alpha1 = alpha1, alpha2 = alpha2,
      beta1 = beta1, beta2 = beta2,
      n1 = n1, n2 = n2,
      de_type = "mean", p = p
    )
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", e$message))
    NULL
  })
}, mc.cores = 4) #parallel::detectCores() - 1
pvVec_margin <- sapply(output, function(x) x$pval_1)
gg_qqplot(unlist(pvVec_margin))
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Poisson Gamma Model")

#  ---- ZERO-INFLATED NEGATIVE BINOMIAL  ----
nCPUS = 10
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
p = 2

output <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    simulateZINB(
      repID = repID, idx = i,
      mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
      theta = theta, b_theta = 1,
      pi = pi, pi2 = pi2,
      n1 = n1, n2 = n2,
      de_type = de_type, p = p
    )
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
    return(NULL)
  })
}, mc.cores = 4)
pvVec_margin = sapply(output, function(x) x$pval_1)
gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB Model")
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB Model")
print(realistic_qqplot)



########## POWER ANALYSIS ##########
#  ---- POISSON GAMMA MEAN ----
nCPUS = 10
DataType = "pg"
de_type = "mean"
maxIter = 300 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 2
repID = 77 

mu1 = 10
V1 = 15
params1 = gamma_params(mu1, V1)
alpha1 = params1$alpha
beta1  = params1$beta
target_means = seq(10, 12, length.out = 15)
target_variance = 15 # must be > any target mean
if(any(target_means >= target_variance)) {
  stop("Each target mean must be less than the target variance.")
}
alpha2s <- target_means^2 / (target_variance - target_means)
beta2s  <- target_means / (target_variance - target_means)
result_df <- data.frame(
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
print(output$binwidth)
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_means - mu1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw = df_pw[1:10,]
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
  labs(x = expression("Effect Size (" ~ mu^{(2)} - mu^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./Powerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- POISSON GAMMA VARIANCE ----
nCPUS = 10
DataType = "pg"
de_type = "variance"
maxIter = 300 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 2
repID = 77 
mu1 = 10
V1 = 15
params1 = gamma_params(mu1, V1)
alpha1 = params1$alpha
beta1  = params1$beta
target_means = 10
target_variance = seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s <- target_means^2 / (target_variance - target_means)
beta2s  <- target_means / (target_variance - target_means)
result_df <- data.frame(
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
  alpha2 <- alpha2s[idx]
  beta2 <- beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
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
df_pw = df_pw[1:10,]
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
GraphName = paste0("./Powerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL MEAN  ----
nCPUS = 10
mu = 10
mu2s = seq(10, 13, 0.25)
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_theta = 1
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "mean"
repID = 77
K = 100
p = 2
DataType = "zinb"
maxIter = 300 
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' ) 
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(mu2s))
for (idx_mu2 in 1:length(mu2s)) {
  mu2 = mu2s[idx_mu2]
  print(paste0("mu2: ", mu2))
  
  
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
  bad_test_cases[idx_mu2] = length(output) - length(output1)
  
  if (is.null(pvMat_margin)) {
    pvMat_margin = list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] = sapply(output1, function(x) x$pval_1)
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing
Allmethods = c("margin", "mom", "ks","ttest")
q = 0.05 
margin_sef = sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (mu2s - mu), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))

df_pw = df_pw[1:10,]
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
  labs(x = expression("Effect Size (" ~ mu[2] - mu[1] ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./Powerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL VARIANCE  ----
nCPUS = 10
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
  pvMat_margin[[length(pvMat_margin) + 1]] = sapply(output1, function(x) x$pval_1)
  pvMat_mom = cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest = cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks = cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) 

Allmethods = c("margin", "mom", "ttest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef = sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))

df_pw = df_pw[1:10,]
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

GraphName = paste0("./Powerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 
