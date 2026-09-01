rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)
source("/home/zaqian/sim_par/lpc_simulation_main_5_8_2026_final.R")
source("/home/zaqian/sim_par/lpc_main_NB_4_30_2026_final.R")
# generate sim_obj data so we don't have to regenerate every simulation setting
n_cores = 20
n_sim = 50

ask_int = function(prompt, default) {
  ans = trimws(readline(sprintf("%s [%d]: ", prompt, default)))
  if (ans == "") return(default)
  
  ans = suppressWarnings(as.integer(ans))
  if (is.na(ans) || ans < 1) {
    stop(sprintf("%s must be a positive integer.", prompt))
  }
  
  ans
}

ask_shift = function() {
  repeat {
    ans = tolower(trimws(readline("Shift type, either mean or variance: ")))
    if (ans %in% c("mean", "variance")) return(ans)
    message("Please enter only 'mean' or 'variance'.")
  }
}

n_cores = ask_int("How many cores?", 20)
n_sim = ask_int("How many simulations?", 50)
shift = ask_shift()

shift_tag = ifelse(shift == "mean", "mean", "var")
out.dir = ifelse(
  shift == "mean",
  "/home/zaqian/sim_par/generate_sim_obj_seeds/mean_shift_NB_small_effect",
  "/home/zaqian/sim_par/generate_sim_obj_seeds/var_shift_NB_small_effect"
)

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)


# ----- params -----
n1 = 100; n2 = 100
lower = 500
upper = 1000
n_nonDE_unimodal = 1000 
n_nonDE_mixture = 1000 
n_DE = 200 
nonDE_unimodal_mu_range = c(5, 20)
nonDE_unimodal_theta_range = c(3, 10)
nonDE_mixture_mu_range = c(5, 15)
nonDE_mixture_theta_range = c(3, 10)
nonDE_mixture_fc_range = c(2, 4)
nonDE_mixture_prob_range = c(0.30, 0.70)
DE_mu_range = c(5, 15)
DE_theta_range = c(3, 10)
DE_mixture_fc_range = c(1.25, 2.25)
sigma_sq_nonDE = 0.05
sigma_sq = 0.05 # sets them equal
sigma_sq_DE = 0.05
de_unimodal_group = "group1"
share_mode = F # DE unimodal case in group 1 has mean = (mu1 + mu2)/2 for the mixture DE

if (shift == "mean") {
  DE_mixture_prob_range = c(0.30, 0.70)
} else if (shift == "variance") {
  DE_mixture_prob_range = c(0.5, 0.5)
}
# ----- generate data -----
message(paste0("beginning data generation with shift type: ",shift)) 


generate_data = function(seed){
  sim_obj = generate_NB_gene_data(
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    n_nonDE_unimodal = n_nonDE_unimodal, 
    n_nonDE_mixture = n_nonDE_mixture, 
    n_DE = n_DE, 
    nonDE_unimodal_mu_range = nonDE_unimodal_mu_range, nonDE_unimodal_theta_range = nonDE_unimodal_theta_range,
    nonDE_mixture_mu_range = nonDE_mixture_mu_range, nonDE_mixture_theta_range =nonDE_mixture_theta_range, 
    nonDE_mixture_fc_range = nonDE_mixture_fc_range, nonDE_mixture_prob_range = nonDE_mixture_prob_range,
    DE_mu_range = DE_mu_range,
    DE_theta_range = DE_theta_range,
    DE_mixture_fc_range = DE_mixture_fc_range,
    DE_mixture_prob_range = DE_mixture_prob_range,
    sigma_sq_nonDE = sigma_sq_nonDE,
    sigma_sq_DE = sigma_sq_DE,
    sigma_sq = sigma_sq,
    de_unimodal_group = de_unimodal_group,
    share_mode = share_mode,
    seed = seed
  )
  message(sprintf("seed %02d data generation complete", seed))
  return(sim_obj)
}

save_one_seed = function(seed = 1) {
  sim_obj = generate_data(seed = seed)
  
  out.file = file.path(
    out.dir,
    sprintf("dat_NB_%s_small_effect_seed_%02d.rds", shift_tag, seed)
  )
  
  saveRDS(sim_obj, out.file)
  message(sprintf("seed %02d saved to %s", seed, out.file))
  
  invisible(out.file)
}

saved_files = pbmcapply::pbmclapply(
  X = seq_len(n_sim),
  FUN = save_one_seed,
  mc.cores = n_cores
)


