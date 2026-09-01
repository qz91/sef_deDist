rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)

source("/home/zaqian/sim_par/lpc_NB_modality_analysis_final.R")

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

n_cores = ask_int("How many cores?", 20)
n_sim = ask_int("How many simulations?", 50)

out.dir = "/home/zaqian/sim_par/generate_sim_obj_seeds/modal_NB"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# ----- params -----
n1 = 100; n2 = 100
lower = 500
upper = 1000

n_null_unimodal = 1000
n_null_bimodal = 1000
n_modal = 200

target_mu_range = c(100, 100)
delta_range = c(30, 30)
target_marginal_variance = 1150

theta_mix_range = c(100, 100)
mixing_prob = 0.5
sigma_sq = 0.01

modal_unimodal_group = "group1"
return_donor_effects = TRUE
message("beginning modality-only NB data generation")

generate_data = function(seed) {
  sim_obj = simulate_NB_modality_analysis(
    seed = seed,
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
    return_donor_effects = return_donor_effects
  )
  
  message(sprintf("seed %02d modality data generation complete", seed))
  return(sim_obj)
}

save_one_seed = function(seed = 1) {
  sim_obj = generate_data(seed = seed)
  
  out.file = file.path(
    out.dir,
    sprintf("dat_NB_modal_seed_%02d.rds", seed)
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