# used for the cluster to run SEF sims in parallel

rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)
source("lpc_simulation_main_5_8_2026_final.R")
source("lpc_revisions_simulation_discrete_6_7_2026_final.R")
n_cores = 25
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

n_cores = ask_int("How many cores?", 25)
shift = ask_shift()
shift_tag = ifelse(shift == "mean", "mean", "var")

ask_choice = function(prompt, choices, default) {
  repeat {
    ans = tolower(trimws(readline(sprintf("%s [%s]: ", prompt, default))))
    if (ans == "") ans = default
    if (ans %in% choices) return(ans)
    message("Please enter one of: ", paste(choices, collapse = ", "))
  }
}

ask_data_size = function() {
  ans = ask_choice(
    "Data size/effect: small, strong, std, or new_15_3_fc",
    c("small", "strong", "std", "standard", "new_15_3_fc"),
    "std"
  )
  if (ans == "standard") ans = "std"
  ans
}

ask_carrier = function() {
  ask_choice("Carrier density type: discrete or cts", c("discrete", "cts"), "cts")
}

data_size = ask_data_size()
carrier_type = ask_carrier()
use_discrete = carrier_type == "discrete"

by_shift = function(mean, variance) list(mean = mean, variance = variance)

mk_case = function(pretty, dat_mean, dat_var,
                   raw_mean_dir, raw_var_dir,
                   scaled_mean_dir, scaled_var_dir,
                   raw_prefix, scaled_prefix) {
  list(
    pretty = pretty,
    dat = by_shift(dat_mean, dat_var),
    raw_dir = by_shift(raw_mean_dir, raw_var_dir),
    scaled_dir = by_shift(scaled_mean_dir, scaled_var_dir),
    raw_seed = by_shift(
      paste0(raw_prefix, "_mean_seed_%02d.rds"),
      paste0(raw_prefix, "_var_seed_%02d.rds")
    ),
    scaled_seed = by_shift(
      paste0(scaled_prefix, "_mean_seed_%02d.rds"),
      paste0(scaled_prefix, "_var_seed_%02d.rds")
    ),
    raw_metrics = by_shift(
      paste0(raw_prefix, "_mean_sef_metrics_list.RDS"),
      paste0(raw_prefix, "_var_sef_metrics_list.RDS")
    ),
    scaled_metrics = by_shift(
      paste0(scaled_prefix, "_mean_sef_metrics_list.RDS"),
      paste0(scaled_prefix, "_var_sef_metrics_list.RDS")
    )
  )
}

cases = list(
  std = list(
    pretty = "standard effect",
    dat = by_shift(
      "/home/zaqian/sim_par/generate_sim_obj_seeds/mean_shift_NB",
      "/home/zaqian/sim_par/generate_sim_obj_seeds/var_shift_NB"
    ),
    raw_dir = by_shift(
      "/home/zaqian/sim_par/ls_sef_mean",
      "/home/zaqian/sim_par/ls_sef"
    ),
    scaled_dir = by_shift(
      "/home/zaqian/sim_par/ls_sef_scaled_mean",
      "/home/zaqian/sim_par/ls_sef_scaled"
    ),
    raw_seed = by_shift("ls_NB_mean_seed_%02d.rds", "ls_NB_large_seed_%02d.rds"),
    scaled_seed = by_shift("ls_NB_scaled_mean_seed_%02d.rds", "ls_NB_large_seed_%02d.rds"),
    raw_metrics = by_shift("ls_NB_sef_mean_metrics_list.RDS", "ls_NB_large_sef_metrics_list.RDS"),
    scaled_metrics = by_shift("ls_NB_scaled_mean_sef_metrics_list.RDS", "ls_NB_large_scaled_sef_metrics_list.RDS")
  ),
  
  small = mk_case(
    "small effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/mean_shift_NB_small_effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/var_shift_NB_small_effect",
    "/home/zaqian/sim_par/ls_small_effect/ls_small_effect_mean",
    "/home/zaqian/sim_par/ls_small_effect/ls_small_effect_var",
    "/home/zaqian/sim_par/ls_small_effect/ls_scaled_small_effect_mean",
    "/home/zaqian/sim_par/ls_small_effect/ls_scaled_small_effect_var",
    "ls_NB_small_effect",
    "ls_NB_scaled_small_effect"
  ),
  
  strong = mk_case(
    "strong effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/strong_mean_shift_NB",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/strong_var_shift_NB",
    "/home/zaqian/sim_par/strong_effect/strong_sef_mean",
    "/home/zaqian/sim_par/strong_effect/strong_sef",
    "/home/zaqian/sim_par/strong_effect/strong_sef_scaled_mean",
    "/home/zaqian/sim_par/strong_effect/strong_sef_scaled",
    "NB_strong_effect",
    "NB_scaled_strong_effect"
  ),
  
  new_15_3_fc = mk_case(
    "new_15_3_fc effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_15_3_fc_mean_shift_NB",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_15_3_fc_var_shift_NB",
    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_sef_mean",
    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_sef",
    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_scaled_mean",
    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_sef_scaled",
    "NB_new_15_3_fc_effect",
    "NB_scaled_new_15_3_fc_effect"
  )
)

cfg = cases[[data_size]]

add_discrete_folder = function(path) {
  if (!use_discrete) return(path)
  
  dir_name = basename(path)
  
  if (!grepl("_(mean|var)$", dir_name)) {
    dir_name = paste0(dir_name, "_", shift_tag)
  }
  
  file.path("/home/zaqian/sim_par/discrete", paste0(dir_name, "_discrete"))
}

add_discrete_to_seed_file = function(x) {
  if (!use_discrete) return(x)
  
  if (grepl("_(mean|var)_seed_%02d\\.rds$", x)) {
    return(sub("_(mean|var)_seed_%02d\\.rds$", "_\\1_discrete_seed_%02d.rds", x))
  }
  
  sub("_seed_%02d\\.rds$", paste0("_", shift_tag, "_discrete_seed_%02d.rds"), x)
}

add_discrete_to_metrics_file = function(x) {
  if (!use_discrete) return(x)
  
  if (grepl("_(mean|var)_sef_metrics_list\\.RDS$", x)) {
    return(sub("_(mean|var)_sef_metrics_list\\.RDS$", "_\\1_discrete_sef_metrics_list.RDS", x))
  }
  
  if (grepl("_sef_(mean|var)_metrics_list\\.RDS$", x)) {
    return(sub("_sef_(mean|var)_metrics_list\\.RDS$", "_\\1_discrete_sef_metrics_list.RDS", x))
  }
  
  sub("_sef_metrics_list\\.RDS$", paste0("_", shift_tag, "_discrete_sef_metrics_list.RDS"), x)
}

raw_seed_file = add_discrete_to_seed_file(cfg$raw_seed[[shift]])
scaled_seed_file = add_discrete_to_seed_file(cfg$scaled_seed[[shift]]) 
raw_metrics_file = add_discrete_to_metrics_file(cfg$raw_metrics[[shift]])
scaled_metrics_file = add_discrete_to_metrics_file(cfg$scaled_metrics[[shift]])

dat.dir = cfg$dat[[shift]]
out.dir.sef = add_discrete_folder(cfg$raw_dir[[shift]])
out.dir.scaled = add_discrete_folder(cfg$scaled_dir[[shift]])
stopifnot(dir.exists(dat.dir))

raw_sef_fun = if (use_discrete) run_sef_test_discrete_covfix else run_sef_test_covfix
scaled_sef_fun = if (use_discrete) run_sef_test_discrete_covfix_scaled_basis else run_sef_test_covfix_scaled_basis
bw_run = if (use_discrete) NULL else 2

message(sprintf(
  "beginning %s %s shift SEF tests using %s carrier",
  cfg$pretty, shift, carrier_type
))

seed_from_path = function(dat.name) {
  as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.name)))
}

make_metrics = function(p_adj, de_status, alphas) {
  setNames(lapply(alphas, function(alpha) {
    m = get_FDR_TPR_metrics(p_adj, de_status, alpha = alpha)
    m$FDP = m$FDR
    m
  }), sprintf("alpha_%0.2f", alphas))
}

run_one_sef_preload = function(dat.name = NULL) {
  seed = seed_from_path(dat.name)
  message(sprintf("seed %02d NB data loading...", seed))
  sim_obj = readRDS(dat.name)
  
  ps = c(2, 3)
  alphas = c(0.05, 0.10)
  
  results_by_p = setNames(lapply(ps, function(p) {
    test = raw_sef_fun(sim_obj = sim_obj, p = p, verbose = FALSE, bw = bw_run)
    test$p_adj = p.adjust(test$pVec, method = "BH")
    
    metrics = make_metrics(test$p_adj, test$de_status, alphas)
    message(sprintf("seed %02d finished p = %d", seed, p))
    
    list(test = test, metrics = metrics)
  }), paste0("p", ps))
  
  list(seed = seed, data_size = data_size, carrier_type = carrier_type, results = results_by_p)
}

run_one_scaled_sef_preload = function(dat.name = NULL, ridge_run = 1e-8) {
  seed = seed_from_path(dat.name)
  message(sprintf("seed %02d NB data loading...", seed))
  sim_obj = readRDS(dat.name)
  
  ps = c(2, 3, 4)
  alphas = c(0.05, 0.10)
  
  results_by_p = setNames(lapply(ps, function(p) {
    test = scaled_sef_fun(
      sim_obj = sim_obj,
      p = p,
      verbose = FALSE,
      bw = bw_run,
      ridge = ridge_run
    )
    
    test$failed = is.na(test$pVec)
    p_for_bh = test$pVec
    p_for_bh[test$failed] = 1
    test$p_adj = p.adjust(p_for_bh, method = "BH")
    
    metrics = make_metrics(test$p_adj, test$de_status, alphas)
    
    message(sprintf(
      "seed %02d finished p = %d scaled basis, failed genes = %d",
      seed, p, sum(test$failed)
    ))
    
    list(test = test, metrics = metrics)
  }), paste0("p", ps))
  
  list(seed = seed, data_size = data_size, carrier_type = carrier_type, ridge = ridge_run, results = results_by_p)
}


dat.files = list.files(dat.dir, pattern = "\\.rds$", full.names = TRUE)
dat.seeds = seed_from_path(dat.files)

stopifnot(!any(is.na(dat.seeds)))

ord = order(dat.seeds)
dat.files = dat.files[ord]
dat.seeds = dat.seeds[ord]

stopifnot(length(dat.files) == n_sim)

message(sprintf("found %02d pre-generated NB data files", length(dat.files)))
flush.console()

run_preloaded_set = function(out.dir, seed_file, metrics_file, runner, label) {
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  
  message(sprintf("running %02d preloaded %s files with %02d cores", n_sim, label, n_cores))
  flush.console()
  
  metrics_list = pbmcapply::pbmclapply(
    dat.files,
    function(dat.name) {
      seed = seed_from_path(dat.name)
      path = file.path(out.dir, sprintf(seed_file, seed))
      
      if (file.exists(path)) {
        message(sprintf("seed %02d %s already exists, loading...", seed, label))
        out = readRDS(path)
      } else {
        out = runner(dat.name)
        saveRDS(out, path)
      }
      
      lapply(out$results, `[[`, "metrics")
    },
    mc.cores = n_cores
  )
  
  saveRDS(metrics_list, file.path(out.dir, metrics_file))
  metrics_list
}

metrics_list = run_preloaded_set(
  out.dir.sef,
  raw_seed_file,
  raw_metrics_file,
  run_one_sef_preload,
  "regular SEF"
)

scaled_metrics_list = run_preloaded_set(
  out.dir.scaled,
  scaled_seed_file,
  scaled_metrics_file,
  run_one_scaled_sef_preload,
  "scaled SEF"
)


message("done")

