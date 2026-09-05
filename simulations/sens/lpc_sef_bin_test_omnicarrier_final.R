rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)
source("lpc_simulation_main_5_8_2026_final.R")
source("lpc_main_NB_4_30_2026_final.R")
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

ask_K = function(prompt, default = "FD") {
  repeat {
    ans = trimws(readline(sprintf("%s [%s]: ", prompt, default)))
    
    if (ans == "") ans = default
    
    if (toupper(ans) == "FD") {
      return(NULL)
    }
    
    K_value = suppressWarnings(as.integer(ans))
    
    if (!is.na(K_value) && K_value >= 10) {
      return(K_value)
    }
    
    message("Please enter 'FD' or an integer greater than or equal to 10.")
  }
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
  ask_choice(
    "Data size/effect: small, strong, new_2_4_fc, or new_15_3_fc",
    c("small", "strong", "new_2_4_fc", "new_15_3_fc"),
    "small"
  )
}

ask_carrier = function() {
  ask_choice("Carrier density type: discrete or cts", c("discrete", "cts"), "cts")
}

data_size = ask_data_size()
carrier_type = ask_carrier()
use_discrete = carrier_type == "discrete"

K = ask_K("Bin resolution: FD or a bin count", "FD")
if (use_discrete && is.null(K)) {
  message("Discrete carrier does not use FD; using full observed integer support.")
}

K_label = if (is.null(K)) {
  if (use_discrete) "observed_support" else "FD"
} else {
  as.character(K)
}


by_shift = function(mean, variance) list(mean = mean, variance = variance)

mk_case = function(pretty, dat_mean, dat_var,
                   scaled_mean_dir, scaled_var_dir, scaled_prefix) {
  list(
    pretty = pretty,
    dat = by_shift(dat_mean, dat_var),
    scaled_dir = by_shift(scaled_mean_dir, scaled_var_dir),
    scaled_seed = by_shift(
      paste0(scaled_prefix, "_mean_seed_%02d.rds"),
      paste0(scaled_prefix, "_var_seed_%02d.rds")
    ),
    scaled_metrics = by_shift(
      paste0(scaled_prefix, "_mean_sef_metrics_list.RDS"),
      paste0(scaled_prefix, "_var_sef_metrics_list.RDS")
    )
  )
}

cases = list(
  small = mk_case(
    "small effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/mean_shift_NB_small_effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/var_shift_NB_small_effect",
    
    "/home/zaqian/sim_par/ls_small_effect/ls_scaled_small_effect_mean",
    "/home/zaqian/sim_par/ls_small_effect/ls_scaled_small_effect_var",
    "NB_scaled_small_effect"
  ),
  
  strong = mk_case(
    "strong effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/strong_mean_shift_NB",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/strong_var_shift_NB",

    "/home/zaqian/sim_par/strong_effect/strong_sef_scaled_mean",
    "/home/zaqian/sim_par/strong_effect/strong_sef_scaled",
    "NB_scaled_strong_effect"
  ),
  
  new_2_4_fc = mk_case(
    "new_2_4_fc effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_2_4_fc_mean_shift_NB",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_2_4_fc_var_shift_NB",

    "/home/zaqian/sim_par/strong_effect/new_2_4_fc_scaled_mean",
    "/home/zaqian/sim_par/strong_effect/new_2_4_fc_sef_scaled",
    "NB_scaled_new_2_4_fc_effect"
  ),
  
  new_15_3_fc = mk_case(
    "new_15_3_fc effect",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_15_3_fc_mean_shift_NB",
    "/home/zaqian/sim_par/generate_sim_obj_seeds/new_15_3_fc_var_shift_NB",

    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_scaled_mean",
    "/home/zaqian/sim_par/new_15_3_fc_effect/new_15_3_fc_sef_scaled",
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

scaled_seed_file = add_discrete_to_seed_file(cfg$scaled_seed[[shift]]) 
scaled_metrics_file = paste0("K_", K_label, "_", add_discrete_to_metrics_file(cfg$scaled_metrics[[shift]]))

dat.dir = cfg$dat[[shift]]
out.dir.scaled = file.path(
  add_discrete_folder(cfg$scaled_dir[[shift]]),
  sprintf("K_%s", K_label)
)
stopifnot(dir.exists(dat.dir))

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

run_one_scaled_sef_preload = function(dat.name = NULL, ridge_run = 1e-8) {
  seed = seed_from_path(dat.name)
  message(sprintf("seed %02d NB data loading...", seed))
  sim_obj = readRDS(dat.name)
  
  ps = c(2, 3, 4)
  alphas = c(0.05, 0.10)
  
  results_by_p = setNames(lapply(ps, function(p) {
    test = scaled_sef_fun(
      sim_obj = sim_obj,
      K = K,
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
  
  list(
    seed = seed,
    data_size = data_size,
    carrier_type = carrier_type,
    K = K,
    K_label = K_label,
    ridge = ridge_run,
    results = results_by_p
  )
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


scaled_metrics_list = run_preloaded_set(
  out.dir.scaled,
  scaled_seed_file,
  scaled_metrics_file,
  run_one_scaled_sef_preload,
  sprintf("scaled SEF (K = %s)", K_label)
)

message("done")

