rm(list = ls())
library(pbmcapply)
library(BiocParallel)
library(ggplot2)
library(moments)
source("/home/zaqian/sim_par/lpc_simulation_main_5_8_2026.R")
source("/home/zaqian/sim_par/lpc_NB_modality_analysis.R")
source("/home/zaqian/sim_par/lpc_revisions_simulation_discrete_6_7_2026.R")

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

ask_choice = function(prompt, choices, default) {
  repeat {
    ans = tolower(trimws(readline(sprintf("%s [%s]: ", prompt, default))))
    if (ans == "") ans = default
    if (ans %in% choices) return(ans)
    message("Please enter one of: ", paste(choices, collapse = ", "))
  }
}

ask_carrier = function() {
  ask_choice("Carrier density type: discrete or cts", c("discrete", "cts"), "cts")
}

n_cores = ask_int("How many cores?", 25)
carrier_type = ask_carrier()
use_discrete = carrier_type == "discrete"

dat.dir = "/home/zaqian/sim_par/generate_sim_obj_seeds/modal_NB"
out.dir.scaled = if (use_discrete) {
  "/home/zaqian/sim_par/NB_modal/sef_scaled_modal_discrete"
} else {
  "/home/zaqian/sim_par/NB_modal/sef_scaled_modal"
}

scaled_seed_file = if (use_discrete) {
  "NB_modal_sef_scaled_modal_discrete_seed_%02d.rds"
} else {
  "NB_modal_sef_scaled_modal_seed_%02d.rds"
}
scaled_metrics_file = if (use_discrete) {
  "NB_modal_sef_scaled_modal_discrete_metrics_list.RDS"
} else {
  "NB_modal_sef_scaled_modal_metrics_list.RDS"
}

scaled_sef_fun = if (use_discrete) {
  run_sef_test_discrete_covfix_scaled_basis
} else {
  run_sef_test_covfix_scaled_basis
}
bw_run = if (use_discrete) NULL else 2
ps = c(2, 3, 4, 5)
alphas = c(0.05, 0.10)

stopifnot(dir.exists(dat.dir))

message(sprintf("beginning modality-only scaled SEF tests using %s carrier", carrier_type))
message(sprintf("data directory: %s", dat.dir))
message(sprintf("output directory: %s", out.dir.scaled))
message(sprintf("p values: %s", paste(ps, collapse = ", ")))

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
  message(sprintf("seed %02d modal data loading...", seed))
  sim_obj = readRDS(dat.name)

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

  list(
    seed = seed,
    analysis = "modality",
    carrier_type = carrier_type,
    ridge = ridge_run,
    p_values = ps,
    results = results_by_p
  )
}

dat.files = list.files(
  dat.dir,
  pattern = "^dat_NB_modal_seed_[0-9]+\\.rds$",
  full.names = TRUE
)
dat.seeds = seed_from_path(dat.files)

stopifnot(!any(is.na(dat.seeds)))

ord = order(dat.seeds)
dat.files = dat.files[ord]
dat.seeds = dat.seeds[ord]

stopifnot(length(dat.files) == n_sim)

message(sprintf("found %02d pre-generated NB modality data files", length(dat.files)))
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
  "scaled modality SEF"
)

message("done")
