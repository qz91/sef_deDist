# generalized version of scDD: uses KS test due to extreme computational constraints, still has issues with FDR control
rm(list = ls()) 
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(scDD)
library(scater)
library(SingleCellExperiment)
library(moments)
source("/home/zaqian/sim_par/lpc_simulation_main_5_8_2026.R")
source("/home/zaqian/sim_par/lpc_main_NB_4_30_2026.R")

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
    ans = tolower(trimws(readline(sprintf(
      "%s: %s [%s]: ",
      prompt,
      paste(choices, collapse = ", "),
      default
    ))))
    
    if (ans == "") ans = default
    if (ans %in% choices) return(ans)
    
    message("Please enter one of: ", paste(choices, collapse = ", "))
  }
}

ask_analysis_type = function() {
  ask_choice("Analysis type", c("shift", "modal"), "shift")
}

ask_shift = function() {
  repeat {
    ans = tolower(trimws(readline("Shift type, either mean or variance: ")))
    if (ans %in% c("mean", "variance")) return(ans)
    message("Please enter only 'mean' or 'variance'.")
  }
}


effect_configs = list(
  standard = list(
    mean_dir = "mean_shift_NB",
    var_dir = "var_shift_NB"
  ),
  small = list(
    mean_dir = "mean_shift_NB_small_effect",
    var_dir = "var_shift_NB_small_effect"
  ),
  strong = list(
    mean_dir = "strong_mean_shift_NB",
    var_dir = "strong_var_shift_NB"
  ),
  new_2_4_fc = list(
    mean_dir = "new_2_4_fc_mean_shift_NB",
    var_dir = "new_2_4_fc_var_shift_NB"
  ),
  new_15_3_fc = list(
    mean_dir = "new_15_3_fc_mean_shift_NB",
    var_dir = "new_15_3_fc_var_shift_NB"
  )
)
effect_choices = names(effect_configs)

ask_effect_size = function() {
  repeat {
    ans = tolower(trimws(readline(sprintf(
      "Effect size: %s [standard]: ",
      paste(effect_choices, collapse = ", ")
    ))))
    
    if (ans == "") return("standard")
    if (ans %in% effect_choices) return(ans)
    
    message(sprintf(
      "Please enter only one of: %s.",
      paste(effect_choices, collapse = ", ")
    ))
  }
}

n_sim = ask_int("How many simulations?", 50)
n_cores = ask_int("How many cores?", 5)
analysis_type = ask_analysis_type()

base_dat_dir = "/home/zaqian/sim_par/generate_sim_obj_seeds"

if (analysis_type == "shift") {
  shift = ask_shift()
  shift_tag = ifelse(shift == "mean", "mean", "var")
  effect_size = ask_effect_size()
  effect_tag = effect_size
  
  effect_config = effect_configs[[effect_size]]
  
  dat_subdir = if (shift == "mean") {
    effect_config$mean_dir
  } else {
    effect_config$var_dir
  }
  
  dat.dir = file.path(base_dat_dir, dat_subdir)
  dat_pattern = "\\.rds$"
  
  run_tag = sprintf("%s_effect_%s", effect_tag, shift_tag)
  seed_file_pattern = "NB_scDD_%s_seed_%02d.rds"
  err_file_pattern = "NB_scDD_%s_seed_%02d_ERROR.txt"
  metrics_file = sprintf("NB_scDD_%s_metrics.RDS", run_tag)
  
  scDD_dir = file.path(
    "/home/zaqian/sim_par",
    sprintf("%s_effect", effect_tag),
    sprintf("%s_scDD_%s", effect_tag, shift_tag)
  )
  
} else {
  dat.dir = file.path(base_dat_dir, "modal_NB")
  dat_pattern = "^dat_NB_modal_seed_[0-9]+\\.rds$"
  
  effect_tag = "modal"
  shift_tag = "modal"
  run_tag = "modal"
  
  seed_file_pattern = "NB_modal_scDD_seed_%02d.rds"
  err_file_pattern = "NB_modal_scDD_seed_%02d_ERROR.txt"
  metrics_file = "NB_modal_scDD_metrics.RDS"
  
  scDD_dir = "/home/zaqian/sim_par/NB_modal/scDD_modal"
}

message("using data directory: ", dat.dir)
stopifnot(dir.exists(dat.dir))

# ----- run scDD once function -----
message(sprintf("beginning %s scDD tests", analysis_type))

register(SnowParam(workers = n_cores))  # start with fewer workers
run_one_scDD = function(dat.name) {
  seed = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.name)))
  if (is.na(seed)) stop(sprintf("could not parse seed from %s", dat.name))
  
  message(sprintf("seed %02d %s data loading for scDD...", seed, analysis_type))
  sim_obj = readRDS(dat.name)
  flush.console()
  
  gene_info = sim_obj$gene_info
  
  metadata = sim_obj$metadata
  rownames(metadata) = metadata$cell_id
  metadata$cell_id = NULL
  metadata$group = factor(metadata$group)
  
  sce = SingleCellExperiment::SingleCellExperiment(
    assays = list(normcounts = t(sim_obj$count_matrix)),
    colData = metadata
  )
  
  rm(sim_obj, metadata)
  gc()
  
  prior_param = list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01) # if using DP test
  scDD_model = scDD(sce, testZeroes = FALSE, condition = "group", categorize = FALSE)
  scDD_res = results(scDD_model)
  
  labels = gene_info$de_status[match(scDD_res$gene, gene_info$gene)]
  metrics = get_FDR_TPR_metrics(
    pvals = scDD_res$nonzero.pvalue.adj,
    labels = labels
  )
  
  out = c(
    list(seed = seed),
    as.list(metrics)
  )
  
  return(out)
}


# ----- load data file names -----
dat.files = list.files(
  dat.dir,
  pattern = dat_pattern,
  full.names = TRUE
)
dat.seeds = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.files)))
stopifnot(!any(is.na(dat.seeds)))
dat.files = dat.files[order(dat.seeds)]
dat.seeds = dat.seeds[order(dat.seeds)]

# need to do this because we aren't looking at all 50 datasets necessarily
keep = dat.seeds %in% seq_len(n_sim)

dat.files = dat.files[keep]
dat.seeds = dat.seeds[keep]

stopifnot(length(dat.files) == n_sim)
stopifnot(all(dat.seeds == seq_len(n_sim)))

message(sprintf("found %02d pre-generated NB data files", length(dat.files)))
flush.console()


# ----- output directories -----
dir.create(scDD_dir, recursive = TRUE, showWarnings = FALSE)
scDD_fname = metrics_file

# ----- run scDD in parallel for selected/missing seeds -----
message(sprintf("checking scDD outputs for seeds 1 to %02d", n_sim))
flush.console()

if (analysis_type == "shift") {
  scDD_paths = file.path(
    scDD_dir,
    sprintf(seed_file_pattern, run_tag, dat.seeds)
  )
} else {
  scDD_paths = file.path(
    scDD_dir,
    sprintf(seed_file_pattern, dat.seeds)
  )
}

missing = !file.exists(scDD_paths) # checks for nonexistent seed results

if (any(missing)) { # avoids running seeds that have been compiled before
  message("running missing scDD seeds: ", paste(dat.seeds[missing], collapse = ", "))
} else {
  message("all requested scDD seed outputs already exist")
}

for (i in which(missing)) { # use for loop to look for corrupted cases
  dat.name = dat.files[i]
  seed = dat.seeds[i]
  path = scDD_paths[i]
  
  if (analysis_type == "shift") {
    err_path = file.path(
      scDD_dir,
      sprintf(err_file_pattern, run_tag, seed)
    )
  } else {
    err_path = file.path(
      scDD_dir,
      sprintf(err_file_pattern, seed)
    )
  }
  
  out = tryCatch({
    out = run_one_scDD(dat.name) # running here
    
    tmp_path = paste0(path, ".tmp")
    saveRDS(out, tmp_path)
    file.rename(tmp_path, path)
    
    if (file.exists(err_path)) file.remove(err_path)
    
    out
  }, error = function(e) {
    msg = sprintf("seed %02d scDD failed: %s", seed, conditionMessage(e))
    message(msg)
    writeLines(msg, err_path)
    NULL
  })
  
  rm(out)
  gc()
}

# ----- combine completed seeds for scDD results -----
results_scDD_list = lapply(seq_along(scDD_paths), function(i) {
  seed = dat.seeds[i]
  path = scDD_paths[i]
  
  if (!file.exists(path)) {
    message(sprintf("seed %02d scDD output missing, skipping in combined results", seed))
    return(NULL)
  }
  
  out = tryCatch(readRDS(path), error = function(e) {
    message(sprintf("seed %02d scDD output unreadable, skipping in combined results", seed))
    NULL
  })
  
  out
})

names(results_scDD_list) = sprintf("seed_%02d", dat.seeds)

failed = vapply(results_scDD_list, is.null, logical(1))

if (any(failed)) {
  message("not included in combined results: ", paste(dat.seeds[failed], collapse = ", "))
}

results_scDD_ok = results_scDD_list[!failed]

results_scDD = do.call(
  rbind,
  lapply(results_scDD_ok, function(x) as.data.frame(as.list(x)))
)

rownames(results_scDD) = NULL
results_scDD = cbind(method = "scDD", results_scDD)

saveRDS(
  results_scDD,
  file.path(scDD_dir, scDD_fname)
)