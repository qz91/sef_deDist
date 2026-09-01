# lpc_NB_muscat_generalized
rm(list = ls())
library(pbmcapply)
library(pbapply)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(variancePartition)
library(muscat)
library(SingleCellExperiment)
source("/home/zaqian/sim_par/lpc_simulation_main_5_8_2026_final.R")

# ----- prompts -----
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

ask_ddf = function() {
  repeat {
    ans = trimws(readline(
      "DDF method: Satterthwaite or KR (Kenward-Roger) [Satterthwaite]: "
    ))
    
    if (ans == "") return("Satterthwaite")
    
    ans_lower = tolower(ans)
    
    if (ans_lower %in% c("satterthwaite", "satt")) {
      return("Satterthwaite")
    }
    
    if (ans_lower %in% c("kr", "kenward-roger", "kenward roger")) {
      return("Kenward-Roger")
    }
    
    message("Please enter Satterthwaite or KR / Kenward-Roger.")
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

ddf_choice = ask_ddf()
ddf_tag = ifelse(ddf_choice == "Kenward-Roger", "KR", "Satterthwaite")
analysis_type = ask_analysis_type()

# ----- load data files -----
base_dat_dir = "/home/zaqian/sim_par/generate_sim_obj_seeds"

if (analysis_type == "shift") {
  effect_size = ask_effect_size()
  effect_tag = effect_size
  
  shift = ask_shift()
  shift_tag = ifelse(shift == "mean", "mean", "var")
  
  effect_config = effect_configs[[effect_size]]
  
  dat_subdir = if (shift == "mean") {
    effect_config$mean_dir
  } else {
    effect_config$var_dir
  }
  
  dat.dir = file.path(base_dat_dir, dat_subdir)
  dat_pattern = "\\.rds$"
  
  muscat_dir = file.path(
    "/home/zaqian/sim_par",
    sprintf("%s_effect", effect_tag),
    sprintf("%s_muscat_%s_ddf_%s", effect_tag, shift_tag, ddf_tag)
  )
  
  muscat_fname = sprintf(
    "NB_muscat_%s_effect_%s_ddf_%s_metrics.RDS",
    effect_tag,
    shift_tag,
    ddf_tag
  )
  
  muscat_seed_file = function(seed) {
    sprintf(
      "NB_muscat_%s_effect_%s_seed_%02d_ddf_%s.rds",
      effect_tag,
      shift_tag,
      seed,
      ddf_tag
    )
  }
  
  muscat_error_file = function(seed) {
    sprintf(
      "NB_muscat_%s_effect_%s_seed_%02d_ddf_%s_ERROR.txt",
      effect_tag,
      shift_tag,
      seed,
      ddf_tag
    )
  }
  
} else {
  dat.dir = file.path(base_dat_dir, "modal_NB")
  dat_pattern = "^dat_NB_modal_seed_[0-9]+\\.rds$"
  
  muscat_dir = file.path(
    "/home/zaqian/sim_par/NB_modal",
    sprintf("muscat_modal_ddf_%s", ddf_tag)
  )
  
  muscat_fname = sprintf(
    "NB_modal_muscat_ddf_%s_metrics.RDS",
    ddf_tag
  )
  
  muscat_seed_file = function(seed) {
    sprintf("NB_modal_muscat_seed_%02d_ddf_%s.rds", seed, ddf_tag)
  }
  
  muscat_error_file = function(seed) {
    sprintf("NB_modal_muscat_seed_%02d_ddf_%s_ERROR.txt", seed, ddf_tag)
  }
}

message("using data directory: ", dat.dir)
stopifnot(dir.exists(dat.dir))

dat.files = list.files(
  dat.dir,
  pattern = dat_pattern,
  full.names = TRUE
)

dat.seeds = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.files)))
stopifnot(!any(is.na(dat.seeds)))

dat.files = dat.files[order(dat.seeds)]
dat.seeds = dat.seeds[order(dat.seeds)]

keep = dat.seeds %in% seq_len(n_sim)

dat.files = dat.files[keep]
dat.seeds = dat.seeds[keep]

stopifnot(length(dat.files) == n_sim)
stopifnot(all(dat.seeds == seq_len(n_sim)))

message(sprintf("found %02d pre-generated NB data files", length(dat.files)))
flush.console()

# ----- run one muscat iteration based on data file name -----
run_one_muscat = function(dat.name, verbose = TRUE, parallel = TRUE, ddf_choice = "Satterthwaite") {
  seed = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.name)))
  if (is.na(seed)) stop(sprintf("could not parse seed from %s", dat.name))
  
  message(sprintf("seed %02d %s data loading for muscat...", seed, analysis_type))
  flush.console()
  
  sim_obj = readRDS(dat.name)
  sim_obj$metadata$cell_type = "ct1"
  gene_info = sim_obj$gene_info
  
  sce = SingleCellExperiment(
    assays = list(counts = t(sim_obj$count_matrix)),
    colData = sim_obj$metadata
  )
  
  sce = muscat::prepSCE(
    sce,
    kid = "cell_type",
    gid = "group",
    sid = "donor_id",
    drop = TRUE
  )
  
  rm(sim_obj)
  gc()
  
  if (parallel == TRUE) {
    bp = BiocParallel::SnowParam(
      workers = n_cores,
      type = "SOCK",
      progressbar = verbose
    )
    
    mm = muscat::mmDS(
      sce,
      method = "dream",
      n_cells = 10,
      n_samples = 2,
      min_count = 1,
      min_cells = 20,
      ddf = ddf_choice,
      BPPARAM = bp
    )
  } else {
    mm = muscat::mmDS(
      sce,
      method = "dream",
      n_cells = 10,
      n_samples = 2,
      min_count = 1,
      min_cells = 20,
      ddf = ddf_choice
    )
  }
  
  result_pb = mm[["ct1"]]
  
  labels = gene_info$de_status[match(result_pb$gene, gene_info$gene)]
  
  metrics = get_FDR_TPR_metrics(
    pvals = result_pb$p_adj.loc,
    labels = labels
  )
  
  out = c(
    list(seed = seed),
    as.list(metrics)
  )
  
  return(out)
}

# ----- output directories -----
# referenced above
dir.create(muscat_dir, recursive = TRUE, showWarnings = FALSE)

# ----- run muscat sequentially for selected/missing seeds -----
use_parallel = ddf_choice != "Kenward-Roger"

if (!use_parallel) {
  message("Kenward-Roger selected; running without SnowParam parallelization.")
}

message(sprintf("checking muscat outputs for seeds 1 to %02d", n_sim))
flush.console()

muscat_paths = file.path(
  muscat_dir,
  muscat_seed_file(dat.seeds)
)

missing = !file.exists(muscat_paths)

if (any(missing)) {
  message("running missing muscat seeds: ", paste(dat.seeds[missing], collapse = ", "))
} else {
  message("all requested muscat seed outputs already exist")
}

for (i in which(missing)) {
  dat.name = dat.files[i]
  seed = dat.seeds[i]
  path = muscat_paths[i]
  
  err_path = file.path(
    muscat_dir,
    muscat_error_file(seed)
  )
  
  out = tryCatch({
    out = run_one_muscat(
      dat.name = dat.name,
      ddf_choice = ddf_choice,
      parallel = use_parallel
    )
    
    tmp_path = paste0(path, ".tmp")
    saveRDS(out, tmp_path)
    file.rename(tmp_path, path)
    
    if (file.exists(err_path)) file.remove(err_path)
    
    out
  }, error = function(e) {
    msg = sprintf("seed %02d muscat failed: %s", seed, conditionMessage(e))
    message(msg)
    writeLines(msg, err_path)
    NULL
  })
  
  rm(out)
  gc()
}

# ----- combine completed seeds for muscat results -----

results_muscat_list = lapply(seq_along(muscat_paths), function(i) {
  seed = dat.seeds[i]
  path = muscat_paths[i]
  
  if (!file.exists(path)) {
    message(sprintf("seed %02d muscat output missing, skipping in combined results", seed))
    return(NULL)
  }
  
  out = tryCatch(readRDS(path), error = function(e) {
    message(sprintf("seed %02d muscat output unreadable, skipping in combined results", seed))
    NULL
  })
  
  out
})

names(results_muscat_list) = sprintf("seed_%02d", dat.seeds)

failed = vapply(results_muscat_list, is.null, logical(1))

if (any(failed)) {
  message("not included in combined results: ", paste(dat.seeds[failed], collapse = ", "))
}

results_muscat_ok = results_muscat_list[!failed]

results_muscat = do.call(
  rbind,
  lapply(results_muscat_ok, function(x) as.data.frame(as.list(x)))
)

rownames(results_muscat) = NULL
results_muscat = cbind(method = "muscat", results_muscat)

saveRDS(
  results_muscat,
  file.path(muscat_dir, muscat_fname)
)

message("done")


