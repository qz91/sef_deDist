# lpc_NB_wilcox_cell_generalized.R
rm(list = ls())
library(pbmcapply)
library(pbapply)
library(presto)
library(Seurat)
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
  
  wilcox_dir = file.path(
    "/home/zaqian/sim_par",
    sprintf("%s_effect", effect_tag),
    sprintf("%s_wilcox_cell_%s", effect_tag, shift_tag)
  )
  
  wilcox_fname = sprintf(
    "NB_wilcox_cell_%s_effect_%s_metrics.RDS",
    effect_tag,
    shift_tag
  )
  
  wilcox_seed_file = function(seed) {
    sprintf("NB_wilcox_cell_%s_effect_%s_seed_%02d.rds", effect_tag, shift_tag, seed)
  }
  
} else {
  dat.dir = file.path(base_dat_dir, "modal_NB")
  dat_pattern = "^dat_NB_modal_seed_[0-9]+\\.rds$"
  
  wilcox_dir = "/home/zaqian/sim_par/NB_modal/wilcox_cell_modal"
  wilcox_fname = "NB_modal_wilcox_cell_metrics.RDS"
  
  wilcox_seed_file = function(seed) {
    sprintf("NB_modal_wilcox_cell_seed_%02d.rds", seed)
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




# ----- run wilcox cell -----
run_one_wilcox = function(dat.name) {
  seed = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.name)))
  if (is.na(seed)) stop(sprintf("could not parse seed from %s", dat.name))
  
  message(sprintf("seed %02d %s data loading for wilcox cell...", seed, analysis_type))
  sim_obj = readRDS(dat.name)
  flush.console()
  
  gene_info = sim_obj$gene_info
  rownames(sim_obj$metadata) = sim_obj$metadata$cell_id
  sim_obj$metadata$cell_id = NULL
  
  test_sObj = Seurat::CreateSeuratObject(
    counts = t(sim_obj$count_matrix),
    meta.data = sim_obj$metadata
  )
  
  rm(sim_obj)
  gc()
  
  Idents(test_sObj) = "group"
  test_sObj = Seurat::NormalizeData(test_sObj)
  
  test_wilcox = Seurat::FindMarkers(
    test_sObj,
    ident.1 = "group1",
    ident.2 = "group2",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = -Inf,
    verbose = FALSE
  )
  
  rm(test_sObj)
  gc()
  
  test_wilcox$p_val_fdr = p.adjust(test_wilcox$p_val, method = "BH")
  rownames(test_wilcox) = gsub("-", "_", rownames(test_wilcox))
  
  labels = gene_info$de_status[match(rownames(test_wilcox), gene_info$gene)]
  
  metrics = get_FDR_TPR_metrics(
    pvals = test_wilcox$p_val_fdr,
    labels = labels
  )
  
  out = c(
    list(seed = seed),
    as.list(metrics)
  )
  
  return(out)
}


dir.create(wilcox_dir, recursive = TRUE, showWarnings = FALSE)

# ----- run in wilcox cell parallel -----
message(sprintf("running %02d preloaded wilcox cell files with %02d cores", n_sim, n_cores))
flush.console()

results_wilcox_list = pbmcapply::pbmclapply(
  dat.files,
  function(dat.name) {
    seed = as.integer(sub(".*_seed_([0-9]+)\\.rds$", "\\1", basename(dat.name)))
    
    path = file.path(
      wilcox_dir,
      wilcox_seed_file(seed)
    )
    
    if (file.exists(path)) {
      message(sprintf("seed %02d wilcox cell already exists, loading...", seed))
      out = readRDS(path)
    } else {
      out = run_one_wilcox(dat.name)
      saveRDS(out, path)
    }
    
    out
  },
  mc.cores = n_cores
)

names(results_wilcox_list) = sprintf("seed_%02d", dat.seeds)

results_wilcox = do.call(
  rbind,
  lapply(results_wilcox_list, function(x) as.data.frame(as.list(x)))
)

rownames(results_wilcox) = NULL
results_wilcox = cbind(method = "wilcox cell-level", results_wilcox)

saveRDS(
  results_wilcox,
  file.path(wilcox_dir, wilcox_fname)
)

