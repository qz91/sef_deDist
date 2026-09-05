rm(list = ls())
library(pbapply)
library(pbmcapply)
library(BiocParallel)
library(moments)
library(future)
library(future.apply)

source("lpc_simulation_main_5_8_2026_final.R")
source("revisions_p_value_combination_sef_final.R")

`%||%` = function(x, y) if (is.null(x)) y else x

# simulate Gaussian hierarchical for comparison to direct-moment test
generate_gauss_gene_data = function(
    n1, n2,
    lower = 100, upper = 200,
    k1 = NULL, k2 = NULL,
    n_nonDE_unimodal = 500,
    n_DE = 50,
    nonDE_params = list(mu = 0, sigma_sq = 0.2, sigmaC_sq = 0.5),
    DE_G1_params = list(mu = 0, sigma_sq = 0.2, sigmaC_sq = 0.5),
    DE_G2_params = list(mu = 0, sigma_sq = 0.9, sigmaC_sq = 0.5),
    seed = NULL,
    unimodal_fun = generate_unimodal) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  stopifnot(length(lower) == 1, length(upper) == 1, lower > 0, upper >= lower)
  stopifnot(n_nonDE_unimodal >= 0, n_DE >= 0)
  
  if (is.null(k1)) {
    k1 = sample(lower:upper, n1, replace = TRUE)
  }
  if (is.null(k2)) {
    k2 = sample(lower:upper, n2, replace = TRUE)
  }
  
  stopifnot(length(k1) == n1, length(k2) == n2)
  stopifnot(all(k1 > 0), all(k2 > 0), all(k1 %% 1 == 0), all(k2 %% 1 == 0))
  
  donor_ids_1 = sprintf("G1_D%03d", seq_len(n1))
  donor_ids_2 = sprintf("G2_D%03d", seq_len(n2))
  total_cells = sum(k1) + sum(k2)
  
  metadata = data.frame(
    cell_id = sprintf("Cell_%07d", seq_len(total_cells)),
    donor_id = c(rep(donor_ids_1, times = k1), rep(donor_ids_2, times = k2)),
    donor_index_within_group = c(rep(seq_len(n1), times = k1), rep(seq_len(n2), times = k2)),
    group = c(rep("group1", sum(k1)), rep("group2", sum(k2))),
    group_indicator = c(rep(1L, sum(k1)), rep(2L, sum(k2))),
    stringsAsFactors = FALSE
  )
  
  expand_param = function(x, n_genes, label) {
    if (n_genes == 0) {
      return(numeric(0))
    }
    if (length(x) == 1L) {
      return(rep(as.numeric(x), n_genes))
    }
    if (length(x) == n_genes) {
      return(as.numeric(x))
    }
    stop(label, " must have length 1 or length ", n_genes, ".")
  }
  
  prepare_params = function(params, n_genes, label) {
    # validate and expand the simulation parameters to one value per gene
    required = c("mu", "sigma_sq", "sigmaC_sq")
    missing = setdiff(required, names(params))
    if (length(missing) > 0) {
      stop(label, " is missing: ", paste(missing, collapse = ", "))
    }
    
    out = list(
      mu = expand_param(params$mu, n_genes, paste0(label, "$mu")),
      sigma_sq = expand_param(params$sigma_sq, n_genes, paste0(label, "$sigma_sq")),
      sigmaC_sq = expand_param(params$sigmaC_sq, n_genes, paste0(label, "$sigmaC_sq"))
    )
    
    if (any(!is.finite(out$mu)) || #check valiid params for mean and variances
        any(!is.finite(out$sigma_sq)) ||
        any(!is.finite(out$sigmaC_sq))) {
      stop(label, " contains non-finite parameter values.")
    }
    if (any(out$sigma_sq < 0) || any(out$sigmaC_sq < 0)) { 
      stop(label, " variance parameters must be nonnegative.")
    }
    
    return(out) # should return list of validi parameters
  }
  
  simulate_block = function(n_genes, g1_params, g2_params, prefix) {
    if (n_genes <= 0) {
      return(matrix(numeric(0), nrow = total_cells, ncol = 0))
    }
    
    out = vapply(
      seq_len(n_genes),
      function(g) {
        y1 = unimodal_fun(n = n1, k = k1, mu = g1_params$mu[g], sigma_sq = g1_params$sigma_sq[g], sigmaC_sq = g1_params$sigmaC_sq[g])
        
        y2 = unimodal_fun(n = n2, k = k2, mu = g2_params$mu[g], sigma_sq = g2_params$sigma_sq[g], sigmaC_sq = g2_params$sigmaC_sq[g])
        
        c(as.numeric(unlist(y1, use.names = FALSE)), as.numeric(unlist(y2, use.names = FALSE))) # flatten subject-level draws
      },
      numeric(total_cells)
    )
    
    colnames(out) = sprintf("%s_%03d", prefix, seq_len(n_genes))
    out
  }
  
  nonDE_params = prepare_params(nonDE_params, n_nonDE_unimodal, "nonDE_params")
  DE_G1_params = prepare_params(DE_G1_params, n_DE, "DE_G1_params")
  DE_G2_params = prepare_params(DE_G2_params, n_DE, "DE_G2_params")
  
  nonDE_unimodal_mat = simulate_block(n_genes = n_nonDE_unimodal, g1_params = nonDE_params, g2_params = nonDE_params, prefix = "nonDE_unimodal")
  
  DE_mat = simulate_block(n_genes = n_DE, g1_params = DE_G1_params, g2_params = DE_G2_params, prefix = "DE")
  
  count_matrix = cbind(nonDE_unimodal_mat, DE_mat)
  rownames(count_matrix) = metadata$cell_id
  
  gene_ids = colnames(count_matrix)
  
  gene_info = data.frame(
    gene = gene_ids,
    modality = rep("unimodal", length(gene_ids)),
    de_status = c(
      rep("nonDE", n_nonDE_unimodal),
      rep("DE", n_DE)
    ),
    gene_class = c(
      rep("nonDE_unimodal", n_nonDE_unimodal),
      rep("DE_unimodal_vs_unimodal", n_DE)
    ),
    mu_g1 = c(nonDE_params$mu, DE_G1_params$mu),
    mu_g2 = c(nonDE_params$mu, DE_G2_params$mu),
    sigma_sq_g1 = c(nonDE_params$sigma_sq, DE_G1_params$sigma_sq),
    sigma_sq_g2 = c(nonDE_params$sigma_sq, DE_G2_params$sigma_sq),
    sigmaC_sq_g1 = c(nonDE_params$sigmaC_sq, DE_G1_params$sigmaC_sq),
    sigmaC_sq_g2 = c(nonDE_params$sigmaC_sq, DE_G2_params$sigmaC_sq),
    stringsAsFactors = FALSE
  )
  
  return(list(
    count_matrix = count_matrix,
    metadata = metadata,
    gene_info = gene_info,
    k1 = k1,
    k2 = k2,
    settings = list(
      n1 = n1,
      n2 = n2,
      lower = lower,
      upper = upper,
      k1 = k1,
      k2 = k2,
      generator = "generate_unimodal",
      model = "gaussian_hierarchical"
    )
  ))
}

ask_int = function(prompt, default) {
  ans = trimws(readline(sprintf("%s [%d]: ", prompt, default)))
  if (ans == "") return(default)
  
  ans = suppressWarnings(as.integer(ans))
  if (is.na(ans) || ans < 1) {
    stop(sprintf("%s must be a positive integer.", prompt))
  }
  
  return(ans)
}

format_ridge_tag = function(ridge) {
  if (length(ridge) != 1L || !is.numeric(ridge) || !is.finite(ridge) || ridge <= 0) {
    stop("ridge must be a positive finite numeric value.")
  }
  
  scientific = format(ridge, scientific = TRUE, trim = TRUE, digits = 15)
  pieces = strsplit(scientific, "e", fixed = TRUE)[[1]]
  mantissa = chartr(".", "p", pieces[1])
  exponent = as.integer(pieces[2])
  exponent_tag = if (exponent < 0L) {
    paste0("en", abs(exponent))
  } else {
    paste0("e", exponent)
  }
  
  return(paste0("ridge_", mantissa, exponent_tag))
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

format_p_specific_label = function(p_specific) {
  if (is.null(p_specific)) {
    return(NA_character_)
  }
  return(paste(as.integer(p_specific), collapse = "_"))
}

build_total_configs = function(p_values) {
  lapply(p_values, function(p) {
    list(
      label = sprintf("p_total_%d", p),
      p_total = p,
      p_specific = NULL,
      p_specific_label = "1:p"
    )
  })
}

prepare_sim_obj_for_hotelling = function(sim_obj) {
  counts = sim_obj$count_matrix
  meta = sim_obj$metadata
  gene_info = sim_obj$gene_info
  
  if (!all(c("cell_id", "donor_id", "group") %in% names(meta))) {
    stop("metadata must contain cell_id, donor_id, and group.")
  }
  if (anyDuplicated(meta$cell_id)) {
    stop("metadata$cell_id must be unique.")
  }
  
  if (all(meta$cell_id %in% rownames(counts))) {
    counts_use = counts[meta$cell_id, , drop = FALSE]
  } else if (all(meta$cell_id %in% colnames(counts))) {
    counts_use = t(counts[, meta$cell_id, drop = FALSE])
  } else {
    stop("cell IDs in metadata do not match counts rownames or colnames.")
  }
  
  genes = as.character(gene_info$gene)
  keep_gene = genes %in% colnames(counts_use)
  if (!any(keep_gene)) {
    stop("No overlap between gene_info$gene and count_matrix columns.")
  }
  
  gene_info = gene_info[keep_gene, , drop = FALSE]
  genes = as.character(gene_info$gene)
  counts_use = counts_use[, genes, drop = FALSE]
  
  sim_obj$count_matrix = counts_use
  sim_obj$metadata = meta[match(rownames(counts_use), meta$cell_id), , drop = FALSE]
  sim_obj$gene_info = gene_info
  return(sim_obj)
}

make_test_metrics = function(test_obj, alpha = 0.05) {
  test_obj$failed = is.na(test_obj$pVec)
  
  p_for_bh = test_obj$pVec
  p_for_bh[test_obj$failed] = 1
  
  test_obj$p_adj = p.adjust(p_for_bh, method = "BH")
  metrics = get_FDR_TPR_metrics(test_obj$p_adj, test_obj$de_status, alpha = alpha)
  metrics$FDP = metrics$FDR
  metrics$n_failed = sum(test_obj$failed)
  
  return(list(test = test_obj, metrics = metrics))
}

extract_metrics_df = function(out) {
  metrics_df = do.call(rbind, lapply(names(out$results), function(config_name) {
    cfg_result = out$results[[config_name]]
    metrics = data.frame(cfg_result$metrics)
    p_specific_label = cfg_result$p_specific_label
    if (is.null(p_specific_label) || is.na(p_specific_label) || !nzchar(p_specific_label)) {
      p_specific_label = if (is.null(cfg_result$p_specific)) {
        "1:p"
      } else {
        format_p_specific_label(cfg_result$p_specific)
      }
    }
    
    return(cbind(
      data.frame(
        case_id = out$case_id,
        n_samplesize = out$n_samplesize %||% NA_integer_,
        n1 = out$n1 %||% NA_integer_,
        n2 = out$n2 %||% NA_integer_,
        sample_size_tag = out$sample_size_tag %||% NA_character_,
        sigma_sq_DE_G2 = out$sigma_sq_DE_G2,
        sigmaC_sq_DE_G2 = out$sigmaC_sq_DE_G2,
        seed = out$seed,
        method = cfg_result$method,
        config = cfg_result$config_label %||% config_name,
        result_key = config_name,
        test_family = cfg_result$test_family %||% NA_character_,
        p_combination_method = cfg_result$p_combination_method %||% NA_character_,
        moment_set = cfg_result$moment_set %||% p_specific_label,
        data_type = out$data_type %||% NA_character_,
        sef_carrier_type = out$sef_carrier_type %||% NA_character_,
        p = cfg_result$p_total,
        p_total = cfg_result$p_total,
        p_specific = p_specific_label,
        alpha = out$alpha,
        stringsAsFactors = FALSE
      ),
      metrics
    ))
  }))
  return(metrics_df)
}

ask_case_selection = function(prompt, n_cases, default = "all") {
  repeat {
    ans = trimws(readline(sprintf("%s [%s]: ", prompt, default)))
    if (ans == "") ans = default
    
    ans_lower = tolower(ans)
    if (ans_lower == "all") {
      return(seq_len(n_cases))
    }
    
    pieces = unlist(strsplit(ans, "[,[:space:]]+"))
    pieces = pieces[pieces != ""]
    selected = suppressWarnings(as.integer(pieces))
    
    if (
      length(selected) > 0L &&
      !any(is.na(selected)) &&
      all(selected >= 1L) &&
      all(selected <= n_cases)
    ) {
      return(unique(selected))
    }
    
    message(sprintf(
      "Please enter all, or case numbers between 1 and %d separated by spaces or commas.",
      n_cases
    ))
  }
}

safe_num = function(x) {
  return(gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE)))
}

# ----- downstream loading/setting up for simulation -----

default_n_samplesize = 50
n_samplesize = ask_int(
  "Group-specific sample size, used for both n1 and n2",
  default_n_samplesize
)
n1 = n_samplesize
n2 = n_samplesize
sample_size_tag = sprintf("n_%d", n_samplesize)
lower = 100
upper = 200
k1 = NULL
k2 = NULL
n_nonDE_unimodal = 500
n_DE = 50
DE_params1 = list(mu = 0, sigma_sq = 0.20, sigmaC_sq = 1)

make_param_case = function(sigma_sq_DE_G2, sigmaC_sq_shared) {
  case_id = sprintf(
    "DE2_sigmasq_%s_shared_sigmaC_%s",
    safe_num(sigma_sq_DE_G2),
    safe_num(sigmaC_sq_shared)
  )
  
  DE_G1_params = list(
    mu = DE_params1$mu,
    sigma_sq = DE_params1$sigma_sq,
    sigmaC_sq = sigmaC_sq_shared
  )
  
  list(
    case_id = case_id,
    sigma_sq_DE_G2 = sigma_sq_DE_G2,
    sigmaC_sq_DE_G2 = sigmaC_sq_shared,
    sigmaC_sq_shared = sigmaC_sq_shared,
    nonDE_params = list(
      mu = DE_G1_params$mu,
      sigma_sq = DE_G1_params$sigma_sq,
      sigmaC_sq = sigmaC_sq_shared
    ),
    DE_params1 = DE_G1_params,
    DE_params2 = list(
      mu = DE_G1_params$mu,
      sigma_sq = sigma_sq_DE_G2,
      sigmaC_sq = sigmaC_sq_shared
    )
  )
}

sigma_sq_DE_G2_values  = c(0.7, 0.9)
sigmaC_sq_DE_G2_values = c(0.7, 0.8)

param_grid = expand.grid(
  sigma_sq_DE_G2 = sigma_sq_DE_G2_values,
  sigmaC_sq_DE_G2 = sigmaC_sq_DE_G2_values
)

param_cases = lapply(seq_len(nrow(param_grid)), function(i) {
  make_param_case(
    sigma_sq_DE_G2 = param_grid$sigma_sq_DE_G2[i],
    sigmaC_sq_shared = param_grid$sigmaC_sq_DE_G2[i]
  )
})

n_sims = ask_int("Number of simulations", 100)
n_cores = ask_int("Number of cores", 20)
data_type = ask_choice("Data type", c("gaussian"), "gaussian")
sef_carrier_type = "continuous"
p_max = ask_int("Maximum moment to test", 4)
if (p_max < 2L) {
  stop("Maximum moment must be at least 2 for p-value combination.")
}
p_values = 2:p_max
pcombo_methods = c("bonferroni")
ridge = 1e-8
ridge_tag = format_ridge_tag(ridge)
config_mode = "momentwise_pvalue_combination"
weight_method = "unit"
weight_col = "adj_weight"
analysis_dir = sprintf("%s_pcombo_1_to_%d", data_type, p_max)

message("Available Gaussian special-case parameter settings:")
for (i in seq_along(param_cases)) {
  message(sprintf(
    "%02d: %s, sigma_sq_DE_G2 = %.2f, sigmaC_sq_DE_G2 = %.2f",
    i,
    param_cases[[i]]$case_id,
    param_cases[[i]]$sigma_sq_DE_G2,
    param_cases[[i]]$sigmaC_sq_DE_G2
  ))
}

selected_cases = ask_case_selection(
  "Parameter set numbers, separated by spaces or commas, or all",
  length(param_cases),
  "all"
)
param_cases = param_cases[selected_cases]

out.root = "/pval_combo"
sample_size_root = file.path(out.root, sample_size_tag)
dir.create(sample_size_root, recursive = TRUE, showWarnings = FALSE)

make_case_dir = function(case) {
  file.path(sample_size_root, case$case_id)
}

make_combo_dir = function(case) {
  file.path(make_case_dir(case), analysis_dir)
}

make_seed_file = function(case, seed) {
  file.path(
    make_combo_dir(case),
    sprintf(
      "pval_combo_%s_%s_%s_%s_seed_%03d.rds",
      sample_size_tag,
      case$case_id,
      analysis_dir,
      ridge_tag,
      seed
    )
  )
}

make_err_file = function(case, seed) {
  file.path(
    make_combo_dir(case),
    sprintf(
      "pval_combo_%s_%s_%s_%s_seed_%03d_ERROR.txt",
      sample_size_tag,
      case$case_id,
      analysis_dir,
      ridge_tag,
      seed
    )
  )
}

make_metrics_file = function(case) {
  file.path(
    make_combo_dir(case),
    sprintf(
      "pval_combo_%s_%s_%s_%s_metrics.RDS",
      sample_size_tag,
      case$case_id,
      analysis_dir,
      ridge_tag
    )
  )
}

make_metrics_list_file = function(case) {
  file.path(
    make_combo_dir(case),
    sprintf(
      "pval_combo_%s_%s_%s_%s_metrics_list.RDS",
      sample_size_tag,
      case$case_id,
      analysis_dir,
      ridge_tag
    )
  )
}

generate_one_seed = function(seed, case) {
  generate_gauss_gene_data(
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    k1 = k1,
    k2 = k2,
    n_nonDE_unimodal = n_nonDE_unimodal,
    n_DE = n_DE,
    nonDE_params = case$nonDE_params,
    DE_G1_params = case$DE_params1,
    DE_G2_params = case$DE_params2,
    seed = seed
  )
}

analysis_configs = build_total_configs(p_values)
alpha = 0.05

message("beginning Gaussian special-case moment-wise p-value combination tests")
message("Configs: ", paste(vapply(analysis_configs, `[[`, character(1), "label"), collapse = ", "))
message("p-value combination methods: ", paste(pcombo_methods, collapse = ", "))
message("data_type: ", data_type)
message("SEF carrier type: ", sef_carrier_type)
message("weight_method: ", weight_method)
message("group-specific sample size: ", n_samplesize, " (", sample_size_tag, ")")
message("ridge: ", ridge, " (", ridge_tag, ")")
message("output root: ", sample_size_root)
flush.console()

make_result_entry = function(method, test, p_total, p_specific = NULL,
                             p_specific_label = NULL, test_family,
                             p_combination_method = NA_character_,
                             moment_set = NULL, config_label = NULL,
                             details = NULL) {
  test_metrics = make_test_metrics(test, alpha = alpha)
  list(
    method = method,
    test = test_metrics$test,
    details = details,
    metrics = test_metrics$metrics,
    p_total = p_total,
    p_specific = p_specific,
    p_specific_label = p_specific_label %||% if (is.null(p_specific)) {
      "1:p"
    } else {
      format_p_specific_label(p_specific)
    },
    test_family = test_family,
    p_combination_method = p_combination_method,
    moment_set = moment_set %||% if (is.null(p_specific)) {
      paste0("1:", p_total)
    } else {
      format_p_specific_label(p_specific)
    },
    config_label = config_label
  )
}

subset_momentwise_test = function(momentwise_df, moment) {
  p_col = moment_pvalue_colnames(moment)
  pVec = setNames(as.numeric(momentwise_df[[p_col]]), as.character(momentwise_df$gene))
  gene_cols = !grepl("^X_[0-9]+$", names(momentwise_df))
  make_pvec_df(momentwise_df[, gene_cols, drop = FALSE], pVec)
}

run_one_combo = function(seed, case) {
  message(sprintf("case %s seed %03d data generating...", case$case_id, seed))
  sim_obj = prepare_sim_obj_for_hotelling(generate_one_seed(seed, case))
  flush.console()

  all_moments = seq_len(p_max)

  message(sprintf("case %s seed %03d running SEF moment-wise tests", case$case_id, seed))
  sef_momentwise = run_sef_test_covfix_momentwise(
    sim_obj = sim_obj,
    p_total = p_max,
    moment_idx = all_moments,
    carrier_type = sef_carrier_type,
    scaled_basis = TRUE,
    verbose = FALSE,
    bw = 0.5,
    ridge = ridge
  )

  message(sprintf("case %s seed %03d running Hotelling component-wise tests", case$case_id, seed))
  hotelling_componentwise = run_hotelling_mom_test_covfix_componentwise(
    sim_obj = sim_obj,
    p = p_max,
    p_total = p_max,
    moment_idx = all_moments,
    verbose = FALSE,
    ridge = ridge,
    parallel = FALSE,
    n_cores = 1L
  )

  run_hotelling_omnibus_for_cfg = function(cfg) {
    test = run_hotelling_mom_test_covfix_combination(
      sim_obj = sim_obj,
      p = cfg$p_total,
      p_total = cfg$p_total,
      p_specific = NULL,
      verbose = FALSE,
      ridge = ridge,
      parallel = FALSE,
      n_cores = 1L
    )
    details = attr(test, "details")
    attr(test, "details") = NULL

    make_result_entry(
      method = "hotelling",
      test = test,
      details = details,
      p_total = cfg$p_total,
      p_specific = NULL,
      p_specific_label = "1:p",
      test_family = "omnibus",
      p_combination_method = "omnibus",
      moment_set = paste0("1:", cfg$p_total),
      config_label = cfg$label
    )
  }

  run_sef_omnibus_for_cfg = function(cfg) {
    test = run_sef_test_covfix_combination(
      sim_obj = sim_obj,
      p = cfg$p_total,
      p_total = cfg$p_total,
      p_specific = NULL,
      carrier_type = sef_carrier_type,
      scaled_basis = TRUE,
      verbose = FALSE,
      bw = 0.5,
      ridge = ridge
    )
    details = attr(test, "details")
    attr(test, "details") = NULL

    make_result_entry(
      method = "sef",
      test = test,
      details = details,
      p_total = cfg$p_total,
      p_specific = NULL,
      p_specific_label = "1:p",
      test_family = "omnibus",
      p_combination_method = "omnibus",
      moment_set = paste0("1:", cfg$p_total),
      config_label = cfg$label
    )
  }

  results_by_config = list()

  for (moment in all_moments) {
    moment_label = as.character(moment)
    sef_single = subset_momentwise_test(sef_momentwise, moment)
    hotelling_single = subset_momentwise_test(hotelling_componentwise, moment)

    results_by_config[[sprintf("sef_momentwise_%s", moment_label)]] = make_result_entry(
      method = "sef",
      test = sef_single,
      p_total = p_max,
      p_specific = moment,
      p_specific_label = moment_label,
      test_family = "momentwise",
      p_combination_method = "none",
      moment_set = moment_label,
      config_label = paste0("p_specific_", moment_label)
    )

    results_by_config[[sprintf("hotelling_momentwise_%s", moment_label)]] = make_result_entry(
      method = "hotelling",
      test = hotelling_single,
      p_total = p_max,
      p_specific = moment,
      p_specific_label = moment_label,
      test_family = "momentwise",
      p_combination_method = "none",
      moment_set = moment_label,
      config_label = paste0("p_specific_", moment_label)
    )
  }

  for (cfg in analysis_configs) {
    moment_idx = seq_len(cfg$p_total)
    moment_set_label = paste0("1:", cfg$p_total)
    pcombo_config_label = paste0("pcombo_1_to_", cfg$p_total)
    p_specific_label = format_p_specific_label(moment_idx)

    results_by_config[[paste("sef_omnibus", cfg$label, sep = "_")]] =
      run_sef_omnibus_for_cfg(cfg)
    results_by_config[[paste("hotelling_omnibus", cfg$label, sep = "_")]] =
      run_hotelling_omnibus_for_cfg(cfg)

    sef_pmat = as.matrix(sef_momentwise[, moment_pvalue_colnames(moment_idx), drop = FALSE])
    hotelling_pmat = as.matrix(hotelling_componentwise[, moment_pvalue_colnames(moment_idx), drop = FALSE])

    sef_combo = combine_momentwise_pmat(
      gene_info = sim_obj$gene_info,
      pmat = sef_pmat,
      moment_idx = moment_idx,
      methods = pcombo_methods,
      na_rm = FALSE
    )
    hotelling_combo = combine_momentwise_pmat(
      gene_info = sim_obj$gene_info,
      pmat = hotelling_pmat,
      moment_idx = moment_idx,
      methods = pcombo_methods,
      na_rm = FALSE
    )

    for (combo_method in pcombo_methods) {
      results_by_config[[sprintf("sef_pcombo_%s_%s", combo_method, cfg$label)]] =
        make_result_entry(
          method = "sef",
          test = sef_combo[[combo_method]],
          p_total = cfg$p_total,
          p_specific = moment_idx,
          p_specific_label = p_specific_label,
          test_family = "pvalue_combination",
          p_combination_method = combo_method,
          moment_set = moment_set_label,
          config_label = pcombo_config_label
        )

      results_by_config[[sprintf("hotelling_pcombo_%s_%s", combo_method, cfg$label)]] =
        make_result_entry(
          method = "hotelling",
          test = hotelling_combo[[combo_method]],
          p_total = cfg$p_total,
          p_specific = moment_idx,
          p_specific_label = p_specific_label,
          test_family = "pvalue_combination",
          p_combination_method = combo_method,
          moment_set = moment_set_label,
          config_label = pcombo_config_label
        )
    }

    message(sprintf("case %s seed %03d finished %s", case$case_id, seed, cfg$label))
    flush.console()
  }
  
  list(
    seed = seed,
    case_id = case$case_id,
    n_samplesize = n_samplesize,
    n1 = n1,
    n2 = n2,
    sample_size_tag = sample_size_tag,
    sigma_sq_DE_G2 = case$sigma_sq_DE_G2,
    sigmaC_sq_DE_G2 = case$sigmaC_sq_DE_G2,
    sigmaC_sq_shared = case$sigmaC_sq_shared,
    results = results_by_config,
    weight_method = weight_method,
    weight_col = weight_col,
    ridge = ridge,
    config_mode = config_mode,
    pcombo_methods = pcombo_methods,
    data_type = data_type,
    sef_carrier_type = sef_carrier_type,
    p_max = p_max,
    alpha = alpha
  )
}

summarize_metrics = function(metrics_df) {
  if (is.null(metrics_df) || !nrow(metrics_df)) {
    return(data.frame())
  }
  
  value_cols = intersect(c("FDR", "TPR", "FDP", "n_failed"), names(metrics_df))
  by_cols = intersect(
    c(
      "case_id",
      "n_samplesize",
      "n1",
      "n2",
      "sample_size_tag",
      "sigma_sq_DE_G2",
      "sigmaC_sq_DE_G2",
      "method",
      "config",
      "test_family",
      "p_combination_method",
      "moment_set",
      "data_type",
      "sef_carrier_type",
      "p",
      "p_total",
      "p_specific",
      "alpha"
    ),
    names(metrics_df)
  )
  
  aggregate(
    metrics_df[, value_cols, drop = FALSE],
    metrics_df[, by_cols, drop = FALSE],
    FUN = function(x) mean(x, na.rm = TRUE)
  )
}

all_case_metrics = list()
case_seeds = seq_len(n_sims)

for (case in param_cases) {
  combo_dir = make_combo_dir(case)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Running Gaussian special-case setting: ", case$case_id)
  message("output directory: ", combo_dir)
  flush.console()
  
  seed_paths = vapply(case_seeds, function(seed) make_seed_file(case, seed), character(1))
  missing = !file.exists(seed_paths)
  
  if (any(missing)) {
    message("running missing seeds: ", paste(case_seeds[missing], collapse = ", "))
  } else {
    message("all requested seed outputs already exist")
  }
  
  run_missing_seed = function(seed) {
    path = make_seed_file(case, seed)
    err_path = make_err_file(case, seed)
    
    ok = tryCatch({
      out = run_one_combo(seed, case)
      
      tmp_path = paste0(path, ".tmp")
      saveRDS(out, tmp_path)
      file.rename(tmp_path, path)
      
      if (file.exists(err_path)) file.remove(err_path)
      
      rm(out)
      gc()
      TRUE
    }, error = function(e) {
      msg = sprintf("case %s seed %03d failed: %s", case$case_id, seed, conditionMessage(e))
      message(msg)
      writeLines(msg, err_path)
      FALSE
    })
    
    ok
  }
  
  missing_seeds = case_seeds[missing]
  if (length(missing_seeds) > 0L) {
    run_status = pbmcapply::pbmclapply(
      missing_seeds,
      run_missing_seed,
      mc.cores = min(n_cores, length(missing_seeds))
    )
    names(run_status) = sprintf("seed_%03d", missing_seeds)
    failed_run = !unlist(run_status, use.names = FALSE)
    if (any(failed_run)) {
      message(
        "Seeds failed during this run for case ",
        case$case_id,
        ": ",
        paste(missing_seeds[failed_run], collapse = ", ")
      )
    }
    gc()
  }
  
  results_list = lapply(seq_along(seed_paths), function(i) {
    seed = case_seeds[i]
    path = seed_paths[i]
    
    if (!file.exists(path)) {
      message(sprintf("case %s seed %03d output missing, skipping in combined results", case$case_id, seed))
      return(NULL)
    }
    
    tryCatch(readRDS(path), error = function(e) {
      message(sprintf("case %s seed %03d output unreadable, skipping in combined results", case$case_id, seed))
      NULL
    })
  })
  
  names(results_list) = sprintf("seed_%03d", case_seeds)
  failed = vapply(results_list, is.null, logical(1))
  
  if (any(failed)) {
    message("not included for case ", case$case_id, ": ", paste(case_seeds[failed], collapse = ", "))
  }
  
  results_ok = results_list[!failed]
  
  case_metrics = if (length(results_ok) > 0L) {
    do.call(rbind, lapply(results_ok, extract_metrics_df))
  } else {
    data.frame()
  }
  rownames(case_metrics) = NULL
  
  saveRDS(case_metrics, make_metrics_list_file(case))
  saveRDS(case_metrics, make_metrics_file(case))
  
  all_case_metrics[[case$case_id]] = case_metrics
  
  message("saved metrics to: ", make_metrics_file(case))
}

all_metrics = do.call(rbind, all_case_metrics)
rownames(all_metrics) = NULL

combined_metrics_file = file.path(
  sample_size_root,
  sprintf(
    "pval_combo_%s_%s_%s_all_cases_metrics.RDS",
    sample_size_tag,
    analysis_dir,
    ridge_tag
  )
)
combined_summary_file = file.path(
  sample_size_root,
  sprintf(
    "pval_combo_%s_%s_%s_all_cases_summary.RDS",
    sample_size_tag,
    analysis_dir,
    ridge_tag
  )
)
combined_summary_csv = file.path(
  sample_size_root,
  sprintf(
    "pval_combo_%s_%s_%s_all_cases_summary.csv",
    sample_size_tag,
    analysis_dir,
    ridge_tag
  )
)

saveRDS(all_metrics, combined_metrics_file)

summary_metrics = summarize_metrics(all_metrics)
saveRDS(summary_metrics, combined_summary_file)
write.csv(summary_metrics, combined_summary_csv, row.names = FALSE)

message("saved combined metrics to: ", combined_metrics_file)
message("saved combined summary to: ", combined_summary_file)
