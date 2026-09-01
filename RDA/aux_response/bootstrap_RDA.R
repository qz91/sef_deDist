# Wider pool of genes tested, more donors, criteria dropped to nonzero counts total of at least 5%?
rm(list = ls())  
library(dplyr)
library(ggplot2)
library(Seurat)
source_dir = file.path("/Users/zaqian/Desktop/merge_sim_RDA")
bootstrap_dir = file.path(source_dir,"bootstrap_results")
source(file.path(source_dir, "merged_self_contained_RDA.R"))

#----- load data ----- 
cell_type = "cd8"

## ----- otherwise, load data -----
if(cell_type == "cd4"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd4_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
  
} else if (cell_type == "cd8"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cd8_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
  
} else if (cell_type == "cM"){
  lcf_matrix = readRDS("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/cM_lcf_matrix.RDS")
  meta_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/", cell_type,"_new_metadata.RDS")
  sObj_meta = readRDS(meta_name)
}
cell_counts = sObj_meta %>%
  dplyr::count(donor_id, disease, disease_state, ind_cov, name = "n_cells") %>%
  dplyr::arrange(desc(n_cells))
cell_counts = cell_counts[order(cell_counts$n_cells, decreasing = F), ]
counts_dir = file.path(source_dir, "counts_by_donor")
# dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
# write.csv(cell_counts, file.path(counts_dir, paste0(cell_type,"_cell_counts_by_donor.csv")))
head(cell_counts,10)

dup_donors = names(which(table(cell_counts$ind_cov) > 1))
stopifnot(
  "Duplicate ind_cov values found across donors; please remove flare samples." = length(dup_donors) == 0
)

donor_ids = sObj_meta$donor_id
length(donor_ids)
dim(lcf_matrix)
stopifnot(
  "Column names of library-corrected count matrix do not match rownames of sObj_meta" = 
    all(colnames(lcf_matrix) == rownames(sObj_meta)),
  "Length of donor_ids does not match number of columns in library-corrected count matrix" = 
    length(donor_ids) == ncol(lcf_matrix)
)


# Within each cell-type group, we filter out donors with fewer than 300 cells, 
# and filter out genes with donor-level cell-type-specific nonzero expression rates below 3%. 

lcf_filtered = filter_SLE(lcf_matrix, donor_ids, min_cells = 300, min_prop = 0.03)
genes_of_interest = lcf_filtered$genes
donors_keep = lcf_filtered$donors
# donors_name = paste0("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/new_seurat/comparisons_4_20_2026/", cell_type,"_donors_used.RDS")
# saveRDS(donors_keep, donors_name)
sObj_meta = sObj_meta[sObj_meta$donor_id %in% lcf_filtered$donors, ]
lcf2= lcf_filtered$matrix
bootstrap_dir = file.path(source_dir,"bootstrap_results")
saveRDS(lcf2, file.path(bootstrap_dir,paste0(cell_type,"_init_lcf2.RDS")))
saveRDS(sObj_meta,file.path(bootstrap_dir,paste0(cell_type,"_init_metadata.RDS")))
dim(lcf2); dim(sObj_meta)
rm(donor_ids) # save memory
rm(lcf_matrix)
rm(lcf_filtered)
gc()


# ----- load things -----
lcf2 = readRDS(file.path(bootstrap_dir,paste0(cell_type,"_init_lcf2.RDS")))
sObj_meta = readRDS(file.path(bootstrap_dir,paste0(cell_type,"_init_metadata.RDS")))
genes_of_interest = row.names(lcf2)
# ----- internal checking -----

donor_key = sObj_meta %>% # one row per eligible donor
  transmute(
    source_donor_id = as.character(donor_id),
    condition = if_else(
      disease == "normal",
      "normal",
      "SLE"
    )
  ) %>%
  distinct()

# Confirm that every donor maps to exactly one condition.
stopifnot(nrow(donor_key) == dplyr::n_distinct(donor_key$source_donor_id) )

# Cell positions belonging to each original donor.
donor_cells = split(
  seq_len(nrow(sObj_meta)),
  as.character(sObj_meta$donor_id)
)

# Donor IDs separated by condition.
donors_by_condition = split(
  donor_key$source_donor_id,
  donor_key$condition
)

print(lengths(donors_by_condition))

# ----- bootstrapping internals -----
B = 20
bootstrap_seed = 20260828
make_one_bootstrap = function(bootstrap_id, seed = bootstrap_seed) {
  # An independent, reproducible seed for each bootstrap.
  set.seed(seed + bootstrap_id)
  
  draw_table = bind_rows(
    lapply(names(donors_by_condition), function(group_name) {
      donor_pool = donors_by_condition[[group_name]]
      data.frame(
        condition = group_name,
        source_donor_id = sample(
          donor_pool,
          size = length(donor_pool),
          replace = TRUE
        ),
        draw_in_condition = seq_along(donor_pool),
        stringsAsFactors = FALSE
      )
    })
  )
  
  # A donor selected multiple times must have a different ID
  # for every draw.
  draw_table$bootstrap_donor_id <- sprintf(
    "boot%02d_%s_draw%03d",
    bootstrap_id,
    draw_table$condition,
    draw_table$draw_in_condition
  )
  
  # Retrieve all cells belonging to each selected donor.
  selected_cell_lists <- unname(
    donor_cells[draw_table$source_donor_id]
  )
  
  cells_per_draw <- lengths(selected_cell_lists)
  
  draw_table$n_cells <- cells_per_draw
  
  list(
    bootstrap_id = bootstrap_id,
    draw_table = draw_table,
    cell_index = unlist(
      selected_cell_lists,
      use.names = FALSE
    ),
    cells_per_draw = cells_per_draw
  )
}

# ----- generate B = 20 bootstrap draws -----
bootstrap_draws = lapply(seq_len(B), make_one_bootstrap)
saveRDS(bootstrap_draws,
        file.path(bootstrap_dir,paste0(cell_type,"_donor_bootstrap_draws_B",B,".RDS"))
)

# ----- manual inspections -----
boot1 = bootstrap_draws[[1]] # look at first draw

table(boot1$draw_table$condition) # same donor sample size within each condition as original RDA

sum(duplicated(boot1$draw_table$source_donor_id)) # print number of repeated original donors.

head(boot1$draw_table)

# ----- make bootstrap metadata -----
make_bootstrap_metadata = function(bootstrap_draw, original_metadata) {
  boot_meta = original_metadata[bootstrap_draw$cell_index, , drop = FALSE] # cell_index likely contains repeats via bootstrap
  
  boot_meta$source_donor_id = as.character(boot_meta$donor_id) # get OG donor_ids
  
  synthetic_ids = rep(bootstrap_draw$draw_table$bootstrap_donor_id, bootstrap_draw$cells_per_draw)
  
  boot_meta$donor_id = synthetic_ids # SEF now treats repeated draws as separate bootstrap donors.
  
  if ("ind_cov" %in% colnames(boot_meta)) { # used for checking to make sure no ICC controls
    boot_meta$source_ind_cov = as.character(boot_meta$ind_cov)
    boot_meta$ind_cov = synthetic_ids
  }
  
  
  rownames(boot_meta) = make.unique( # make bootstrap IDs for cells
    paste0( # repeat donors contain repeated OG cell IDs so we construct unique bootstrap cell names
      synthetic_ids,
      "__",
      rownames(boot_meta)
    )
  )
  
  return(boot_meta)
}

boot1_meta = make_bootstrap_metadata( # test, it's pretty fast
  bootstrap_draw = boot1,
  original_metadata = sObj_meta
)

gene <- genes_of_interest[1]

exprVec_boot1 <- lcf2[
  gene,
  boot1$cell_index,
  drop = FALSE
]

fit_boot1 <- run_sef_procedure(
  exprVec = exprVec_boot1,
  sObj_meta = boot1_meta,
  p = 2,
  bandwidth = 0.5,
  plot_flag = FALSE
)

fit_boot1$pval_1

# ------ bootstrapping SEF -----
B = 20
p = 2
bw = 0.5
alpha = 0.05
analysis_seed = 20260828

options(future.globals.maxSize = 5 * 1024^3)  # 3 GiB
plan(multisession, workers = 2)
handlers(global = TRUE)  # enables progress bars

res_dir = file.path(bootstrap_dir, paste0(cell_type, "_res_objects"))
pval_dir = file.path(bootstrap_dir, paste0(cell_type, "_bootstrap_pvals"))

# load existing bootstrap draws in files if not already
if (!exists("bootstrap_draws")) {
  bootstrap_draws = readRDS(file.path(bootstrap_dir, paste0(cell_type,"_donor_bootstrap_draws_B",B,".RDS")))
}

n_tests = length(genes_of_interest)

for (b in seq_len(B)) {
  message("\nStarting bootstrap ",b," of ", B)
  
  bootstrap_draw = bootstrap_draws[[b]]
  
  boot_meta = make_bootstrap_metadata( # generate bootstrapped metadata for b of B
    bootstrap_draw = bootstrap_draw,
    original_metadata = sObj_meta
  )
  
  stopifnot(nrow(boot_meta) == length(bootstrap_draw$cell_index),length(unique(boot_meta$donor_id)) == nrow(bootstrap_draw$draw_table))
  
  bootstrap_tag = sprintf("bootstrap_%03d", b)
  
  res_file = file.path(res_dir, paste0("res_", cell_type, "_", bootstrap_tag, "_p_", p, "_sef.RDS"))
  pval_file = file.path(pval_dir, paste0(cell_type, "_", bootstrap_tag, "_p_", p,"_bonferroni_pvals.csv"))
  sig_file = file.path(pval_dir, paste0(cell_type, "_", bootstrap_tag,"_p_", p, "_bonferroni_significant.csv"))
  
  if (file.exists(res_file) && file.exists(pval_file)) { # makes the loop resumable
    message(bootstrap_tag," already completed; skipping.")
    next
  }
  
  
  if (file.exists(res_file)) { # If the result object already exists but the p-value, file does not, load it without rerunning SEF.
    message("Loading existing result object for ", bootstrap_tag)
    new_sef = readRDS(res_file)
  } else {
    set.seed(analysis_seed + b) # Reproducible RNG for this bootstrap analysis.
    
    new_sef = with_progress({
      p_bar = progressor(steps = length(genes_of_interest))
      
      future_lapply(
        genes_of_interest,
        function(gene) {
          p_bar(gene)
          
          tryCatch({
            # Duplicate only this gene's expression vector,
            # rather than the entire expression matrix.
            exprVec = lcf2[gene, bootstrap_draw$cell_index, drop = FALSE]
            
            run_sef_procedure(
              exprVec = exprVec,
              sObj_meta = boot_meta,
              p = p,
              bandwidth = bw,
              plot_flag = FALSE
            )
            
          }, error = function(e) {
            list(
              gene = gene,
              pval_1 = NA_real_,
              error = conditionMessage(e)
            )
          })
        },
        future.seed = TRUE,
        
        # Use one chunk per worker, avoiding excessive
        # repeated serialization of lcf2.
        future.scheduling = 1
      )
    })
    
    new_sef = setNames(new_sef, genes_of_interest)
    
    # Save immediately so this bootstrap does not need
    # to be rerun if a later replicate fails.
    saveRDS(new_sef, res_file)
    
    message("Saved result object: ", res_file)
  }
  
  # Extract exactly one p-value per tested gene.
  pvals = vapply(
    new_sef,
    function(x) {
      if (!is.list(x) || is.null(x$pval_1) || length(x$pval_1) != 1L) { # check if it didn't run properly
        return(NA_real_)
      }
      
      value = as.numeric(x$pval_1)
      
      if (!is.finite(value)) { #if the pval is NA
        return(NA_real_)
      }
      
      value
    },
    numeric(1)
  )
  
  # Record errors alongside the p-values.
  error_messages = vapply(
    new_sef,
    function(x) {if (is.list(x) && !is.null(x$error)) {
        return(as.character(x$error))
      }
      
      NA_character_
    },
    character(1)
  )
  
  pval_df = data.frame(
    bootstrap_id = b,
    gene = names(pvals),
    p_value = unname(pvals),
    error = unname(error_messages),
    stringsAsFactors = FALSE
  )
  
  # Correct across the complete, fixed gene universe
  # for this bootstrap.
  pval_df$bonferroni_pval = p.adjust(
    pval_df$p_value,
    method = "bonferroni" # optionally , n = n_tests
  )
  
  pval_df$significant_bonferroni =
    !is.na(pval_df$bonferroni_pval) &
    pval_df$bonferroni_pval < 0.05
  
  pval_df = pval_df %>%
    arrange(
      bonferroni_pval,
      p_value
    )
  
  # Save all raw and corrected p-values.
  write.csv(
    pval_df,
    pval_file,
    row.names = FALSE,
    na = ""
  )
  
  # Also save only the significant genes.
  write.csv(subset(pval_df, significant_bonferroni),
    sig_file,
    row.names = FALSE,
    na = ""
  )
  
  message(
    "Completed ",
    bootstrap_tag,
    ": ",
    sum(
      pval_df$significant_bonferroni,
      na.rm = TRUE
    ),
    " Bonferroni-significant genes; ",
    sum(is.na(pval_df$p_value)),
    " failed gene fits."
  )
  
  # Release replicate-specific objects before the
  # next sequential bootstrap.
  rm(
    new_sef,
    pval_df,
    pvals,
    error_messages,
    boot_meta,
    bootstrap_draw
  )
  gc()
}
plan(sequential)  # shut down sequential session (by-seed) that at every seed runs in parallel

# ----- metrics of bootstrapping -----
library(dplyr)

cell_type = "cd8"
B = 20
p = 2
alpha = 0.05
top_k_values = c(50, 100, 200)

bootstrap_dir = file.path(source_dir, "bootstrap_results")

pval_dir = file.path(bootstrap_dir, paste0(cell_type, "_bootstrap_pvals"))

stability_dir = file.path(bootstrap_dir,paste0(cell_type, "_bootstrap_stability"))

bootstrap_files = sort(
  list.files(
    pval_dir,
    pattern = paste0("^",cell_type,
      "_bootstrap_[0-9]{3}_p_",p,
      "_bonferroni_pvals[.]csv$"
    ),
    full.names = TRUE
  )
)

stopifnot(length(bootstrap_files) == B)

bootstrap_tables = lapply(bootstrap_files,
  function(file) {
    read.csv(file,stringsAsFactors = FALSE)
  }
)

# Use the first bootstrap to define the gene universe.
gene_universe = bootstrap_tables[[1]]$gene

stopifnot(!anyDuplicated(gene_universe))

# Verify and align every table by gene name.
bootstrap_tables = lapply(bootstrap_tables,
  function(x) {
    stopifnot(setequal(x$gene, gene_universe))
    
    x[match(gene_universe, x$gene),, drop = FALSE]
  }
)

bootstrap_tables

#  ----- discovery tables -----

discovery_counts = vapply(bootstrap_tables,
  function(x) {
    sum(!is.na(x$bonferroni_pval) &  x$bonferroni_pval < alpha)
  },
  integer(1)
)

discovery_count_df = data.frame(
  bootstrap_id = seq_len(B),
  discoveries = discovery_counts
)

discovery_quantiles = quantile(
  discovery_counts,
  probs = c(0.25, 0.50, 0.75),
  names = FALSE
)

discovery_summary = data.frame(
  B = B,
  minimum = min(discovery_counts),
  Q1 = discovery_quantiles[1],
  median = discovery_quantiles[2],
  Q3 = discovery_quantiles[3],
  IQR = IQR(discovery_counts),
  maximum = max(discovery_counts)
)

print(discovery_count_df)
print(discovery_summary)

# ----- bootstrap selection frequencies -----
selection_matrix = vapply(bootstrap_tables,
  function(x) {
    as.integer(!is.na(x$bonferroni_pval) & x$bonferroni_pval < alpha)
  },
  integer(length(gene_universe))
)

rownames(selection_matrix) = gene_universe
colnames(selection_matrix) = sprintf("bootstrap_%03d", seq_len(B))

selection_count = rowSums(selection_matrix)

selection_frequency = selection_count/B

selection_frequency_df = data.frame(
  gene = gene_universe,
  selected_n = selection_count,
  selection_frequency = selection_frequency,
  stringsAsFactors = FALSE
) %>%
  arrange(
    desc(selection_frequency),
    gene
  )

head(selection_frequency_df)

write.csv(selection_frequency_df,
  file.path(stability_dir, paste0(cell_type, "_bootstrap_selection_frequencies.csv")),
  row.names = FALSE
)

# ----- frequency threshold summary -----
selection_threshold_summary <- data.frame(
  threshold = c(
    "20/20",
    "18/20 or more",
    "16/20 or more",
    "10/20 or more",
    "At least once"
  ),
  n_genes = c(
    sum(selection_count == 20),
    sum(selection_count >= 18),
    sum(selection_count >= 16),
    sum(selection_count >= 10),
    sum(selection_count >= 1)
  )
)
print(selection_threshold_summary)

# ----- get raw-pvals for stability ranking
full_result_file = paste0("/Users/zaqian/Desktop/density_estimation/", "JASA_revision_3_09_2026/", "with_revised_covariance_RDA_5_1_2026/", "res_cd8_p_2_sef_5_1_2026.RDS")
full_results = readRDS(full_result_file)

stopifnot(!is.null(names(full_results)),!anyDuplicated(names(full_results)))

full_pvalues = vapply(full_results,
  function(x) {
    if (!is.list(x) || is.null(x$pval_1) || length(x$pval_1) != 1L) {
      return(NA_real_)
    }
    
    as.numeric(x$pval_1)
  },
  numeric(1)
)

names(full_pvalues) = names(full_results)

stopifnot(setequal(names(full_pvalues), gene_universe))

full_pvalues <- full_pvalues[gene_universe]  #put OG full-data p-values in the same gene order

stopifnot(!anyNA(full_pvalues), all(is.finite(full_pvalues)))

# ----- get ranks across entire datasets  -----
full_ranks = rank(full_pvalues,ties.method = "average") # get ranks  and use avg for the pvals as tiebreaker
names(full_ranks) = gene_universe

get_top_k = function(genes, pvalues, k) {
  stopifnot(k <= length(genes))
  # Gene name is a reproducible secondary ordering if a tie occurs exactly at the top-k boundary
  ordering = order(
    pvalues,
    genes,
    na.last = NA
  )
  return(genes[ordering[seq_len(k)]])
}

full_top_sets = lapply(top_k_values,
  function(k) {
    get_top_k(genes = gene_universe, pvalues = full_pvalues, k = k)
  }
)

names(full_top_sets) = paste0("top_",top_k_values)

# ----- rank Spearman correlation and top-k overlap for every bootstrap -----
ranking_results = vector("list", B)

for (b in seq_len(B)) {
  bootstrap_pvalues = bootstrap_tables[[b]]$p_value
  
  stopifnot(!anyNA(bootstrap_pvalues), all(is.finite(bootstrap_pvalues))) # no NA or infinite vals
  
  bootstrap_ranks = rank(
    bootstrap_pvalues,
    ties.method = "average"
  )
  

  rho_b = cor( # pearson correlation on the RANK vectors is identical to Spearman rank corr
    full_ranks,
    bootstrap_ranks,
    method = "pearson"
  )
  
  result_b = data.frame(
    bootstrap_id = b,
    spearman_rho = rho_b
  )
  
  for (k in top_k_values) {
    bootstrap_top = get_top_k(
      genes = gene_universe,
      pvalues = bootstrap_pvalues,
      k = k
    )
    
    overlap = length(intersect(full_top_sets[[paste0("top_", k)]], bootstrap_top)) / k
    
    result_b[[paste0("top_", k, "_overlap")]] = overlap
  }
  
  ranking_results[[b]] = result_b
}

ranking_results = bind_rows(ranking_results)

print(ranking_results)

summarize_metric = function(x) {
  q = quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    names = FALSE
  )
  
  return(data.frame(
    median = q[2],
    Q1 = q[1],
    Q3 = q[3],
    IQR = q[3] - q[1],
    minimum = min(x),
    maximum = max(x)
  ))
}

rank_correlation_summary = cbind(
  metric = "Spearman Rank Correlation",
  summarize_metric(ranking_results$spearman_rho)
)

print(rank_correlation_summary)

# ----- top-k summaries -----
top_k_summary = bind_rows(
  lapply(
    top_k_values,
    function(k) {
      metric_name = paste0("top_",k,"_overlap")
      cbind(
        k = k,
        summarize_metric(ranking_results[[metric_name]])
      )
    }
  )
)

print(top_k_summary)



# ----- save remaining results -----
write.csv(
  ranking_results,
  file.path(stability_dir, paste0(cell_type,"_bootstrap_ranking_results.csv")),
  row.names = FALSE
)

write.csv(
  rank_correlation_summary,
  file.path(stability_dir,  paste0(cell_type, "_bootstrap_rank_correlation_summary.csv")),
  row.names = FALSE
)

write.csv(
  top_k_summary,
  file.path(stability_dir,  paste0(cell_type, "_bootstrap_top_k_overlap_summary.csv")),
  row.names = FALSE
)
