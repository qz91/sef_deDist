library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
# ----- functions to compute averages -----
get_avg_metrics = function(x, p_vals, alpha_vals = c("alpha_0.05", "alpha_0.10"), metrics = c("TPR", "FDR")) {
  do.call(rbind, lapply(p_vals, function(p) {
    do.call(rbind, lapply(alpha_vals, function(a) {
      seed_df <- do.call(rbind, lapply(x, function(seed) {
        data.frame(seed[[p]][[a]])
      }))
      
      avgs <- colMeans(seed_df[, metrics, drop = FALSE], na.rm = TRUE)
      
      data.frame(
        moment = p,
        alpha = sub("alpha_", "", a),
        t(avgs),
        row.names = NULL
      )
    }))
  }))
}

get_avg_pb_metrics = function(x, methods = NULL, metrics = c("TPR", "FDR")) {
  if (is.null(methods)) {
    methods <- names(x[[1]]$results)
  }
  
  out <- do.call(rbind, lapply(methods, function(method) {
    seed_df <- do.call(rbind, lapply(x, function(seed) {
      df <- data.frame(seed$results[[method]]$metrics)
      df[, metrics, drop = FALSE]
    }))
    
    avgs <- colMeans(seed_df, na.rm = TRUE)
    
    data.frame(
      method = method,
      t(avgs),
      row.names = NULL
    )
  }))
  
  rownames(out) <- NULL
  out
}


# ----- mean shift (strong) -----
load.dir = file.path("/Users/zaqian/Desktop/finalSims/NB_strong")
## ----- load all results -----
sef_scaled = readRDS(file.path(load.dir,"sef_strong/NB_scaled_strong_effect_mean_sef_metrics_list.RDS"))
mast_res = readRDS(file.path(load.dir,"NB_MAST_cell_strong_effect_mean_metrics.RDS"))
muscat_res = readRDS(file.path(load.dir,"muscat_strong_mean/NB_muscat_strong_effect_mean_ddf_Satterthwaite_metrics.RDS"))
scDD_res = readRDS(file.path(load.dir,"NB_scDD_strong_effect_mean_metrics.RDS"))
pb_res = readRDS(file.path(load.dir,"NB_strong_mean_shift_avg_pb_metrics_5_22_2026.RDS"))
wilcox_res = readRDS(file.path(load.dir,"wilcox_cell_strong/NB_strong_mean_shift_wilcox_cell_metrics_5_22_2026.RDS"))
wilcox_row = cbind(method = "Wilcox (cell-level)", data.frame(t(colMeans(wilcox_res, na.rm = TRUE))))
scDD_row = cbind(method = "scDD", data.frame(t(colMeans(scDD_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
mast_row = cbind(method = "MAST", data.frame(t(colMeans(mast_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
muscat_row = cbind(method = "muscat (cell-level)", data.frame(t(colMeans(muscat_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))

sef_scaled_avg_metrics = get_avg_metrics(sef_scaled, c("p2", "p3", "p4"))
sef_scaled_avg_metrics = cbind(method = paste0("SEF (p = ", sub("^p", "", sef_scaled_avg_metrics$moment),")"), sef_scaled_avg_metrics)
sef_scaled_avg_metrics = subset(sef_scaled_avg_metrics, alpha == 0.05)

pb_res$method[pb_res$method == "wilcox"] = "Wilcox"
pb_res = pb_res %>%
  mutate(method = paste0(method, " (PB)"))

## ----- make manuscript table -----
strong_mean_avg_combined = dplyr::bind_rows(
  sef_scaled_avg_metrics %>% dplyr::select(method, TPR, FDR),
  muscat_row %>% dplyr::select(method, TPR, FDR),
  scDD_row %>% dplyr::select(method, TPR, FDR),
  mast_row %>% dplyr::select(method, TPR, FDR),
  wilcox_row %>% dplyr::select(method, TPR, FDR),
  pb_res
)

out.dir = file.path("/Users/zaqian/Desktop/finalSims/table_figs")
write.csv(strong_mean_avg_combined, file = file.path(out.dir, "NB_strong_mean_avg_metrics.csv"))

## ----- make plot -----
strong_mean_avg_combined2 = strong_mean_avg_combined

method_order = unique(strong_mean_avg_combined2$method) # keep same order
df_long = strong_mean_avg_combined2 %>%
  mutate(method = factor(method, levels = method_order)) %>%
  pivot_longer(
    cols = c("TPR", "FDR"),
    names_to = "Metric",
    values_to = "Value"
  )

plt = ggplot(df_long, aes(x = method, y = Value, fill = Metric)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0.8,
    linetype = "dashed",
    color = "seagreen3",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  labs(x = "", y = "", fill = "Metric") +
  scale_fill_manual(
    values = c("TPR" = "royalblue1", "FDR" = "orange2")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plt

plt_name = "NB_manu_strong_mean_res.pdf"

# ggsave(
#   file.path(out.dir, plt_name),
#   plot = plt,
#   width = 8,
#   height = 5
# )


# ----- var shift (smaller) -----
## ----- load all results -----
load.dir2 = file.path("/Users/zaqian/Desktop/finalSims/NB_smaller")
sef_scaled = readRDS(file.path(load.dir2,"sef_scaled_smaller/ls_NB_scaled_small_effect_var_sef_metrics_list.RDS"))
muscat_res = readRDS(file.path(load.dir2,"cell_level_smaller/muscat_smaller_var/NB_muscat_small_effect_var_ddf_Satterthwaite_metrics.RDS"))
mast_res = readRDS(file.path(load.dir2,"cell_level_smaller/NB_MAST_cell_small_effect_var_metrics.RDS"))
scDD_res = readRDS(file.path(load.dir2,"cell_level_smaller/NB_scDD_small_effect_var_metrics.RDS"))
pb_res = readRDS(file.path(load.dir2,"pseudobulk_smaller/NB_pseudobulk_small_effect_result_5_19_2026.RDS"))
wilcox_res = readRDS(file.path(load.dir2,"cell_level_smaller/NB_wilcox_cell_small_effect_metrics_5_26_2026.RDS"))
wilcox_row = cbind(method = "Wilcox (cell-level)", data.frame(t(colMeans(wilcox_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
scDD_row = cbind(method = "scDD", data.frame(t(colMeans(scDD_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
mast_row = cbind(method = "MAST", data.frame(t(colMeans(mast_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
muscat_row = cbind(method = "muscat (cell-level)", data.frame(t(colMeans(muscat_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
vartest_metrics = readRDS(file.path(load.dir2, "vartest_smaller/NB_vartest_small_effect_result_5_26_2026.RDS")) # not run on cluster
vartest_avg_metrics = lapply(vartest_metrics, function(x) {
  data.frame(
    TPR = x$metrics$TPR,
    FDR = x$metrics$FDR,
    TP = x$metrics$TP,
    FP = x$metrics$FP,
    FN = x$metrics$FN
  )
}) %>%
  bind_rows() %>%
  summarise(
    TPR = mean(TPR, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    TP = mean(TP, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE)
  )
vartest_avg_metrics = cbind(method = "Wilcox (variance)", vartest_avg_metrics)

sef_scaled_avg_metrics = get_avg_metrics(sef_scaled, c("p2", "p3", "p4"))
sef_scaled_avg_metrics = cbind(method = paste0("SEF (p = ", sub("^p", "", sef_scaled_avg_metrics$moment),")"), sef_scaled_avg_metrics)
sef_scaled_avg_metrics = subset(sef_scaled_avg_metrics, alpha == 0.05)

pb_res = get_avg_pb_metrics(pb_res)
pb_res$method[pb_res$method == "wilcox"] = "Wilcox"
pb_res = pb_res %>%
  mutate(method = paste0(method, " (PB)"))

## ----- make manuscript table -----
small_var_avg_combined = dplyr::bind_rows(
  sef_scaled_avg_metrics %>% dplyr::select(method, TPR, FDR),
  muscat_row %>% dplyr::select(method, TPR, FDR),
  scDD_row %>% dplyr::select(method, TPR, FDR),
  mast_row %>% dplyr::select(method, TPR, FDR),
  wilcox_row %>% dplyr::select(method, TPR, FDR),
  vartest_avg_metrics %>% dplyr::select(method, TPR, FDR),
  pb_res
)

out.dir = file.path("/Users/zaqian/Desktop/finalSims/table_figs")
write.csv(small_var_avg_combined, file = file.path(out.dir, "NB_small_var_avg_metrics.csv"))
## ----- make plot -----
small_var_avg_combined2 = small_var_avg_combined

method_order = unique(small_var_avg_combined2$method) # keep same order
df_long = small_var_avg_combined2 %>%
  mutate(method = factor(method, levels = method_order)) %>%
  pivot_longer(
    cols = c("TPR", "FDR"),
    names_to = "Metric",
    values_to = "Value"
  )

plt = ggplot(df_long, aes(x = method, y = Value, fill = Metric)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_hline(
    yintercept = 0.8,
    linetype = "dashed",
    color = "seagreen3",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  labs(x = "", y = "", fill = "Metric") +
  scale_fill_manual(
    values = c("TPR" = "royalblue1", "FDR" = "orange2")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plt

plt_name = "NB_manu_small_var_res.pdf"

# ggsave(
#   file.path(out.dir, plt_name),
#   plot = plt,
#   width = 8,
#   height = 5
# )
