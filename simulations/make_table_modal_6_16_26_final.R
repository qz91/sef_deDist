rm(list = ls()) 

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

# ----- modal directory -----
# make table for modality analysis
load.dir = file.path("/Users/zaqian/Desktop/finalSims/modal")
## ----- load all results -----
sef_scaled = readRDS(file.path(load.dir,"NB_modal_sef_scaled_modal_metrics_list.RDS"))
mast_res = readRDS(file.path(load.dir,"modal_mast/NB_modal_MAST_cell_metrics.RDS"))
muscat_res = readRDS(file.path(load.dir,"modal_muscat/NB_modal_muscat_ddf_Satterthwaite_metrics.RDS"))
scDD_res = readRDS(file.path(load.dir,"modal_scDD/NB_modal_scDD_metrics.RDS"))
pb_res = readRDS(file.path(load.dir,"NB_modal_avg_pb_metrics_6_16_2026.RDS"))
wilcox_res = readRDS(file.path(load.dir,"modal_wilcox_cell/NB_modal_wilcox_cell_metrics.RDS"))
wilcox_row = cbind(method = "Wilcox (cell-level)", data.frame(t(colMeans(wilcox_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
scDD_row = cbind(method = "scDD", data.frame(t(colMeans(scDD_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
mast_row = cbind(method = "MAST", data.frame(t(colMeans(mast_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))
muscat_row = cbind(method = "muscat (cell-level)", data.frame(t(colMeans(muscat_res[, c("TPR", "FDR", "TP", "FP", "FN")], na.rm = TRUE))))



sef_scaled_avg_metrics = get_avg_metrics(sef_scaled, c("p2","p3", "p4", "p5"))
sef_scaled_avg_metrics = cbind(method = paste0("SEF (p = ", sub("^p", "", sef_scaled_avg_metrics$moment),")"), sef_scaled_avg_metrics)
sef_scaled_avg_metrics = subset(sef_scaled_avg_metrics, alpha == 0.05)


pb_res$method[pb_res$method == "wilcox"] = "Wilcox"
pb_res = pb_res %>%
  mutate(method = paste0(method, " (PB)"))


modal_avg_combined = dplyr::bind_rows(
  sef_scaled_avg_metrics %>% dplyr::select(method, TPR, FDR),
  muscat_row %>% dplyr::select(method, TPR, FDR),
  scDD_row %>% dplyr::select(method, TPR, FDR),
  mast_row %>% dplyr::select(method, TPR, FDR),
  wilcox_row %>% dplyr::select(method, TPR, FDR),
  pb_res
)

out.dir = file.path("/Users/zaqian/Desktop/finalSims/table_figs")
write.csv(modal_avg_combined, file = file.path(out.dir, "NB_modal_avg_metrics.csv"))

## ----- make manuscript-ready plot -----
modal_avg_combined2 = modal_avg_combined

method_order = unique(modal_avg_combined2$method) # keep same order
df_long = modal_avg_combined2 %>%
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

plt_name = "NB_manu_modal_res.pdf"

# ggsave(
#   file.path(out.dir, plt_name),
#   plot = plt,
#   width = 8,
#   height = 5
# )


# ----- make table for p = 2, 3, 4, 5 -----
load.dir = file.path("/Users/zaqian/Desktop/finalSims/modal")
sef_scaled = readRDS(file.path(load.dir,"NB_modal_sef_scaled_modal_metrics_list.RDS"))
sef_discrete = readRDS(file.path(load.dir,"NB_modal_sef_scaled_modal_discrete_metrics_list.RDS"))
sef_scaled_avg_metrics = get_avg_metrics(sef_scaled, c("p2","p3", "p4", "p5"))
sef_scaled_avg_metrics = cbind(method = paste0("SEF (p = ", sub("^p", "", sef_scaled_avg_metrics$moment),")"), sef_scaled_avg_metrics)
sef_scaled_avg_metrics = subset(sef_scaled_avg_metrics, alpha == 0.05)
sef_discrete_avg_metrics = get_avg_metrics(sef_discrete, c("p2","p3", "p4", "p5"))
sef_discrete_avg_metrics = cbind(method = paste0("Discrete SEF (p = ", sub("^p", "", sef_discrete_avg_metrics$moment),")"), sef_discrete_avg_metrics)
sef_discrete_avg_metrics = subset(sef_discrete_avg_metrics, alpha == 0.05)
write.csv(sef_scaled_avg_metrics, "sef_scaled_modal_p2345.csv")
write.csv(sef_discrete_avg_metrics, "sef_discrete_modal_p2345.csv")



