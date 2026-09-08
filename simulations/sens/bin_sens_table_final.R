rm(list = ls())  
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
# ----- directory -----
out.dir = file.path("")

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



# ----- get original outputs -----
fd_mean = readRDS(file.path("NB_scaled_strong_effect_mean_sef_metrics_list.RDS")) # mean-shift results from IS source file
fd_var = readRDS(file.path("ls_NB_scaled_small_effect_var_sef_metrics_list.RDS")) # variance-shift results from IS source file
fd_modal = readRDS(file.path("NB_modal_sef_scaled_modal_metrics_list.RDS")) # modality-shift results from IS source file

fd_mean_avg  = get_avg_metrics(fd_mean, c("p2", "p3", "p4"))
fd_mean_avg = cbind(method = paste0("SEF (p = ", sub("^p", "", fd_mean_avg$moment),")"), fd_mean_avg)
fd_mean_avg = subset(fd_mean_avg, alpha == 0.05)

fd_var_avg = get_avg_metrics(fd_var, c("p2", "p3", "p4"))
fd_var_avg = cbind(method = paste0("SEF (p = ", sub("^p", "", fd_var_avg$moment),")"), fd_var_avg)
fd_var_avg = subset(fd_var_avg, alpha == 0.05)

fd_modal_avg = get_avg_metrics(fd_modal, c("p2","p3", "p4", "p5"))
fd_modal_avg = cbind(method = paste0("SEF (p = ", sub("^p", "", fd_modal_avg$moment),")"), fd_modal_avg)
fd_modal_avg = subset(fd_modal_avg, alpha == 0.05)

# ----- get outputs for varying K -----
Ks =  c(50,100,200)
file_templates = c(
  modal = "K_%d_NB_modal_sef_scaled_modal_metrics_list.RDS",
  variance = "K_%d_NB_scaled_small_effect_var_sef_metrics_list.RDS",
  mean = "K_%d_NB_scaled_strong_effect_mean_sef_metrics_list.RDS"
)
results = setNames(
  lapply(Ks, function(K) {
    lapply(file_templates, function(template) {
      readRDS(file.path(out.dir, sprintf(template, K)))
    })
  }),
  paste0("K_", Ks)
)
result_names = c(
  mean = "mean",
  variance = "variance",
  modality = "modal"
)

p_tested = list(
  mean = c("p2", "p3", "p4"),
  variance = c("p2", "p3", "p4"),
  modality = c("p2", "p3", "p4", "p5")
)

combined_by_K = list()

for (K in Ks) {
  K_name = paste0("K_", K)
  
  metrics_by_type = lapply(names(result_names), function(shift_type) {
    result_name = result_names[[shift_type]]
    temp = subset(
      get_avg_metrics(results[[K_name]][[result_name]],p_tested[[shift_type]]),
      alpha == 0.05 # we only care about alpha of 0.05
    )
    
    data.frame(
      K = K,
      shift_type = shift_type,
      method = paste0("SEF (p = ", sub("^p", "", temp$moment),")"),
      temp,
      row.names = NULL
    )
  })
  
  combined = do.call(rbind, metrics_by_type)
  rownames(combined) = NULL
  
  combined_by_K[[K_name]] = combined
  
  write.csv(combined, file.path(out.dir, sprintf("K_%d_avg_metrics.csv", K)), row.names = FALSE)
}

# ----- combine and save FD and fixed-K averages -----

fd_combined = bind_rows(
  fd_mean_avg |>
    mutate(K = "FD", shift_type = "mean", .before = 1),
  
  fd_var_avg |>
    mutate(K = "FD", shift_type = "variance", .before = 1),
  
  fd_modal_avg |>
    mutate(K = "FD", shift_type = "modality", .before = 1)
)

fixed_K_combined = bind_rows(combined_by_K) |>
  mutate(K = as.character(K))

all_K_metrics = bind_rows(
  fd_combined,
  fixed_K_combined
) |>
  mutate(
    K_setting = if_else(K == "FD", "FD", paste0("K_", K)),
    K_setting = factor(
      K_setting,
      levels = c("FD", paste0("K_", Ks))
    ),
    shift_type = factor( # add shift column
      shift_type,
      levels = c("mean", "variance", "modality")
    ),
    method = factor( # add the neat method version
      method,
      levels = paste0("SEF (p = ", 2:5, ")")
    )
  )

write.csv(
  all_K_metrics,
  file.path(out.dir, "FD_K_50_K_100_K_200_avg_metrics.csv"),
  row.names = FALSE
)

# ----- fig -----
# colors metric determines hue
# bin resolution K determines shade
fill_colors = c(
  "TPR__K_50"  = "#234EA5",
  "TPR__K_100" = "#3568D4",
  "TPR__K_200" = "#4B7BEC",
  "TPR__FD"    = "#78A0F5",
  
  "FDR__K_50"  = "#C86A00",
  "FDR__K_100" = "#E38200",
  "FDR__K_200" = "#F39C12",
  "FDR__FD"    = "#F7B955"
)

fill_labels = c(
  "TPR__K_50"  = "TPR: K = 50",
  "TPR__K_100" = "TPR: K = 100",
  "TPR__K_200" = "TPR: K = 200",
  "TPR__FD"    = "TPR: FD",
  
  "FDR__K_50"  = "FDR: K = 50",
  "FDR__K_100" = "FDR: K = 100",
  "FDR__K_200" = "FDR: K = 200",
  "FDR__FD"    = "FDR: FD"
)

fill_levels = c(  # forces an interleaved order to make sure that FDR is on left and TPR is on right within each K group
  "FDR__K_50",  "TPR__K_50",
  "FDR__K_100", "TPR__K_100",
  "FDR__K_200", "TPR__K_200",
  "FDR__FD",    "TPR__FD"
)

legend_order = c(
  "TPR__K_50",
  "TPR__K_100",
  "TPR__K_200",
  "TPR__FD",
  "FDR__K_50",
  "FDR__K_100",
  "FDR__K_200",
  "FDR__FD"
)

make_shift_plot = function(
    data,
    selected_shift,
    selected_moments,
    plot_title
) {
  
  moment_labels = paste0(
    "p = ",
    sub("^p", "", selected_moments)
  )
  
  temp = data |>
    filter(
      shift_type == selected_shift,
      moment %in% selected_moments
    ) |>
    pivot_longer(
      cols = c(TPR, FDR),
      names_to = "Metric",
      values_to = "Value"
    ) |>
    mutate(
      Metric = factor(
        Metric,
        levels = c("FDR", "TPR")
      ),
      moment_label = factor(
        moment,
        levels = selected_moments,
        labels = moment_labels
      ),
      fill_key = factor(
        paste(Metric, K_setting, sep = "__"),
        levels = fill_levels
      )
    )
  
  ggplot(
    temp,
    aes(
      x = K_setting,
      y = Value,
      fill = fill_key
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.82),
      width = 0.72
    ) +
    geom_hline(
      yintercept = 0.80,
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
    facet_wrap(
      ~ moment_label,
      nrow = 1
    ) +
    scale_x_discrete(
      labels = c(
        "K_50" = "50",
        "K_100" = "100",
        "K_200" = "200",
        "FD" = "FD"
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_fill_manual(
      values = fill_colors,
      breaks = legend_order,
      labels = fill_labels[legend_order],
      drop = FALSE,
      guide = guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    ) +
    labs(
      title = plot_title,
      x = "Bin Resolution",
      y = NULL,
      fill = "Metric and Bin Resolution"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

mean_plot = make_shift_plot(
  data = all_K_metrics,
  selected_shift = "mean",
  selected_moments = c("p2", "p3", "p4"),
  plot_title = "Mean Shift"
)

variance_plot = make_shift_plot(
  data = all_K_metrics,
  selected_shift = "variance",
  selected_moments = c("p2", "p3", "p4"),
  plot_title = "Variance Shift"
)

modality_plot = make_shift_plot(
  data = all_K_metrics,
  selected_shift = "modality",
  selected_moments = c("p3", "p4", "p5"),
  plot_title = "Modality Shift"
)

combined_plot = wrap_plots(
  mean_plot,
  variance_plot,
  modality_plot,
  ncol = 1,
  guides = "collect"
) &
  theme(legend.position = "bottom")

combined_plot

ggsave(
  filename = file.path(
    out.dir,
    "K_sensitivity_combined_TPR_FDR.pdf"
  ),
  plot = combined_plot,
  width = 14,
  height = 12
)
