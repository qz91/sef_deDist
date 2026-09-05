# qq plot for sim_obj across different K bin counts
rm(list = ls()) 
library(pbapply)
library(pbmcapply)
source("revisions_simulation_main_5_8_2026.R")
source("/Users/zaqian/Desktop/density_estimation/JASA_revision_3_09_2026/simulation_revisions/modality/main_NB_4_30_2026.R")

# ----- params -----
n1 = 100
n2 = 100
lower = 500
upper = 1000
n_nonDE_unimodal = 1
n_nonDE_mixture = 1
n_DE = 0
nonDE_unimodal_mu_range = c(5, 20)
nonDE_unimodal_theta_range = c(3, 10)
nonDE_mixture_mu_range = c(5, 15)
nonDE_mixture_theta_range = c(3, 10)
nonDE_mixture_fc_range = c(2, 4)
nonDE_mixture_prob_range = c(0.30, 0.70)
DE_mu_range = c(5, 15)
DE_theta_range = c(3, 10)
DE_mixture_fc_range = c(2, 4)
DE_mixture_prob_range = c(0.5, 0.5) # enforces that there is a mean matching
sigma_sq_nonDE = 0.05
sigma_sq = 0.05 # sets them equal
sigma_sq_DE = 0.05
de_unimodal_group = "group1"
share_mode = F # DE unimodal case in group 1 has mean = (mu1 + mu2)/2 for the mixture DE

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")

# ----- computational functions used -----
empirical_reject_rate = function(qq_out, alpha = 0.05) {
  do.call(rbind, lapply(names(qq_out$pVecs), function(p_name) {
    pvals <- unlist(qq_out$pVecs[[p_name]], use.names = FALSE)
    
    data.frame(
      p = p_name,
      alpha = alpha,
      n_tested = sum(!is.na(pvals)),
      n_failed = sum(is.na(pvals)),
      reject_rate = mean(pvals < alpha, na.rm = TRUE)
    )
  }))
}

qq_run_one_sef_sens = function(seed = 1, ridge_run = 1e-8, ps = c(2, 3, 4), K = NULL, bw = 2, qq_type = c("scaled", "raw")) {
  qq_type = match.arg(qq_type)
  sef_test_fun = switch(
    qq_type,
    scaled = run_sef_test_covfix_scaled_basis,
    raw = run_sef_test_covfix
  )
  sim_obj <- generate_NB_gene_data(
    n1 = n1,
    n2 = n2,
    lower = lower,
    upper = upper,
    n_nonDE_unimodal = n_nonDE_unimodal,
    n_nonDE_mixture = n_nonDE_mixture,
    n_DE = n_DE,
    nonDE_unimodal_mu_range = nonDE_unimodal_mu_range,
    nonDE_unimodal_theta_range = nonDE_unimodal_theta_range,
    nonDE_mixture_mu_range = nonDE_mixture_mu_range,
    nonDE_mixture_theta_range = nonDE_mixture_theta_range,
    nonDE_mixture_fc_range = nonDE_mixture_fc_range,
    nonDE_mixture_prob_range = nonDE_mixture_prob_range,
    DE_mu_range = DE_mu_range,
    DE_theta_range = DE_theta_range,
    DE_mixture_fc_range = DE_mixture_fc_range,
    DE_mixture_prob_range = DE_mixture_prob_range,
    sigma_sq_nonDE = sigma_sq_nonDE,
    sigma_sq_DE = sigma_sq_DE,
    sigma_sq = sigma_sq,
    de_unimodal_group = de_unimodal_group,
    share_mode = share_mode,
    seed = seed
  )
  
  results_by_p = setNames(lapply(ps, function(p) {
    test = sef_test_fun(sim_obj = sim_obj, K = K, p = p, verbose = FALSE, bw = bw, ridge = ridge_run)
    
    list(
      pVec = test$pVec,
      failed = is.na(test$pVec)
    )
  }), paste0("p", ps))
  
  list(
    seed = seed,
    K = K,
    qq_type = qq_type,
    results = results_by_p
  )
}

run_qq_sims_parallel_sens = function(n_sim, ridge_run = 1e-8, K = NULL, bw = 2,
                                      ps = c(2, 3, 4),
                                      n_cores = 3,
                                      qq_type = c("scaled", "raw")) {
  qq_type = match.arg(qq_type)
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  future::plan(future::multisession, workers = n_cores)
  
  seeds <- seq_len(n_sim)
  p_names <- paste0("p", ps)
  
  p <- progressr::progressor(along = seeds)
  
  sim_results <- future.apply::future_lapply(
    seeds,
    function(seed) {
      out <- qq_run_one_sef_sens(
        seed = seed,
        ridge_run = ridge_run,
        ps = ps,
        K = K,
        bw = bw,
        qq_type = qq_type
      )
      
      p(sprintf("seed %d", seed))
      
      out
    },
    future.seed = TRUE
  )
  
  pVecs_by_p <- setNames(lapply(p_names, function(p_name) {
    lapply(sim_results, function(x) x$results[[p_name]]$pVec)
  }), p_names)
  
  failed_by_p <- setNames(lapply(p_names, function(p_name) {
    lapply(sim_results, function(x) x$results[[p_name]]$failed)
  }), p_names)
  
  list(
    seeds = seeds,
    K = K,
    bw = bw,
    qq_type = qq_type,
    pVecs = pVecs_by_p,
    failed = failed_by_p
  )
}
# ----- figure generating functions used -----
make_header_plot = function(label, angle = 0, size = 5, parse = FALSE) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      fontface = "bold",
      size = size,
      angle = angle,
      parse = parse
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}
make_axis_label_plot = function(label, angle = 0, size = 4) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = size,
      angle = angle
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}
make_qq_panel = function(qq_out, p_name) {
  pvals = unlist(qq_out$pVecs[[p_name]], use.names = FALSE)
  pvals = pvals[is.finite(pvals) & !is.na(pvals)]
  
  gg_qqplot(pvals) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}
make_K_sens_qq_patchwork = function(qq_by_K, ps = c(2, 3, 4)) {
  p_names = paste0("p", ps)
  K_names = names(qq_by_K)
  K_labels = c(
    K50 = "K == 50",
    K100 = "K == 100",
    K200 = "K == 200",
    FD = "K[FD]"
  )
  
  header_row = wrap_plots(
    c(
      list(make_header_plot("")),
      lapply(ps, function(p) {
        make_header_plot(
          paste0("p == ", p),
          size = 5,
          parse = TRUE
        )
      })
    ),
    nrow = 1,
    widths = c(0.18, rep(1, length(ps)))
  )
  
  plot_rows = lapply(seq_along(qq_by_K), function(i) {
    K_label = K_names[i]
    
    row_header = make_header_plot(
      K_labels[[K_label]],
      angle = 90,
      size = 5,
      parse = TRUE
    )
    
    panels = lapply(p_names, function(p_name) {
      make_qq_panel(
        qq_out = qq_by_K[[i]],
        p_name = p_name
      )
    })
    
    wrap_plots(
      c(list(row_header), panels),
      nrow = 1,
      widths = c(0.18, rep(1, length(ps)))
    )
  })
  
  core_grid = wrap_plots(
    c(list(header_row), plot_rows),
    ncol = 1,
    heights = c(0.16, rep(1, length(plot_rows)))
  )
  
  y_label = make_axis_label_plot(
    expression(paste("Observed -log"[10], plain(P))),
    angle = 90,
    size = 5
  )
  
  x_label = make_axis_label_plot(
    expression(paste("Expected -log"[10], plain(P))),
    angle = 0,
    size = 5
  )
  
  main = y_label | core_grid
  main = main + plot_layout(widths = c(0.06, 1))
  
  bottom = plot_spacer() | x_label
  bottom = bottom + plot_layout(widths = c(0.06, 1))
  
  main / bottom +
    plot_layout(heights = c(1, 0.06))
}
qq_scaled_grid = make_K_sens_qq_patchwork(
  qq_by_K = qq_by_K,
  ps = c(2, 3, 4)
)
# ----- simulate gene density estimates with various K -----
sens_dir = file.path("/Users/zaqian/Desktop/finalSims/NB_sens")
Ks = list( K50 = 50, K100 = 100, K200 = 200, FD = NULL)
n_sims = 500
qq_by_K = lapply(names(Ks), function(K_name) {
  message("Running ", K_name)
  
  progressr::with_progress({
    run_qq_sims_parallel_sens(
      n_sim = n_sims,
      ps = c(2, 3, 4),
      n_cores = 3,
      K = Ks[[K_name]],
      bw = 2,
      qq_type = "scaled"
    )
  })
})

names(qq_by_K) = names(Ks)
# saveRDS(qq_by_K, file.path(sens_dir,"Ktest_pvecs_scaled_NB.rds"))
qq_by_K = readRDS(file.path(sens_dir, "Ktest_pvecs_scaled_NB.rds"))
type_i_by_K <- do.call(rbind, lapply(qq_by_K, empirical_reject_rate))

## ----- empirical type i figs using patchwork -----
ggplot2::ggsave(
  file.path(sens_dir,"qq_Ktest_scaled_NB.pdf"),
  qq_scaled_grid,
  width = 15,
  height = 10,
  dpi = 300
)
qq_scaled_grid

# ----- p = 2, 3 moments -----
qq_by_K_raw = lapply(names(Ks), function(K_name) {
  message("Running ", K_name)
  
  progressr::with_progress({
    run_qq_sims_parallel_sens(
      n_sim = n_sims,
      ps = c(2, 3),
      n_cores = 3,
      K = Ks[[K_name]],
      bw = 2,
      qq_type = "raw"
    )
  })
})

names(qq_by_K_raw) = names(Ks)
saveRDS(qq_by_K_raw, file.path(sens_dir,"Ktest_pvecs_raw_NB.rds"))
qq_raw_grid = make_K_sens_qq_patchwork(
  qq_by_K = qq_by_K_raw,
  ps = c(2, 3)
)
## ----- empirical type i figs using patchwork -----
ggplot2::ggsave(
  file.path(sens_dir,"qq_Ktest_raw_NB.pdf"),
  qq_raw_grid,
  width = 15,
  height = 8,
  dpi = 300
)
qq_scaled_grid


