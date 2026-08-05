#!/usr/bin/env Rscript
# partE.01.spatial_clusters_profiling.R
#

library(data.table)
library(ggplot2)
library(zoo)
library(parallel)
library(ggpubr) # 改用 ggpubr 规避拼接报错

# ================= 0. 配置 =================
PROJECT <- "/data/cai801/data/HKUcas"
WET_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partE/01.spatial_clusters")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

WINDOWS_TO_TEST <- seq(9, 36, by = 3)
N_ITERATIONS <- 1000

set.seed(42)
log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

roll_mean_partial <- function(x, w) {
  rollapply(x, width = w, FUN = mean, na.rm = TRUE, align = "center", partial = TRUE)
}

# ================= 1. 加载数据 =================
log_msg("Step 1: Loading raw wet-lab data...")
df_TTR <- fread(file.path(WET_DIR, "TTR_shRNAone_clean.csv"))[, .(Gene = "TTR", start_pos, log2FC = mean_log2FC)]
df_PCSK9 <- fread(file.path(WET_DIR, "PCSK9_shRNAone_3nt_clean.csv"))[, .(Gene = "PCSK9", start_pos, log2FC = mean_log2FC)]

df_raw <- rbind(df_TTR, df_PCSK9)[!is.na(start_pos) & !is.na(log2FC)]
setorder(df_raw, Gene, start_pos)

# ================= 2. 核心加速引擎 =================
num_cores <- max(1, detectCores() - 2)
log_msg(sprintf("Step 2: Activating parallel permutation engine with %d CPU cores...", num_cores))

run_spatial_permutation <- function(dt, window_nt) {
  res_list <- list()
  
  for (gene in c("TTR", "PCSK9")) {
    gene_dt <- dt[Gene == gene]
    raw_vals <- gene_dt$log2FC
    n <- length(raw_vals)
    if (n == 0) next
    
    # 动态推断靶点坐标间距 (Resolution)
    diffs <- diff(gene_dt$start_pos)
    res <- if (length(diffs[diffs > 0]) > 0) median(diffs[diffs > 0]) else 1
    
    # 将实际物理窗口长度转化为数组步长 (TTR=1nt, PCSK9=3nt)
    window_rows <- max(1, round(window_nt / res))
    log_msg(sprintf("  [%s] Resolution: %d nt | Window: %d nt -> Averaging %d targets", gene, res, window_nt, window_rows))
    
    obs_rolling <- roll_mean_partial(raw_vals, window_rows)
    gene_dt[, Rolling_log2FC := obs_rolling]
    
    # 多线程并发计算 1000 次 Null Distribution
    sim_list <- mclapply(1:N_ITERATIONS, function(i) {
      roll_mean_partial(sample(raw_vals), window_rows)
    }, mc.cores = num_cores)
    
    simulated_means <- do.call(rbind, sim_list)
    
    upper_bound <- apply(simulated_means, 2, quantile, probs = 0.995, na.rm = TRUE)
    lower_bound <- apply(simulated_means, 2, quantile, probs = 0.005, na.rm = TRUE)
    
    gene_dt[, CI_995_Upper := upper_bound]
    gene_dt[, CI_005_Lower := lower_bound]
    
    # 逻辑修正：log2FC 越高切割越好
    gene_dt[, Region_Type := "Background"]
    gene_dt[Rolling_log2FC > upper_bound, Region_Type := "Optimal"]
    gene_dt[Rolling_log2FC < lower_bound, Region_Type := "Poor"]
    
    gene_dt[, Window_NT := window_nt]
    gene_dt[, Window_Rows := window_rows]
    
    res_list[[gene]] <- gene_dt
  }
  return(rbindlist(res_list))
}

# ================= 3. 批量执行与超长画卷输出 =================
all_results <- list()
all_plots <- list()

for (w in WINDOWS_TO_TEST) {
  log_msg(sprintf("Processing Window %d nt...", w))
  res_dt <- run_spatial_permutation(df_raw, w)
  all_results[[as.character(w)]] <- res_dt
  
  for (g in c("TTR", "PCSK9")) {
    plot_dt <- res_dt[Gene == g]
    plot_dt[, Rel_Pos := (start_pos - min(start_pos)) / (max(start_pos) - min(start_pos)) * 100]
    
    # 无缝带状图坐标处理
    plot_dt[, Opt_Ymin := CI_995_Upper]
    plot_dt[, Opt_Ymax := pmax(Rolling_log2FC, CI_995_Upper)]
    plot_dt[, Poor_Ymax := CI_005_Lower]
    plot_dt[, Poor_Ymin := pmin(Rolling_log2FC, CI_005_Lower)]
    
    p <- ggplot(plot_dt, aes(x = Rel_Pos)) +
      geom_ribbon(aes(ymin = CI_005_Lower, ymax = CI_995_Upper), fill = "gray60", alpha = 0.3) +
      geom_line(aes(y = CI_995_Upper), linetype = "dashed", linewidth = 0.4) +
      geom_line(aes(y = CI_005_Lower), linetype = "dashed", linewidth = 0.4) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      geom_ribbon(aes(ymin = Opt_Ymin, ymax = Opt_Ymax), fill = "#2ECC71", alpha = 0.7) +
      geom_ribbon(aes(ymin = Poor_Ymin, ymax = Poor_Ymax), fill = "#C39BD3", alpha = 0.7) +
      geom_line(aes(y = Rolling_log2FC), color = "#1A5276", linewidth = 0.8) +
      scale_x_continuous(limits = c(-2, 102), expand = c(0, 0)) +
      labs(
        title = sprintf("%s (Window = %d nt)", g, w),
        x = "Relative Transcript Position (%)",
        y = bquote("Rolling Mean"~log[2]~"FC")
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")
      )
    # 修复：使用无名索引存入列表，避免底层命名空间冲突
    all_plots[[length(all_plots) + 1]] <- p
  }
}

final_df <- rbindlist(all_results)
fwrite(final_df, file.path(OUT_DIR, "01.Spatial_Clusters_Details.csv"))

log_msg("Step 4: Stitching plots together...")
# 使用 ggarrange 拼接，设定 2 列 (左 TTR, 右 PCSK9)
combined_plot <- ggarrange(plotlist = all_plots, ncol = 2, nrow = length(WINDOWS_TO_TEST))

ggsave(file.path(OUT_DIR, "02.Spatial_Permutation_Tracks_All.pdf"), 
       plot = combined_plot, 
       width = 8, 
       height = 2.5 * length(WINDOWS_TO_TEST),
       limitsize = FALSE)

log_msg("Part E.01 All Done.")