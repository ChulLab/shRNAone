#!/usr/bin/env Rscript
# partE.02.clusters_acc_compare.R
#

library(data.table)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_CLUSTER <- file.path(PROJECT, "result/08.manuscript_03/partE/01.spatial_clusters/01.Spatial_Clusters_Details.csv")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partE/02.clusters_acc")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MAT_PATHS <- list(
  TTR_W35 = file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix/TTR_shRNAone_TTR_based_W35_L20_U7.csv"),
  TTR_W80 = file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix/TTR_shRNAone_RNAxs_based_W80_L40_U7.csv"),
  PCSK9_W35 = file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix/PCSK9_shRNAone_TTR_based_W35_L20_U7.csv"),
  PCSK9_W80 = file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix/PCSK9_shRNAone_RNAxs_based_W80_L40_U7.csv")
)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. 加载数据 =================
log_msg("Loading cluster labels and accessibility matrices...")
cluster_dt <- fread(IN_CLUSTER)

load_acc <- function(path, param_name, gene_name) {
  dt <- fread(path)[, .(start_pos, Accessibility)]
  dt[, Parameter := param_name]
  dt[, Gene := gene_name]
  return(dt)
}

acc_all <- rbindlist(list(
  load_acc(MAT_PATHS$TTR_W35, "W35_L20_U7", "TTR"),
  load_acc(MAT_PATHS$TTR_W80, "W80_L40_U7", "TTR"),
  load_acc(MAT_PATHS$PCSK9_W35, "W35_L20_U7", "PCSK9"),
  load_acc(MAT_PATHS$PCSK9_W80, "W80_L40_U7", "PCSK9")
))

merged_dt <- merge(acc_all, cluster_dt[, .(Gene, start_pos, Region_Type, Window_NT)], 
                   by = c("Gene", "start_pos"), all.x = FALSE, allow.cartesian = TRUE)

# 严格维持 Optimal -> Background -> Poor 逻辑顺序
merged_dt$Region_Type <- factor(merged_dt$Region_Type, levels = c("Optimal", "Background", "Poor"))
fwrite(merged_dt, file.path(OUT_DIR, "01.Merged_Clusters_Acc_AllWindows.csv"))

# ================= 2. 积木式组装与画图 =================
log_msg("Generating synchronized Boxplots across all windows...")

windows <- sort(unique(merged_dt$Window_NT))
all_rows <- list()

# 定义每一行四个分面的顺序
plot_configs <- list(
  list(param = "W35_L20_U7", gene = "TTR"),
  list(param = "W35_L20_U7", gene = "PCSK9"),
  list(param = "W80_L40_U7", gene = "TTR"),
  list(param = "W80_L40_U7", gene = "PCSK9")
)

for (w in windows) {
  window_dt <- merged_dt[Window_NT == w]
  plot_list_for_w <- list()
  
  # 提取当前 Window 下所有亚组的最大非异常值 (Upper Whisker)
  whisker_dt <- window_dt[, .(max_whisker = {
    x <- na.omit(Accessibility)
    if (length(x) < 5) max(c(x, 0)) else boxplot.stats(x)$stats[5]
  }), by = .(Parameter, Gene, Region_Type)]
  
  global_max <- max(whisker_dt$max_whisker, na.rm = TRUE)
  if (is.infinite(global_max) || is.na(global_max) || global_max <= 0) global_max <- 0.1
  
  # 全局统一 Y 轴上限：留出充足空间画三条显著性横线
  y_limit <- global_max * 1.8 
  y_pos_1 <- global_max * 1.15  # 第一条 (Optimal vs Background)
  y_pos_2 <- global_max * 1.35  # 第二条 (Poor vs Background)
  y_pos_3 <- global_max * 1.60  # 第三条 (Optimal vs Poor) 位于最顶端
  
  for (cfg in plot_configs) {
    sub_dt <- window_dt[Parameter == cfg$param & Gene == cfg$gene]
    
    # 动态统计数量 n，拼接标签
    sub_dt[, n := .N, by = Region_Type]
    sub_dt[, Region_Label := sprintf("%s\n(n=%d)", Region_Type, n)]
    
    # 保持因子顺序
    ordered_labels <- unique(sub_dt[order(Region_Type)]$Region_Label)
    sub_dt$Region_Label <- factor(sub_dt$Region_Label, levels = ordered_labels)
    
    # 精准提取名称避免自动匹配错误
    lbl_opt  <- sub_dt[Region_Type == "Optimal", unique(as.character(Region_Label))]
    lbl_bg   <- sub_dt[Region_Type == "Background", unique(as.character(Region_Label))]
    lbl_poor <- sub_dt[Region_Type == "Poor", unique(as.character(Region_Label))]
    
    my_comparisons <- list()
    y_positions <- numeric()
    
    # 动态装载显著性组并强制分配对应的高度
    if (length(lbl_opt) > 0 && length(lbl_bg) > 0 && sum(sub_dt$Region_Type == "Optimal") >= 2 && sum(sub_dt$Region_Type == "Background") >= 2) {
      my_comparisons[[length(my_comparisons) + 1]] <- c(lbl_opt, lbl_bg)
      y_positions <- c(y_positions, y_pos_1)
    }
    if (length(lbl_poor) > 0 && length(lbl_bg) > 0 && sum(sub_dt$Region_Type == "Poor") >= 2 && sum(sub_dt$Region_Type == "Background") >= 2) {
      my_comparisons[[length(my_comparisons) + 1]] <- c(lbl_poor, lbl_bg)
      y_positions <- c(y_positions, y_pos_2)
    }
    # 新增：Optimal vs Poor 比较
    if (length(lbl_opt) > 0 && length(lbl_poor) > 0 && sum(sub_dt$Region_Type == "Optimal") >= 2 && sum(sub_dt$Region_Type == "Poor") >= 2) {
      my_comparisons[[length(my_comparisons) + 1]] <- c(lbl_opt, lbl_poor)
      y_positions <- c(y_positions, y_pos_3)
    }
    
    # 构建基础 Boxplot，应用全局统一的 coord_cartesian
    p <- ggplot(sub_dt, aes(x = Region_Label, y = Accessibility, fill = Region_Type)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 0.4) +
      scale_fill_manual(values = c("Optimal" = "#2ECC71", "Background" = "#95A5A6", "Poor" = "#C39BD3")) +
      coord_cartesian(ylim = c(0, y_limit)) +
      labs(
        title = sprintf("%s | %s\n(Window = %d nt)", cfg$param, cfg$gene, w),
        x = NULL,
        y = if(cfg$param == "W35_L20_U7" && cfg$gene == "TTR") "Accessibility Score" else NULL
      ) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
        axis.text.x = element_text(face = "plain", color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )
    
    # 追加配对检验与显著性星号
    if (length(my_comparisons) > 0) {
      p <- p + stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                                  method = "wilcox.test", size = 3.5, 
                                  y.position = y_positions, tip.length = 0.02)
    }
    
    plot_list_for_w[[length(plot_list_for_w) + 1]] <- p
  }
  
  # 拼接这一行的 4 列
  row_plot <- ggarrange(plotlist = plot_list_for_w, ncol = 4, nrow = 1)
  all_rows[[length(all_rows) + 1]] <- row_plot
}

# ================= 3. 总体拼接输出 =================
log_msg("Stitching all windows vertically into a single seamless PDF...")
final_combined <- ggarrange(plotlist = all_rows, ncol = 1, nrow = length(windows))

out_pdf <- file.path(OUT_DIR, "02.Clusters_Acc_Comparison_All.pdf")
ggsave(out_pdf, plot = final_combined, width = 8.5, height = 2.5 * length(windows), limitsize = FALSE)

log_msg("Part E.02 All Done. Unified scales, clean labels, and complete significance markers applied.")