#!/usr/bin/env Rscript
# partD.01.optimal_U_selection.R
#

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partD/01.plateau_compare")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")
REGIONS_TO_EXTRACT <- c("whole", "CDS", "3'UTR")

# --- 核心约束: 继承 partC.02 中严格定义的人为 Plateau ---
U_PLATEAU <- list(
  shRNAone   = c(6, 14),
  PspCas13b  = c(9, 23),
  CasRx_day5 = c(8, 14)
)

# --- 批量比对序列配置 ---
CANDIDATES <- list(
  c(W = 30, L = 20),
  c(W = 35, L = 20),
  c(W = 40, L = 20),
  c(W = 45, L = 20),
  c(W = 50, L = 20)
)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 1. 数据全局载入与预处理 =================
log_msg("Step 1: Loading all grid results and computing significance...")

# 追加载入 spearman_p，用于计算显著性
cols_to_keep <- c("W", "L", "U", "direction", "region", "spearman_r", "spearman_p")
dt_list <- lapply(TOOLS, function(t) {
  fpath <- file.path(IN_DIR, paste0(t, "_grid.csv"))
  if (file.exists(fpath)) {
    dt <- fread(fpath, select = cols_to_keep)
    dt$tool <- t
    return(dt)
  }
})
all_dt <- rbindlist(dt_list)

# 增加一列 is_sig 判定显著性 (p <= 0.05 且正相关)
all_dt[, is_sig := (spearman_p <= 0.05 & spearman_r > 0)]

# ================= 2. 核心绘图函数 =================
plot_single_tool <- function(tool_name, data, fixed_plateau_list, peak_anno_whole, peak_anno_all, cand_label) {
  
  tool_data <- data[tool == tool_name]
  tool_peak_all <- peak_anno_all[tool == tool_name]
  
  plateau_min <- fixed_plateau_list[[tool_name]][1]
  plateau_max <- fixed_plateau_list[[tool_name]][2]
  
  cand_peak_u <- peak_anno_whole[tool == tool_name & WL_Label == cand_label, Peak_U]
  cand_peak_r <- peak_anno_whole[tool == tool_name & WL_Label == cand_label, Peak_r]
  ref_peak_u <- peak_anno_whole[tool == tool_name & WL_Label == "W80_L40", Peak_U]
  ref_peak_r <- peak_anno_whole[tool == tool_name & WL_Label == "W80_L40", Peak_r]
  
  sub_label <- sprintf(
    "Defined Plateau (whole): [%d, %d]  |  Peak (whole): %s = %d (r=%.2f), W80_L40 = %d (r=%.2f)",
    plateau_min, plateau_max, cand_label, cand_peak_u, cand_peak_r, ref_peak_u, ref_peak_r
  )
  
  shade_df <- data.frame(xmin = plateau_min, xmax = plateau_max, ymin = -Inf, ymax = Inf)
  
  p <- ggplot(tool_data, aes(x = U, y = spearman_r, color = WL_Label, linetype = WL_Label)) +
    geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#3498DB", alpha = 0.12, color = NA, inherit.aes = FALSE) +
    
    geom_line(linewidth = 1) +
    # 将显著性映射到形状 (shape) 上
    geom_point(aes(shape = is_sig), size = 2.5, stroke = 1.2, alpha = 0.8) +
    
    geom_vline(data = tool_peak_all, aes(xintercept = Peak_U, color = WL_Label), 
               linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    
    # 动态在图中每一个 Peak 点上方标注出对应的 r 值
    geom_text(data = tool_peak_all, 
              aes(x = Peak_U, y = Peak_r, label = sprintf("%.2f", Peak_r), color = WL_Label),
              size = 3.5, fontface = "bold", vjust = -0.8, show.legend = FALSE, inherit.aes = FALSE) +
    
    facet_wrap(~ region, scales = "free_y", ncol = 3) +
    
    scale_color_manual(values = setNames(c("#E74C3C", "#7F8C8D"), c(cand_label, "W80_L40"))) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), c(cand_label, "W80_L40"))) +
    # TRUE 为实心圆 (19)，FALSE 为空心圆 (1)
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                       labels = c("TRUE" = "Sig. (p < 0.05 & r > 0)", "FALSE" = "Not Sig. (p > 0.05 or r < 0)"),
                       name = "Significance") +
    
    # 增加顶部的留白比例 (0.15)，防止 geom_text 的数值越界被裁剪
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    
    theme_bw(base_size = 11) +
    labs(
      title = sprintf("Tool: %s", tool_name),
      subtitle = sub_label,
      x = "Target Length (U)",
      y = "Spearman Correlation (r)",
      color = "Parameters", linetype = "Parameters"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "#2C3E50", face = "italic"),
      strip.background = element_rect(fill = "#F0F3F4"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


# ================= 3. 批量循环处理与拼图 =================
log_msg("Step 3: Processing candidate lists and generating merged PDFs...")

for (cand in CANDIDATES) {
  W_cand <- cand["W"]
  L_cand <- cand["L"]
  cand_label <- sprintf("W%d_L%d", W_cand, L_cand)
  
  log_msg(sprintf("--> Processing %s vs W80_L40 ...", cand_label))
  
  target_dt <- all_dt[
    direction == "3prime" & 
    region %in% REGIONS_TO_EXTRACT & 
    ((W == W_cand & L == L_cand) | (W == 80 & L == 40))
  ]
  
  target_dt[, WL_Label := sprintf("W%d_L%d", W, L)]
  target_dt$WL_Label <- factor(target_dt$WL_Label, levels = c(cand_label, "W80_L40"))
  target_dt$region <- factor(target_dt$region, levels = c("whole", "CDS", "3'UTR"))
  
  # 同时提取 Peak U 与对应的 Max r (Peak_r)
  peak_whole <- target_dt[region == "whole", .(
    Peak_U = U[which.max(spearman_r)],
    Peak_r = max(spearman_r, na.rm = TRUE)
  ), by = .(tool, WL_Label)]
  
  peak_all_regions <- target_dt[, .(
    Peak_U = U[which.max(spearman_r)],
    Peak_r = max(spearman_r, na.rm = TRUE)
  ), by = .(tool, region, WL_Label)]
  
  plot_list <- list()
  for (t in TOOLS) {
    plot_list[[t]] <- plot_single_tool(t, target_dt, U_PLATEAU, peak_whole, peak_all_regions, cand_label)
  }
  
  merged_plot <- ggarrange(
    plotlist = plot_list, 
    ncol = 1, 
    nrow = length(TOOLS),
    common.legend = TRUE, 
    legend = "bottom",
    align = "v"
  )
  
  out_pdf <- file.path(OUT_DIR, sprintf("W%dL%d_vs_W80L40.pdf", W_cand, L_cand))
  ggsave(out_pdf, plot = merged_plot, width = 8, height = 2.5 * length(TOOLS))  
  log_msg(sprintf("    Saved merged PDF: %s", out_pdf))
}

log_msg("Module partD.01 Batch Output completed successfully.")