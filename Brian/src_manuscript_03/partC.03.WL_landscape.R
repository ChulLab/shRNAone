# partC.03.WL_landscape.R 
# 
# 核心逻辑：
# 1. 采用滑动窗口算法全面扫描 5x5 区块。

library(ggplot2)
library(data.table)
library(dplyr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
WORK_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/01.Self_WL_list")

REGIONS <- c("whole", "3UTR")
METRICS <- c("Weighted_Score_Mean", "Weighted_Score_Max")

REF_W <- 80
REF_L <- 40
BLOCK_SIZE <- 5

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

heatmap_colors <- c("#313695", "#74ADD1", "#E0F3F8", "#FDAE61", "#F46D43", "#A50026")

# ================= 1. 绘图引擎函数 =================
generate_landscape_plot <- function(data, metric, top_point, top_block, title_prefix, zoom = FALSE) {
  
  # 动态判断当前提取的是 Mean 还是 Max 的各项工具独立跑分
  val_prefix <- ifelse(grepl("Mean", metric), "Mean_r", "Max_r")
  col_psp <- paste0(val_prefix, "_PspCas13b")
  col_shr <- paste0(val_prefix, "_shRNAone")
  col_cas <- paste0(val_prefix, "_CasRx_day5")
  
  # 提取 Ref 点的原始数据
  ref_row <- data[W == REF_W & L == REF_L]
  if(nrow(ref_row) > 0) {
    ref_txt <- sprintf("× Ref (%d,%d)  |  Psp: %.2f  shRNA: %.2f  CasRx: %.2f", 
                       REF_W, REF_L, ref_row[[col_psp]], ref_row[[col_shr]], ref_row[[col_cas]])
  } else {
    ref_txt <- sprintf("× Ref (%d,%d)  |  Data missing", REF_W, REF_L)
  }
  
  # 提取 Peak 标签
  peak_txt <- sprintf("+ Peak (%d,%d)  |  Psp: %.2f  shRNA: %.2f  CasRx: %.2f", 
                      top_point$W, top_point$L,
                      top_point[[col_psp]], top_point[[col_shr]], top_point[[col_cas]])
  
  # 提取 Block 标签
  block_txt <- sprintf("[ ] Top %dx%d Block (W:%d-%d, L:%d-%d)  |  Score: %.3f", 
                       BLOCK_SIZE, BLOCK_SIZE,
                       top_block$Block_W, top_block$Block_W_end,
                       top_block$Block_L, top_block$Block_L_end,
                       top_block$Block_Score)
  
  # 组合为左上角的信息面板文本 (换行排布)
  dashboard_txt <- paste(ref_txt, peak_txt, block_txt, sep = "\n\n")
  
  # 动态设定面板锚点 (左上角)
  anchor_x <- ifelse(zoom, 2, 5)
  anchor_y <- ifelse(zoom, 98, max(data$L) * 0.98)
  
  p <- ggplot(data, aes(x = W, y = L)) +
    geom_raster(aes_string(fill = metric), interpolate = TRUE) +
    geom_contour(aes_string(z = metric), color = "white", alpha = 0.4, bins = 15) +
    
    # 1. 细化的 Block 虚线框 (linewidth 降低，更精致)
    geom_rect(data = top_block,
              aes(xmin = Block_W - 0.5, xmax = Block_W_end + 0.5, 
                  ymin = Block_L - 0.5, ymax = Block_L_end + 0.5),
              fill = NA, color = "black", linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
    
    # 2. 细化的 Peak 最佳点 (+)
    geom_point(data = top_point, aes(x = W, y = L), 
               color = "black", shape = 3, size = 2, stroke = 0.8, inherit.aes = FALSE) +
    
    # 3. 细化的 Ref 经典参考点 (×)
    geom_point(data = data.frame(W = REF_W, L = REF_L), aes(x = W, y = L), 
               color = "black", shape = 4, size = 2, stroke = 0.8, inherit.aes = FALSE) +
    
    # 4. 左上角统一信息面板 (Dashboard)
    geom_text(data = data.frame(x = anchor_x, y = anchor_y, label = dashboard_txt),
              aes(x = x, y = y, label = label),
              color = "black", hjust = 0, vjust = 1, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
    
    scale_fill_gradientn(colors = heatmap_colors, name = metric) +
    labs(
      title = title_prefix,
      x = "Macro Window Size (W)",
      y = "Max Folding Span (L)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      # 取消坐标轴扩张，使其完全贴合
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  if (zoom) {
    p <- p + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100))
  } else {
    p <- p + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }
  
  return(p)
}

# ================= 2. 核心分析循环 =================
for (tgt_region in REGIONS) {
  file_name <- sprintf("00.WL_Performance_Matrix_%s.csv", tgt_region)
  matrix_dt <- fread(file.path(WORK_DIR, file_name))
  plot_data <- matrix_dt[W >= L]
  
  for (tgt_metric in METRICS) {
    log_msg(sprintf("Processing %s for %s region...", tgt_metric, tgt_region))
    
    # 1. 寻找单点绝对 Peak
    top_point <- plot_data[order(-get(tgt_metric))][1]
    
    # 2. 严苛的滑动窗口区块分析
    min_w <- min(plot_data$W)
    max_w <- max(plot_data$W) - (BLOCK_SIZE - 1)
    min_l <- min(plot_data$L)
    max_l <- max(plot_data$L) - (BLOCK_SIZE - 1)
    
    cand_grid <- as.data.table(expand.grid(Block_W = min_w:max_w, Block_L = min_l:max_l))
    cand_grid <- cand_grid[Block_W >= Block_L + (BLOCK_SIZE - 1)]
    cand_grid[, Block_W_end := Block_W + (BLOCK_SIZE - 1)]
    cand_grid[, Block_L_end := Block_L + (BLOCK_SIZE - 1)]
    
    block_stats <- plot_data[cand_grid, 
                             on = .(W >= Block_W, W <= Block_W_end, L >= Block_L, L <= Block_L_end),
                             .(Block_Score = mean(get(tgt_metric), na.rm = TRUE), Point_Count = .N), 
                             by = .EACHI]
    
    setnames(block_stats, 1:6, c("Block_W", "Block_W_end", "Block_L", "Block_L_end", "Block_Score", "Point_Count"))
    block_stats <- block_stats[Point_Count == (BLOCK_SIZE * BLOCK_SIZE)]
    setorder(block_stats, -Block_Score)
    
    block_out_file <- file.path(WORK_DIR, sprintf("02.Block_Stats_Strict_%s_%s.csv", tgt_region, tgt_metric))
    fwrite(block_stats, block_out_file)
    top_block <- block_stats[1]
    
    log_msg(sprintf("Top Strict Block found at W:%d-%d, L:%d-%d (Score: %.4f)", 
                    top_block$Block_W, top_block$Block_W_end,
                    top_block$Block_L, top_block$Block_L_end,
                    top_block$Block_Score))
    
    # 3. 绘制并输出视图
    p_full <- generate_landscape_plot(plot_data, tgt_metric, top_point, top_block, 
                                      sprintf("Landscape (%s): %s", tgt_metric, tgt_region), zoom = FALSE)
    pdf_full <- file.path(WORK_DIR, sprintf("01.WL_Landscape_Full_%s_%s.pdf", tgt_region, tgt_metric))
    pdf(pdf_full, width = 7.5, height = 5)
    print(p_full)
    dev.off()
    
    p_zoom <- generate_landscape_plot(plot_data, tgt_metric, top_point, top_block, 
                                      sprintf("Landscape Zoomed (W<=100): %s", tgt_region), zoom = TRUE)
    pdf_zoom <- file.path(WORK_DIR, sprintf("01.WL_Landscape_Zoom_%s_%s.pdf", tgt_region, tgt_metric))
    pdf(pdf_zoom, width = 6.5, height = 5)
    print(p_zoom)
    dev.off()
  }
}

log_msg("Module partC.03 visualization with aesthetic dashboard finished successfully.")