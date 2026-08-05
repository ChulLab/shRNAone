# partC.01.define_refU_plateau.R 
# 
# 核心逻辑：锚定 W=80, L=40 作为统一基准。

library(ggplot2)
library(dplyr)
library(data.table)

# ================= 0. 配置与数据加载 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/00.Ref_U_plateau")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")
DIRECTIONS <- c("3prime", "downstream")

# ----------------- 核心传参设计 -----------------
# 模式 1: 自动计算 (auto)
PLATEAU_CONFIG <- "auto"
# 定义自动计算时的阈值乘数，例如 0.8 表示取 Peak r 值的 80% 以上作为平台区
PLATEAU_THRESHOLD <- 0.75 

# 模式 2: 手动人为指定区域 (直接解除下方注释即可覆盖 auto 模式)
# PLATEAU_CONFIG <- list(
#   shRNAone   = c(10, 18),
#   PspCas13b  = c(12, 20),
#   CasRx_day5 = c(8, 15)
# )
# ------------------------------------------------

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# 极致优化 I/O：使用 fread 仅读取绘图必须的列，节约内存并极速加载
log_msg("Loading and merging data using optimized data.table I/O...")
cols_to_keep <- c("W", "L", "U", "direction", "region", "spearman_r", "spearman_p")
all_data_list <- lapply(TOOLS, function(t) {
  fpath <- file.path(IN_DIR, paste0(t, "_grid.csv"))
  if (file.exists(fpath)) {
    # 极速读取，并自动打上 tool 标签
    dt <- fread(fpath, select = cols_to_keep)
    dt$tool <- t
    return(dt)
  }
})
ref_data <- rbindlist(all_data_list)

# ================= 1. 数据全局过滤与重构 =================
log_msg("Filtering data for Reference Baseline: W=80, L=40")

# 使用 data.table 语法或 dplyr 均可，这里接管给 dplyr 做后续管道操作
ref_data <- ref_data %>%
  filter(W == 80, L == 40, direction %in% DIRECTIONS) %>%
  mutate(
    region = factor(region, levels = c("whole", "CDS", "3'UTR")),
    tool = factor(tool, levels = TOOLS),
    Is_Significant = spearman_p < 0.05 & spearman_r > 0
  )

# ================= 2. 循环处理并独立输出每个 Direction =================
for (target_dir in DIRECTIONS) {
  log_msg(sprintf("Processing and plotting for direction: %s", target_dir))
  
  sub_data <- ref_data %>% filter(direction == target_dir)
  dir_color <- ifelse(target_dir == "3prime", "#E41A1C", "#984EA3")
  
  # ---------- 极速计算 Peak 和 Plateau 区间 (消灭 rowwise 瓶颈) ----------
  plateau_stats <- sub_data %>%
    group_by(tool, region) %>%
    summarise(
      Max_Cor = max(spearman_r, na.rm = TRUE),
      Peak_U  = U[which.max(spearman_r)],
      Peak_P  = spearman_p[which.max(spearman_r)],
      
      # 内联向量化计算 Plateau_Min
      Plateau_Min = if (max(spearman_r, na.rm = TRUE) > 0) {
        if (is.character(PLATEAU_CONFIG) && PLATEAU_CONFIG == "auto") {
          min(U[spearman_r >= PLATEAU_THRESHOLD * max(spearman_r, na.rm = TRUE) & spearman_r > 0], na.rm = TRUE)
        } else {
          PLATEAU_CONFIG[[as.character(tool[1])]][1]
        }
      } else { NA_real_ },
      
      # 内联向量化计算 Plateau_Max
      Plateau_Max = if (max(spearman_r, na.rm = TRUE) > 0) {
        if (is.character(PLATEAU_CONFIG) && PLATEAU_CONFIG == "auto") {
          max(U[spearman_r >= PLATEAU_THRESHOLD * max(spearman_r, na.rm = TRUE) & spearman_r > 0], na.rm = TRUE)
        } else {
          PLATEAU_CONFIG[[as.character(tool[1])]][2]
        }
      } else { NA_real_ },
      
      .groups = "drop"
    ) %>%
    mutate(
      Label_X = Peak_U,
      Label_Y = Max_Cor + 0.12
    )
  
  csv_name <- sprintf("01.Ref_U_Plateau_Stats_%s.csv", target_dir)
  write.csv(plateau_stats, file.path(OUT_DIR, csv_name), row.names = FALSE)
  
  # ---------- 高维可视化绘图 ----------
  y_max_limit <- max(plateau_stats$Label_Y, na.rm = TRUE) + 0.05
  y_min_limit <- min(sub_data$spearman_r, na.rm = TRUE) - 0.05
  
  p_traj <- ggplot(sub_data, aes(x = U, y = spearman_r)) +
    # 1. 绘制 Plateau 阴影背景
    geom_rect(data = plateau_stats %>% filter(!is.na(Plateau_Min)), 
              aes(xmin = Plateau_Min, xmax = Plateau_Max, ymin = -Inf, ymax = Inf), 
              fill = dir_color, alpha = 0.15, inherit.aes = FALSE, color = NA) +
    
    # 2. Y=0 基准线
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    
    # 3. 轨迹线与显著性点
    geom_line(color = dir_color, alpha = 0.8, linewidth = 1) +
    geom_point(aes(shape = Is_Significant), color = dir_color, size = 2, stroke = 1, fill = "white") +
    
    # 4. 动态箭头
    geom_segment(data = plateau_stats %>% filter(Max_Cor > 0),
                 aes(x = Label_X, y = Label_Y - 0.03, xend = Peak_U, yend = Max_Cor + 0.02),
                 color = dir_color, arrow = arrow(length = unit(0.12, "cm"), type = "closed"), linewidth = 0.5) +
    
    # 5. 信息浓缩标签
    geom_text(data = plateau_stats %>% filter(Max_Cor > 0),
              aes(x = Label_X, y = Label_Y, 
                  label = sprintf("Peak U=%d (r=%.2f)\nPlateau: %d-%d", 
                                  Peak_U, Max_Cor, Plateau_Min, Plateau_Max)),
              color = dir_color, size = 2.8, vjust = 0, lineheight = 1.1, fontface = "bold") +
    
    # 分面与样式
    facet_grid(region ~ tool) +
    scale_y_continuous(limits = c(y_min_limit, y_max_limit)) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1), 
                       labels = c("TRUE" = "Raw p < 0.05", "FALSE" = "p >= 0.05 / Negative")) +
    labs(
      title = sprintf("Reference U Plateau Trajectories: %s (W=80, L=40)", target_dir),
      subtitle = sprintf("Shaded areas indicate U Plateau (>%g%% Max). Solid points denote raw significance.", PLATEAU_THRESHOLD * 100),
      x = "Unpaired Region Length (U)",
      y = "Spearman Correlation (r)"
    ) +
    theme_bw(base_size = 8) + # 适配小尺寸画幅的字体基础大小
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  # 调整 PDF 输出参数以满足 Manuscript 排版标准
  pdf_name <- sprintf("02.Ref_U_Trajectories_%s.pdf", target_dir)
  pdf(file.path(OUT_DIR, pdf_name), width = 7, height = 6)
  print(p_traj)
  dev.off()
}

log_msg("Module partC.01 finished successfully in high-performance mode.")