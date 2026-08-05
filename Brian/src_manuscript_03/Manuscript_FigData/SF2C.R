#!/usr/bin/env Rscript
# derived from: partD.03.spatial_permutation_null.R
#
# SF2C: Spatial Permutation Test for TTR and PCSK9 Accessibility vs log2FC Correlation

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
TTR_MAT <- file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix")
PCSK9_MAT <- file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix")

OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERMUTATIONS <- 1000

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

circular_shift <- function(x, shift_val) {
  n <- length(x)
  if (shift_val == 0 || shift_val == n) return(x)
  return(c(x[(shift_val + 1):n], x[1:shift_val]))
}

# ================= 1. 空间置换核心引擎 =================
run_spatial <- function(file_dir, gene_name) {
  files <- list.files(file_dir, pattern = "\\.csv$", full.names = TRUE)
  stat_list <- list()
  null_list <- list()
  
  for (f in files) {
    name_prefix <- gsub("\\.csv$", "", basename(f))
    
    # 正则清洗：彻底剔除基因名、_based 以及 _RNAxs，最后压缩 W_L_U
    clean_label <- name_prefix
    clean_label <- gsub("TTR_", "", clean_label)
    clean_label <- gsub("PCSK9_", "", clean_label)
    clean_label <- gsub("_based", "", clean_label)
    clean_label <- gsub("_RNAxs", "", clean_label) 
    clean_label <- gsub("_W([0-9]+)_L([0-9]+)_U([0-9]+)", "_W\\1L\\2U\\3", clean_label)
    
    dt <- fread(f)[!is.na(Accessibility) & !is.na(log2FC)]
    n_rows <- nrow(dt)
    if (n_rows < 20) next
    
    true_r <- cor(dt$Accessibility, dt$log2FC, method = "spearman")
    null_dist <- numeric(N_PERMUTATIONS)
    shift_steps <- sample(1:(n_rows - 1), N_PERMUTATIONS, replace = TRUE)
    
    for (i in seq_len(N_PERMUTATIONS)) {
      shifted_log2fc <- circular_shift(dt$log2FC, shift_steps[i])
      null_dist[i] <- cor(dt$Accessibility, shifted_log2fc, method = "spearman")
    }
    
    emp_p <- (sum(abs(null_dist) >= abs(true_r)) + 1) / (N_PERMUTATIONS + 1)
    
    stat_list[[clean_label]] <- data.table(
      Gene = gene_name, 
      Model_Label = clean_label, 
      True_Value = true_r, 
      Empirical_P = emp_p
    )
    
    null_list[[clean_label]] <- data.table(
      Gene = gene_name, 
      Model_Label = clean_label, 
      Null_Value = null_dist
    )
  }
  
  return(list(stat = rbindlist(stat_list), null = rbindlist(null_list)))
}

# ================= 2. 批量运行并合并数据 =================
log_msg("Processing TTR Spatial Permutation...")
res_ttr <- run_spatial(TTR_MAT, "TTR")

log_msg("Processing PCSK9 Spatial Permutation...")
res_pcsk9 <- run_spatial(PCSK9_MAT, "PCSK9")

stat_all <- rbindlist(list(res_ttr$stat, res_pcsk9$stat))
null_all <- rbindlist(list(res_ttr$null, res_pcsk9$null))

# 排序逻辑：按 Gene 和 真实 r 值排序
stat_all <- stat_all[order(Gene, True_Value)]
ordered_models <- unique(stat_all$Model_Label)

stat_all[, Model_Label := factor(Model_Label, levels = ordered_models)]
null_all[, Model_Label := factor(Model_Label, levels = ordered_models)]
stat_all[, Gene := factor(Gene, levels = c("TTR", "PCSK9"))]
null_all[, Gene := factor(Gene, levels = c("TTR", "PCSK9"))]

stat_all[, P_Label := sprintf("P = %.3f", Empirical_P)]
stat_all[Empirical_P < 0.001, P_Label := "P < 0.001"]

# ================= 3. 绘制主图 (小提琴图) =================
log_msg("Rendering primary violin plots...")

p_violin <- ggplot() +
  geom_violin(data = null_all, aes(x = Model_Label, y = Null_Value), 
              fill = "gray85", color = "gray60", linewidth = 0.2, scale = "width", alpha = 0.8) +
  
  geom_point(data = stat_all, aes(x = Model_Label, y = True_Value), 
             color = "#E74C3C", size = 0.8, stroke = 0.4) +
  
  geom_text(data = stat_all, aes(x = Model_Label, y = True_Value, label = P_Label), 
            color = "#E74C3C", fontface = "bold", size = 1.6, vjust = -0.8) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  
  facet_grid(. ~ Gene, scales = "free_x", space = "free_x") +
  
  labs(x = NULL, y = "Spearman (r)") +
  
  theme_bw(base_size = 5, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 5.5, margin = margin(t = 2, b = 2)),
    axis.title.y = element_text(face = "bold", size = 5.5, margin = margin(r = 2)),
    axis.text.x = element_text(face = "bold", color = "black", size = 4.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black", size = 5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    panel.spacing.x = unit(0.2, "lines"), 
    plot.margin = margin(2, 2, 1, 2)
  )

# ================= 4. 构建右侧详细工作流图 =================
log_msg("Constructing detailed textual workflow diagram...")

# 定义更详细的文本框内容与坐标
flow_dt <- data.table(
  x = 0.5,
  y = c(4.2, 3.1, 2.0, 0.9),
  label = c("1. Extract True Pairs\n(Target Acc. vs log2FC)",
            "2. Circular Shift log2FC\n(Break spatial synchrony)",
            "3. Build Null Distribution\n(1000x spatial permutations)",
            "4. Calculate Empirical P\n(True r vs Null background)")
)

# 动态调整箭头的坐标，避免与文本框重叠
arrow_dt <- data.table(
  x = 0.5, xend = 0.5,
  y = c(3.7, 2.6, 1.5), yend = c(3.55, 2.45, 1.35)
)

p_flow <- ggplot() +
  geom_label(data = flow_dt, aes(x = x, y = y, label = label),
             size = 1.5, fontface = "bold", lineheight = 1.2, # 稍微调小字号，增加行高让文字更易读
             fill = "gray98", color = "gray20", 
             label.size = 0.2, label.padding = unit(0.2, "lines")) +
  
  geom_segment(data = arrow_dt, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.04, "inches"), type = "closed"),
               linewidth = 0.3, color = "gray20") +
  
  coord_cartesian(xlim = c(0, 1), ylim = c(0.4, 4.6)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 0))

# ================= 5. 合并并输出极限尺寸成图 =================
log_msg("Merging plot and workflow into final layout...")

# 微微加大右侧宽度比例（从 0.3 提升到 0.35），适配加长的文本
p_combined <- ggarrange(
  p_violin, p_flow, 
  ncol = 2, nrow = 1, 
  widths = c(1, 0.35), 
  align = "h"
)

out_pdf <- file.path(OUT_DIR, "Supple_Spatial_Permutation_TTR_PCSK9_Combined.pdf")
ggsave(out_pdf, plot = p_combined, width = 4, height = 1.7, device = cairo_pdf)

# --- Save permutation stats ---
fwrite(stat_all, file.path(OUT_DIR, "SF2C_Spatial_Permutation_Stats.csv"))
log_msg(sprintf("Detailed workflow layout saved successfully to: %s", out_pdf))