#!/usr/bin/env Rscript
# partD.03.spatial_permutation_null.R
#
# 空间自相关假阳性排查：
# 1. 遍历 partD.02 生成的 PCSK9 和 TTR 验证矩阵。
# 2. 对湿实验数据 (log2FC) 实施 1,000 次环状移位 (Circular Block Permutation)，打破点对点对应关系，但完美保留空间自相关性。
# 3. 计算 Null 分布的 Spearman r，并计算 Empirical P-value。
# 4. 生成高颜值的 Null Distribution 密度图，红线标记 True r。

library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
TTR_MAT <- file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix")
PCSK9_MAT <- file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partD/03.spatial_null")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERMUTATIONS <- 1000

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# 环状移位核心函数
circular_shift <- function(x, shift_val) {
  n <- length(x)
  if (shift_val == 0 || shift_val == n) return(x)
  return(c(x[(shift_val + 1):n], x[1:shift_val]))
}

# ================= 1. 核心检验引擎 =================
run_spatial_permutation <- function(matrix_file, title_prefix) {
  if (!file.exists(matrix_file)) return(NULL)
  
  dt <- fread(matrix_file)
  # 确保剔除 NA
  dt <- dt[!is.na(Accessibility) & !is.na(log2FC)]
  n_rows <- nrow(dt)
  
  if (n_rows < 20) {
    log_msg(sprintf("    [SKIP] Not enough targets in %s", basename(matrix_file)))
    return(NULL)
  }
  
  # 计算真实的 Spearman r
  true_test <- cor.test(dt$Accessibility, dt$log2FC, method = "spearman", exact = FALSE)
  true_r <- true_test$estimate
  
  # 运行 1000 次环状移位 Permutation
  null_r_dist <- numeric(N_PERMUTATIONS)
  # 随机生成移位步长 (1 到 n_rows-1)
  shift_steps <- sample(1:(n_rows - 1), N_PERMUTATIONS, replace = TRUE)
  
  for (i in seq_len(N_PERMUTATIONS)) {
    # 仅移位 log2FC 向量，保留 Accessibility 不变
    shifted_log2fc <- circular_shift(dt$log2FC, shift_steps[i])
    null_r_dist[i] <- cor(dt$Accessibility, shifted_log2fc, method = "spearman")
  }
  
  # 计算 Empirical P-value (双侧检验)
  emp_p <- (sum(abs(null_r_dist) >= abs(true_r)) + 1) / (N_PERMUTATIONS + 1)
  
  # 汇总数据
  stat_res <- data.table(
    Dataset = title_prefix,
    True_r = true_r,
    Parametric_P = true_test$p.value,
    Empirical_P = emp_p,
    Null_Mean = mean(null_r_dist),
    Null_SD = sd(null_r_dist)
  )
  
  # 绘图：Null Distribution 密度图
  null_dt <- data.table(Null_r = null_r_dist)
  p <- ggplot(null_dt, aes(x = Null_r)) +
    geom_density(fill = "gray80", color = "black", alpha = 0.7) +
    geom_vline(xintercept = true_r, color = "#E41A1C", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = true_r, y = Inf, 
             label = sprintf("True r: %.3f\nEmpirical P: %.3f", true_r, emp_p), 
             vjust = 1.5, hjust = ifelse(true_r > 0, -0.1, 1.1), 
             color = "#E41A1C", fontface = "bold", size = 3.5) +
    theme_bw(base_size = 10) +
    labs(
      title = title_prefix,
      subtitle = sprintf("Circular Block Permutation (N=%d)", N_PERMUTATIONS),
      x = "Spearman r (Null Distribution)",
      y = "Density"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(color = "#555555", size = 9)
    )
  
  return(list(stats = stat_res, plot = p))
}

# ================= 2. 批量处理与输出 =================
log_msg("Step 1: Running Spatial Permutation for TTR Models...")
ttr_files <- list.files(TTR_MAT, pattern = "\\.csv$", full.names = TRUE)
ttr_results <- list()
ttr_plots <- list()

for (f in ttr_files) {
  name_prefix <- gsub(".csv", "", basename(f))
  log_msg(sprintf("  -> Permuting %s", name_prefix))
  res <- run_spatial_permutation(f, name_prefix)
  if (!is.null(res)) {
    ttr_results[[name_prefix]] <- res$stats
    ttr_plots[[name_prefix]] <- res$plot
  }
}

log_msg("Step 2: Running Spatial Permutation for PCSK9 Models...")
pcsk9_files <- list.files(PCSK9_MAT, pattern = "\\.csv$", full.names = TRUE)
pcsk9_results <- list()
pcsk9_plots <- list()

for (f in pcsk9_files) {
  name_prefix <- gsub(".csv", "", basename(f))
  log_msg(sprintf("  -> Permuting %s", name_prefix))
  res <- run_spatial_permutation(f, name_prefix)
  if (!is.null(res)) {
    pcsk9_results[[name_prefix]] <- res$stats
    pcsk9_plots[[name_prefix]] <- res$plot
  }
}

# ================= 3. 保存统计表与拼图 =================
log_msg("Step 3: Saving results and generating PDFs...")

all_stats <- rbindlist(c(ttr_results, pcsk9_results))
fwrite(all_stats, file.path(OUT_DIR, "01.Spatial_Null_Statistics.csv"))

# 输出 TTR 的拼图 (按 3 列排布)
if (length(ttr_plots) > 0) {
  ttr_pdf <- ggarrange(plotlist = ttr_plots, ncol = 3, nrow = ceiling(length(ttr_plots) / 3))
  ggsave(file.path(OUT_DIR, "02.TTR_Spatial_Null_Distributions.pdf"), ttr_pdf, width = 12, height = 3 * ceiling(length(ttr_plots) / 3))
}

# 输出 PCSK9 的拼图
if (length(pcsk9_plots) > 0) {
  pcsk9_pdf <- ggarrange(plotlist = pcsk9_plots, ncol = 2, nrow = ceiling(length(pcsk9_plots) / 2))
  ggsave(file.path(OUT_DIR, "03.PCSK9_Spatial_Null_Distributions.pdf"), pcsk9_pdf, width = 8, height = 3 * ceiling(length(pcsk9_plots) / 2))
}

log_msg("Module partD.03 Spatial Permutation finished perfectly.")