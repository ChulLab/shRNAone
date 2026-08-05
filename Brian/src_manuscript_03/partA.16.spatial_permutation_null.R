#!/usr/bin/env Rscript
# partA.16.spatial_permutation_null.R
#

library(data.table)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
WET_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
MODEL_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
STATS_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/05.spatial_null")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERMUTATIONS <- 1000
ALIGN_BY <- "end_pos" 

TOOL_PREFIXES <- list(
  "CasRx_day5" = "02",
  "shRNAone"   = "03",
  "PspCas13b"  = "04"
)

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
run_spatial_permutation <- function(df_merged, true_r, title_prefix) {
  n_rows <- nrow(df_merged)
  if (n_rows < 20) return(NULL)
  
  null_r_dist <- numeric(N_PERMUTATIONS)
  shift_steps <- sample(1:(n_rows - 1), N_PERMUTATIONS, replace = TRUE)
  
  # [Debug] 排查数据异常
  if (any(is.na(df_merged$prediction_score))) {
    log_msg("    [Warning] NA values found in prediction_score.")
  }
  if (var(df_merged$prediction_score, na.rm = TRUE) == 0) {
    log_msg("    [Warning] prediction_score has ZERO variance (all values identical).")
  }

  # 在最佳对齐位置上，保持干实验得分不变，环状移位湿实验 log2FC
  for (i in seq_len(N_PERMUTATIONS)) {
    shifted_log2fc <- circular_shift(df_merged$mean_log2FC, shift_steps[i])
    # 强制忽略 NA 值，防止单点 NA 毁掉整个向量
    null_r_dist[i] <- cor(df_merged$prediction_score, shifted_log2fc, method = "spearman", use = "complete.obs")
  }
  
  # 如果防线被零方差击穿，null_r_dist 全是 NA，需要手动给个假数据防止画图崩溃
  if (all(is.na(null_r_dist))) {
    log_msg("    [Error] All permuted r values are NA. Generating blank plot safeguard.")
    null_r_dist <- rep(0, N_PERMUTATIONS)
  }
  
  # 计算 Empirical P-value (双侧)
  emp_p <- (sum(abs(null_r_dist) >= abs(true_r)) + 1) / (N_PERMUTATIONS + 1)
  
  stat_res <- data.table(
    Model_Comparison = title_prefix,
    True_r = true_r,
    Empirical_P = emp_p,
    Null_Mean = mean(null_r_dist),
    Null_SD = sd(null_r_dist)
  )
  
  null_dt <- data.table(Null_r = null_r_dist)
  p <- ggplot(null_dt, aes(x = Null_r)) +
    geom_density(fill = "gray80", color = "black", alpha = 0.7) +
    geom_vline(xintercept = true_r, color = "#E41A1C", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = true_r, y = Inf, 
             label = sprintf("True r: %.3f\nEmp. P: %.3f", true_r, emp_p), 
             vjust = 1.5, hjust = ifelse(true_r > 0, -0.1, 1.1), 
             color = "#E41A1C", fontface = "bold", size = 3) +
    theme_bw(base_size = 9) +
    labs(
      title = title_prefix,
      x = "Spearman r (Null)", y = "Density"
    ) +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  return(list(stats = stat_res, plot = p))
}

# ================= 2. 批量处理 =================
log_msg("Starting Spatial Permutation for Part A External Models...")

all_stats <- list()
all_plots <- list()

for (tool in names(TOOL_PREFIXES)) {
  prefix <- TOOL_PREFIXES[[tool]]
  stats_file <- file.path(STATS_DIR, paste0(prefix, ".compare_wet_dry_", tool, ".csv"))
  wet_file <- file.path(WET_DIR, paste0("TTR_", tool, "_clean.csv"))
  
  if (!file.exists(stats_file) || !file.exists(wet_file)) next
  
  df_stats <- fread(stats_file)
  df_wet <- fread(wet_file)
  
  for (i in 1:nrow(df_stats)) {
    dry_model <- df_stats$dry_model[i]
    best_off <- df_stats$best_offset[i]
    true_r <- df_stats$cor_r[i]
    
    dry_file <- file.path(MODEL_DIR, paste0(dry_model, "_aligned.csv"))
    if (!file.exists(dry_file)) next
    
    df_dry <- fread(dry_file)
    df_dry_sub <- df_dry[, c(ALIGN_BY, "prediction_score"), with = FALSE]
    setnames(df_dry_sub, ALIGN_BY, "pos_align")
    df_dry_agg <- df_dry_sub[, .(prediction_score = mean(prediction_score)), by = pos_align]
    
    # 根据最佳 offset 对齐坐标
    df_dry_agg[, target_wet_pos := pos_align - best_off]
    
    # 严格的 Inner Join 保证序列无空缺
    df_merged <- merge(df_wet, df_dry_agg, by.x = ALIGN_BY, by.y = "target_wet_pos")
    df_merged <- df_merged[order(get(ALIGN_BY))]
    
    title_prefix <- sprintf("%s vs %s (Off:%d)", tool, dry_model, best_off)
    log_msg(sprintf("  -> Permuting %s", title_prefix))
    
    res <- run_spatial_permutation(df_merged, true_r, title_prefix)
    if (!is.null(res)) {
      all_stats[[title_prefix]] <- res$stats
      all_plots[[title_prefix]] <- res$plot
    }
  }
}

# ================= 3. 保存与拼图 =================
if (length(all_stats) > 0) {
  final_stats <- rbindlist(all_stats)
  fwrite(final_stats, file.path(OUT_DIR, "01.PartA_Spatial_Null_Statistics.csv"))
  
  merged_pdf <- ggarrange(plotlist = all_plots, ncol = 3, nrow = ceiling(length(all_plots) / 3))
  ggsave(file.path(OUT_DIR, "02.PartA_Spatial_Null_Distributions.pdf"), merged_pdf, 
         width = 10, height = 2.5 * ceiling(length(all_plots) / 3))
  
  log_msg("Module partA.16 Spatial Permutation finished perfectly.")
}