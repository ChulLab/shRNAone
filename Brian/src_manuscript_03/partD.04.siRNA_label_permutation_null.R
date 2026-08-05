#!/usr/bin/env Rscript
# partD.04.siRNA_label_permutation_null.R
#
# 最终参数的独立验证集质控：
# 1. 跨模块读取 partC 产出的 siRNA_raw_matrix 和 binary 分类标签。
# 2. 针对本模块 (Part D) 最终敲定的王牌参数（W35_L20_U7 / U11），提取其预测值。
# 3. 对生物学真实标签进行 1000 次全局洗牌 (Label Shuffling)。
# 4. 构建 Spearman r 和 AUC 的 Null 分布

library(data.table)
library(ggplot2)
library(ggpubr)
library(pROC)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
CLEAN_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/03.GT_siRNA_acc/00.clean_siRNA")
MATRIX_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/03.GT_siRNA_acc/02.siRNA_matrix")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partD/04.siRNA_null")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERMUTATIONS <- 1000

# Part D 最终敲定的参数组合
# TARGET_TESTS <- list(
#   list(tool = "shRNAone",  param = "W35_L20_U7",  col_name = "W_35_L_20_U_7"),
#   list(tool = "PspCas13b", param = "W35_L20_U11", col_name = "W_35_L_20_U_11"),
#   list(tool = "CasRx",     param = "W35_L20_U11", col_name = "W_35_L_20_U_11")
# )
TARGET_TESTS <- list(
  list(tool = "shRNAone",  param = "W35_L20_U7",  col_name = "W_35_L_20_U_7"),
  list(tool = "shRNAone",  param = "W80_L40_U7",  col_name = "W_80_L_40_U_7"),
  list(tool = "shRNAone",  param = "W80_L40_U6",  col_name = "W_80_L_40_U_6"),
  list(tool = "shRNAone",  param = "W80_L40_U16",  col_name = "W_80_L_40_U_16")
)
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 1. 加载数据 =================
log_msg("Step 1: Loading raw matrix and binary tags from upstream...")

raw_dt <- fread(file.path(MATRIX_DIR, "00.siRNA_raw_matrix.csv"))
bin_info <- fread(file.path(CLEAN_DIR, "siRNA_binary.csv"), select = c("ID", "group"))

# ================= 2. 核心置换引擎 =================
run_label_permutation <- function(raw_dt, bin_info, pred_col, tool_name, param_name) {
  log_msg(sprintf("  -> Processing %s (%s)", tool_name, param_name))
  
  df_cont <- raw_dt[!is.na(get(pred_col)) & !is.na(i_score), c("ID", "i_score", pred_col), with = FALSE]
  df_bin <- merge(df_cont, bin_info, by = "ID")
  df_bin[, is_func := as.numeric(group == "functional")]
  
  if (nrow(df_cont) < 20 || nrow(df_bin) < 20) {
    log_msg("    [Error] Not enough data points.")
    return(NULL)
  }
  
  true_r <- cor(df_cont[[pred_col]], df_cont$i_score, method = "spearman")
  true_auc <- as.numeric(roc(df_bin$is_func, df_bin[[pred_col]], quiet = TRUE)$auc)
  if (true_auc < 0.5) true_auc <- 1 - true_auc
  
  null_r_dist <- numeric(N_PERMUTATIONS)
  null_auc_dist <- numeric(N_PERMUTATIONS)
  
  for (i in seq_len(N_PERMUTATIONS)) {
    shuffled_i_score <- sample(df_cont$i_score)
    null_r_dist[i] <- cor(df_cont[[pred_col]], shuffled_i_score, method = "spearman")
    
    shuffled_func <- sample(df_bin$is_func)
    curr_auc <- as.numeric(roc(shuffled_func, df_bin[[pred_col]], quiet = TRUE)$auc)
    if (curr_auc < 0.5) curr_auc <- 1 - curr_auc
    null_auc_dist[i] <- curr_auc
  }
  
  emp_p_r <- (sum(abs(null_r_dist) >= abs(true_r)) + 1) / (N_PERMUTATIONS + 1)
  emp_p_auc <- (sum(null_auc_dist >= true_auc) + 1) / (N_PERMUTATIONS + 1)
  
  title_prefix <- sprintf("%s (%s)", tool_name, param_name)
  
  p_r <- ggplot(data.table(Null_r = null_r_dist), aes(x = Null_r)) +
    geom_density(fill = "#4DAF4A", color = "black", alpha = 0.6) +
    geom_vline(xintercept = true_r, color = "#E41A1C", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = true_r, y = Inf, 
             label = sprintf("True r: %.3f\nEmp. P: %.3f", true_r, emp_p_r), 
             vjust = 1.5, hjust = ifelse(true_r > 0, -0.1, 1.1), 
             color = "#E41A1C", fontface = "bold", size = 3) +
    theme_bw(base_size = 9) +
    labs(title = paste0(title_prefix, "\nSpearman r Null Dist."), x = "Permuted r", y = "Density") +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  p_auc <- ggplot(data.table(Null_auc = null_auc_dist), aes(x = Null_auc)) +
    geom_density(fill = "#377EB8", color = "black", alpha = 0.6) +
    geom_vline(xintercept = true_auc, color = "#E41A1C", linewidth = 1.2, linetype = "dashed") +
    annotate("text", x = true_auc, y = Inf, 
             label = sprintf("True AUC: %.3f\nEmp. P: %.3f", true_auc, emp_p_auc), 
             vjust = 1.5, hjust = -0.1, 
             color = "#E41A1C", fontface = "bold", size = 3) +
    theme_bw(base_size = 9) +
    labs(title = paste0(title_prefix, "\nAUC Null Dist."), x = "Permuted AUC", y = "Density") +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  stats_dt <- data.table(
    Tool = tool_name, Parameter = param_name,
    True_r = true_r, Emp_P_r = emp_p_r,
    True_AUC = true_auc, Emp_P_AUC = emp_p_auc
  )
  
  return(list(plot_r = p_r, plot_auc = p_auc, stats = stats_dt))
}

# ================= 3. 批量执行 =================
log_msg("Step 2: Starting permutation tests...")

all_stats <- list()
all_plots <- list()

for (test in TARGET_TESTS) {
  res <- run_label_permutation(raw_dt, bin_info, test$col_name, test$tool, test$param)
  if (!is.null(res)) {
    
    unique_key <- paste0(test$tool, "_", test$param)
    
    all_stats[[unique_key]] <- res$stats
    all_plots[[paste0(unique_key, "_r")]] <- res$plot_r
    all_plots[[paste0(unique_key, "_auc")]] <- res$plot_auc
  }
}

# ================= 4. 输出保存 =================
log_msg("Step 3: Saving results and rendering plots...")
if (length(all_plots) > 0) {
  final_stats <- rbindlist(all_stats)
  fwrite(final_stats, file.path(OUT_DIR, "01.siRNA_Permutation_Stats.csv"))
  
  merged_pdf <- ggarrange(plotlist = all_plots, ncol = 2, nrow = length(TARGET_TESTS))
  ggsave(file.path(OUT_DIR, "02.siRNA_Permutation_Distributions.pdf"), merged_pdf, 
         width = 7.5, height = 3 * length(TARGET_TESTS))
  log_msg("Module partD.04 finished perfectly.")
}