#!/usr/bin/env Rscript
# derived from: partA.16.spatial_permutation_null.R
# 
# SF2F: Spatial Permutation Test for Wet vs Dry Correlation across all models

library(data.table)
library(dplyr)
library(ggplot2)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
WET_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
MODEL_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
STATS_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")

OUT_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/Manuscript_FigData/supple"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_PERMUTATIONS <- 1000
ALIGN_BY <- "end_pos" 

# 【修改 1】：剔除 shRNAone，仅保留 CasRx 和 Psp
TOOL_PREFIXES <- list(
  "CasRx_day5" = "02",
  "PspCas13b"  = "04"
)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

circular_shift <- function(x, shift_val) {
  n <- length(x)
  if (shift_val == 0 || shift_val == n) return(x)
  return(c(x[(shift_val + 1):n], x[1:shift_val]))
}

# ================= 1. 数据收集引擎 =================
log_msg("Starting Spatial Permutation & Data Collection...")

all_null_dt <- list()
all_stat_dt <- list()

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
    
    df_dry_agg[, target_wet_pos := pos_align - best_off]
    
    df_merged <- merge(df_wet, df_dry_agg, by.x = ALIGN_BY, by.y = "target_wet_pos")
    df_merged <- df_merged[order(get(ALIGN_BY))]
    
    n_rows <- nrow(df_merged)
    if (n_rows < 20) next
    
    null_r_dist <- numeric(N_PERMUTATIONS)
    shift_steps <- sample(1:(n_rows - 1), N_PERMUTATIONS, replace = TRUE)
    
    for (k in seq_len(N_PERMUTATIONS)) {
      shifted_log2fc <- circular_shift(df_merged$mean_log2FC, shift_steps[k])
      null_r_dist[k] <- cor(df_merged$prediction_score, shifted_log2fc, method = "spearman", use = "complete.obs")
    }
    
    if (all(is.na(null_r_dist))) null_r_dist <- rep(0, N_PERMUTATIONS)
    
    emp_p <- (sum(abs(null_r_dist) >= abs(true_r)) + 1) / (N_PERMUTATIONS + 1)
    
    model_label <- sprintf("%s (Off:%d)", dry_model, best_off)
    
    all_stat_dt[[model_label]] <- data.table(
      Tool = tool,
      Model_Label = model_label,
      True_r = true_r,
      Empirical_P = emp_p
    )
    
    all_null_dt[[model_label]] <- data.table(
      Tool = tool,
      Model_Label = model_label,
      Null_r = null_r_dist
    )
  }
}

df_stat_merged <- rbindlist(all_stat_dt)
df_null_merged <- rbindlist(all_null_dt)

# ================= 2. 数据排序与因子化 =================
df_stat_merged <- df_stat_merged[order(Tool, True_r)]
ordered_models <- unique(df_stat_merged$Model_Label)

df_stat_merged[, Model_Label := factor(Model_Label, levels = ordered_models)]
df_null_merged[, Model_Label := factor(Model_Label, levels = ordered_models)]
df_stat_merged[, Tool := factor(Tool, levels = names(TOOL_PREFIXES))]
df_null_merged[, Tool := factor(Tool, levels = names(TOOL_PREFIXES))]

df_stat_merged[, P_Label := sprintf("P = %.3f", Empirical_P)]
df_stat_merged[Empirical_P < 0.001, P_Label := "P < 0.001"]

# ================= 3. 高密度上下堆叠出图 (2.3 x 3) =================
log_msg("Rendering 2x1 Stacked Vertical Violin Layout...")

p_final <- ggplot() +
  geom_violin(data = df_null_merged, aes(x = Model_Label, y = Null_r), 
              fill = "gray85", color = "gray60", linewidth = 0.3, scale = "width", alpha = 0.8) +
  
  geom_point(data = df_stat_merged, aes(x = Model_Label, y = True_r), 
             color = "#E74C3C", size = 1.2, stroke = 0.5) +
  
  geom_text(data = df_stat_merged, aes(x = Model_Label, y = True_r, label = P_Label), 
            color = "#E74C3C", fontface = "bold", size = 2, vjust = -0.8) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # 【修改 2】：使用 ncol = 1 让分面上下堆叠
  facet_wrap(~ Tool, ncol = 1, scales = "free_x") +
  
  labs(
    # 为了适应更窄的宽度，稍微精简了标题并允许它居中
    title = NULL,
    x = NULL,
    y = "Spearman Correlation (r)"
  ) +
  
  theme_bw(base_size = 6, base_family = "Arial") +
  theme(
    # 调小了主标题字号以适应 2.3 英寸宽幅
    plot.title = element_text(face = "bold", size = 6.5, hjust = 0.5, margin = margin(b = 4), lineheight = 1.1),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 7, margin = margin(t = 2, b = 2)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 3), size = 7),
    axis.text.x = element_text(face = "bold", color = "black", size = 5.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black", size = 6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.4),
    panel.spacing.y = unit(0.4, "lines"), # 增加上下分面的间距
    plot.margin = margin(4, 4, 2, 2)
  )

out_pdf <- file.path(OUT_DIR, "Supple_Spatial_Permutation.pdf")

# 【修改 3】：严格按照要求的 2.3 x 3 尺寸保存
ggsave(out_pdf, plot = p_final, width = 1.6, height = 3, device = cairo_pdf)

# --- Save permutation stats ---
fwrite(df_stat_merged, file.path(OUT_DIR, "SF2F_External_Spatial_Permutation_Stats.csv"))
log_msg(sprintf("Stacked vertical layout saved successfully to: %s", out_pdf))