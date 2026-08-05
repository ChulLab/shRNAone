# 20.self_other_compare.R — Performance Comparison across Dry Models by Region.
#
# For each wet tool (CasRx_day5, shRNAone, PspCas13b):
#   1. Reads best_offset from 02/03/04.compare_wet_dry_...csv
#   2. Aligned dry models to wet data using best_offset (aligned by end_pos).
#   3. Evaluates Spearman correlation across Whole, CDS, and 3'UTR regions.
#   4. Saves full evaluation stats to CSV and outputs a faceted 3-row bar plot to PDF.
#
# conda env: self_model

library(ggplot2)
library(grid)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush=TRUE)
}

# ================= Configuration =================
PROJECT <- "/data/cai801/data/HKUcas"
WET_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
MODEL_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
FIGURE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

COR_METHOD <- "spearman"

TOOL_PREFIXES <- list(
  "CasRx_day5" = "02",
  "shRNAone"   = "03",
  "PspCas13b"  = "04"
)

# 包含 Whole, CDS, 3'UTR 三个部分
REGIONS <- c("Whole", "CDS", "3'UTR")

get_sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

calc_perf <- function(df, score_col, target_col, method_name, tool_name, region_name) {
  # 根据 region 筛选数据
  if (region_name == "Whole") {
    sub_df <- df
  } else {
    sub_df <- df[df$region == region_name, ]
  }
  
  sub_df <- sub_df[!is.na(sub_df[[score_col]]) & !is.na(sub_df[[target_col]]), ]
  if (nrow(sub_df) < 3) {
    return(data.frame(Tool = tool_name, Region = region_name, Method = method_name, 
                      Correlation = NA, P_value = NA, Sig_Label = "ns", N = nrow(sub_df)))
  }
  
  test <- cor.test(sub_df[[score_col]], sub_df[[target_col]], method = COR_METHOD, exact = FALSE)
  data.frame(
    Tool = tool_name, Region = region_name, Method = method_name,
    Correlation = test$estimate, P_value = test$p.value,
    Sig_Label = get_sig_label(test$p.value), N = nrow(sub_df),
    stringsAsFactors = FALSE
  )
}

log_msg("20. Comparing Dry-models performance across Whole, CDS, and 3'UTR")
all_perf <- data.frame()

for (tool in names(TOOL_PREFIXES)) {
  prefix <- TOOL_PREFIXES[[tool]]
  
  # 读取清洗后的湿实验数据 (TTR_[tool]_clean.csv)
  wet_file <- file.path(WET_DIR, paste0("TTR_", tool, "_clean.csv"))
  if (!file.exists(wet_file)) {
    log_msg(sprintf("  [SKIP] Wet file not found: %s", wet_file))
    next
  }
  df_wet <- read.csv(wet_file)
  
  # 读取 offset offset 比对汇总表
  stats_file <- file.path(FIGURE_DIR, paste0(prefix, ".compare_wet_dry_", tool, ".csv"))
  if (!file.exists(stats_file)) {
    log_msg(sprintf("  [SKIP] Stats file not found: %s", stats_file))
    next
  }
  df_stats <- read.csv(stats_file)
  
  for (i in 1:nrow(df_stats)) {
    dry_model <- df_stats$dry_model[i]
    best_off <- df_stats$best_offset[i]
    dry_file <- file.path(MODEL_DIR, paste0(dry_model, "_aligned.csv"))
    if (!file.exists(dry_file)) next
    
    df_dry <- read.csv(dry_file)
    df_dry_agg <- aggregate(prediction_score ~ end_pos, data = df_dry, FUN = mean)
    df_dry_agg$target_wet_ep <- df_dry_agg$end_pos - best_off
    
    # 将最佳 Offset 处的预测得分与湿实验合并
    df_merged <- merge(df_wet, df_dry_agg, by.x = "end_pos", by.y = "target_wet_ep")
    
    # 分别计算 Whole, CDS, 3'UTR 的相关性
    for (reg in REGIONS) {
      perf <- calc_perf(df_merged, "prediction_score", "mean_log2FC", dry_model, tool, reg)
      all_perf <- rbind(all_perf, perf)
    }
  }
}

# 1. 保存全部分组评估数据到 CSV
out_csv <- file.path(FIGURE_DIR, "05.self_other_compare.csv")
write.csv(all_perf, out_csv, row.names = FALSE)
log_msg(sprintf("  Saved performance stats: %s", basename(out_csv)))

# ================= 2. Plotting =================
all_perf <- all_perf[!is.na(all_perf$Correlation), ]
all_perf$Tool <- factor(all_perf$Tool, levels = names(TOOL_PREFIXES))
all_perf$Region <- factor(all_perf$Region, levels = REGIONS)

methods <- sort(unique(all_perf$Method))
all_perf$Method <- factor(all_perf$Method, levels = methods)

all_perf$label_y <- ifelse(all_perf$Correlation >= 0, 
                           all_perf$Correlation + 0.04, 
                           all_perf$Correlation - 0.04)

p_bar <- ggplot(all_perf, aes(x = Method, y = Correlation)) +
  geom_bar(stat = "identity", fill = "gray70", color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = Sig_Label, y = label_y), size = 2.0, 
            vjust = ifelse(all_perf$Correlation >= 0, 0, 1)) +
  # 行按 Region，列按 Tool 展开
  facet_grid(Region ~ Tool, scales = "free_x", space = "free_x") + 
  # 锁定 Y 轴刻度 0 - 0.6，避免显著性标注超界
  scale_y_continuous(limits = c(-0.1, 0.65), breaks = c(0, 0.25, 0.5)) +
  labs(
    title = sprintf("Dry Models Performance Comparison (%s r)", tools::toTitleCase(COR_METHOD)),
    y = sprintf("%s r", tools::toTitleCase(COR_METHOD)), x = NULL
  ) +
  theme_bw(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 3)),
    strip.text = element_text(size = 6.5, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.3),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 5.5, color = "black"),
    axis.title.y = element_text(size = 6.5, face = "bold", margin = margin(r = 2)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

# 多加了 2 行 Region 图，适度增加 PDF 高度 (5.5 x 4.2 吋)
pdf_path <- file.path(FIGURE_DIR, "05.self_other_compare.pdf")
pdf(pdf_path, width = 5, height = 2.5) 
print(p_bar)
dev.off()
log_msg(sprintf("  Saved multi-region bar plot: %s", basename(pdf_path)))

log_msg("20. Done.")