# partC.12.WL_compare.R
#
# 核心逻辑：
# 1. 载入 Continuous 和 Binary 评估结果 (过滤 U=6~16)。
# 2. 提取 AUC 和 Spearman 的第一四分位数 (Q1) 作为核心鲁棒性指标。
# 3. 使用稳健的 Z-score 平均综合 Block_Score / Spearman_Q1 / AUC_Q1（排除 Reference 对均值/方差的影响）。
# 4. 图表排版全部按 Composite_Score 降序排序。

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置与初始化 =================
BASE_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partC"
EVAL_DIR <- file.path(BASE_DIR, "03.GT_siRNA_acc", "02.siRNA_matrix")
CAND_FILE <- file.path(BASE_DIR, "01.Self_WL_list", "03.selected_WL_Candidates.csv")
OUT_DIR <- file.path(BASE_DIR, "03.GT_siRNA_acc", "03.WL_selection")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# 归一化函数 (Min-Max Scaling)
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# ================= 1. 数据载入与清洗 =================
log_msg("Step 1: Loading evaluation metrics...")

cont_df <- fread(file.path(EVAL_DIR, "01.Evaluation_Continuous.csv"))
bin_df  <- fread(file.path(EVAL_DIR, "02.Evaluation_Binary.csv"))
cands   <- fread(CAND_FILE)

# 限制 U 的范围 6-16
cont_df <- cont_df[U >= 6 & U <= 16]
bin_df  <- bin_df[U >= 6 & U <= 16]

eval_merged <- merge(
  cont_df[, .(W, L, U, Spearman_r, P_value)],
  bin_df[, .(W, L, U, AUC, Cliff_Delta)],
  by = c("W", "L", "U")
)

# 提取参考基准点 (W=80, L=40, U=8 和 U=16)
ref_U8  <- eval_merged[W == 80 & L == 40 & U == 8]
ref_U16 <- eval_merged[W == 80 & L == 40 & U == 16]

# ================= 2. 下限提取与稳健综合打分 =================
log_msg("Step 2: Aggregating using Q1 and computing robust Composite Score...")

# 计算第一四分位数 (Q1, 25%) 作为保底能力的评估指标
agg_df <- eval_merged[, .(
  Spearman_Q1 = quantile(Spearman_r, 0.25, na.rm = TRUE),
  AUC_Q1      = quantile(AUC, 0.25, na.rm = TRUE)
), by = .(W, L)]

final_df <- merge(cands, agg_df, by = c("W", "L"), all.x = TRUE)
final_df <- final_df[!is.na(AUC_Q1)]

# ---------- 给 Reference 填入真实的 Block_Score ----------
# 优先从 whole 区域性能矩阵读取 W=80, L=40 的 Weighted_Score_Mean
ref_block_file <- file.path(BASE_DIR, "01.Self_WL_list", "00.WL_Performance_Matrix_whole.csv")
if (file.exists(ref_block_file)) {
  ref_mat <- fread(ref_block_file)
  ref_score <- ref_mat[W == 80 & L == 40, Weighted_Score_Mean]
  if (length(ref_score) == 1 && !is.na(ref_score)) {
    final_df[Source == "Reference", Block_Score := as.numeric(ref_score)]
    log_msg(sprintf("Reference Block_Score filled from matrix: %.6f", ref_score))
  }
}

# 若仍为 NA，则使用候选中的最小值（不再人为制造极端 outlier）
if (any(is.na(final_df$Block_Score))) {
  min_block <- min(final_df$Block_Score, na.rm = TRUE)
  final_df[is.na(Block_Score), Block_Score := min_block]
  log_msg(sprintf("Remaining NA Block_Score filled with min value: %.6f", min_block))
}

# ---------- 稳健综合得分：Z-score 平均（排除 Reference 对 mean/sd 的影响） ----------
is_ref   <- final_df$Source == "Reference"
train_idx <- !is_ref

# 用真实候选点计算 mean / sd
mean_block <- mean(final_df$Block_Score[train_idx],  na.rm = TRUE)
sd_block   <- sd(final_df$Block_Score[train_idx],    na.rm = TRUE)
mean_spear <- mean(final_df$Spearman_Q1[train_idx],  na.rm = TRUE)
sd_spear   <- sd(final_df$Spearman_Q1[train_idx],    na.rm = TRUE)
mean_auc   <- mean(final_df$AUC_Q1[train_idx],       na.rm = TRUE)
sd_auc     <- sd(final_df$AUC_Q1[train_idx],         na.rm = TRUE)

# 对全部点（含 Reference）做标准化
z_block <- (final_df$Block_Score  - mean_block) / sd_block
z_spear <- (final_df$Spearman_Q1 - mean_spear) / sd_spear
z_auc   <- (final_df$AUC_Q1      - mean_auc)   / sd_auc

# 等权平均（可按需调整权重，例如 0.4*Block + 0.3*Spear + 0.3*AUC）
final_df[, Composite_raw := z_block + z_spear + z_auc]
final_df[, Composite_Score := min_max_scale(Composite_raw)]

setorder(final_df, -Composite_Score)

# 创建用于可视化的统一 Label
final_df[, WL_Label := sprintf("W%d_L%d", W, L)]
eval_merged[, WL_Label := sprintf("W%d_L%d", W, L)]

fwrite(final_df, file.path(OUT_DIR, "01.WL_Final_Ranking.csv"))
log_msg(sprintf("Final ranking saved successfully. Top parameter is %s (Composite=%.4f)",
                final_df$WL_Label[1], final_df$Composite_Score[1]))

# ================= 3. 数据可视化 =================
log_msg("Step 3: Generating scientific visualizations (Width <= 8 inches)...")

# 统一按 Composite_Score 降序排列所有图表的 X 轴
ordered_labels_comp <- final_df$WL_Label
eval_merged$WL_Label <- factor(eval_merged$WL_Label, levels = ordered_labels_comp)

# --- 图 1: AUC 与 Spearman 分布箱线图 (带双基准线) ---
p1 <- ggplot(eval_merged, aes(x = WL_Label, y = AUC)) +
  geom_boxplot(fill = "#5D9BCA", outlier.size = 1) +
  geom_hline(yintercept = ref_U8$AUC,  linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_hline(yintercept = ref_U16$AUC, linetype = "dashed", color = "#27AE60", linewidth = 0.8) +
  annotate("text", x = 1, y = ref_U8$AUC + 0.005,  label = "W80_L40_U8",  color = "#E74C3C", hjust = 0, fontface = "bold", size = 3) +
  annotate("text", x = 1, y = ref_U16$AUC - 0.005, label = "W80_L40_U16", color = "#27AE60", hjust = 0, fontface = "bold", size = 3) +
  theme_bw(base_size = 10) +
  labs(title = "AUC Distribution across U (6-16)", y = "Binary AUC", x = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

p2 <- ggplot(eval_merged, aes(x = WL_Label, y = Spearman_r)) +
  geom_boxplot(fill = "#E67E22", outlier.size = 1) +
  geom_hline(yintercept = ref_U8$Spearman_r,  linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_hline(yintercept = ref_U16$Spearman_r, linetype = "dashed", color = "#27AE60", linewidth = 0.8) +
  annotate("text", x = 1, y = ref_U8$Spearman_r + 0.005,  label = "W80_L40_U8",  color = "#E74C3C", hjust = 0, fontface = "bold", size = 3) +
  annotate("text", x = 1, y = ref_U16$Spearman_r - 0.005, label = "W80_L40_U16", color = "#27AE60", hjust = 0, fontface = "bold", size = 3) +
  theme_bw(base_size = 10) +
  labs(title = "Spearman r Distribution across U (6-16)",
       y = "Spearman r", x = "Candidates Ranked by Composite Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", size = 11))

fig1 <- ggarrange(p1, p2, ncol = 1, align = "v", heights = c(1, 1.2))
ggsave(file.path(OUT_DIR, "Fig1_Performance_Distribution.pdf"), fig1, width = 8, height = 7)


# --- 图 2: 泛化相关性散点图 (基于 Q1 鲁棒性指标) ---
plot_scatter <- function(data, x_var, y_var, x_lab, y_lab, title) {
  ct <- cor.test(data[[x_var]], data[[y_var]], use = "complete.obs")
  anno_text <- sprintf("r: %.3f\np: %.2e", ct$estimate, ct$p.value)
  
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(size = 2.5, alpha = 0.7, shape = 21, fill = "#8E44AD", color = "black") +
    geom_smooth(method = "lm", color = "red", linetype = "dashed",
                fill = "red", alpha = 0.15, linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, label = anno_text,
             hjust = -0.1, vjust = 1.2, fontface = "bold", size = 3) +
    theme_bw(base_size = 10) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

# 绘制散点图（仅用非 Reference 点更干净，也可保留全部）
plot_df <- final_df[!is.na(Block_Score)]
s1 <- plot_scatter(plot_df, "Block_Score", "AUC_Q1",
                   "Block Score", "siRNA AUC (Q1)", "Generalization: Block vs AUC")
s2 <- plot_scatter(plot_df, "Block_Score", "Spearman_Q1",
                   "Block Score", "siRNA Spearman (Q1)", "Generalization: Block vs Spearman")
s3 <- plot_scatter(plot_df, "Spearman_Q1", "AUC_Q1",
                   "siRNA Spearman (Q1)", "siRNA AUC (Q1)", "Internal Consistency")

fig2 <- ggarrange(s1, s2, s3, ncol = 2, nrow = 2)
ggsave(file.path(OUT_DIR, "Fig2_Generalization_Scatter.pdf"), fig2, width = 8, height = 7)


# --- 图 3: U 参数与 (W,L) 表现热图 ---
h1 <- ggplot(eval_merged, aes(x = WL_Label, y = factor(U), fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", AUC)), size = 2) +
  scale_fill_viridis_c(option = "viridis", name = "AUC") +
  theme_minimal(base_size = 10) +
  labs(title = "AUC Heatmap", x = "", y = "Target Length (U)") +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

h2 <- ggplot(eval_merged, aes(x = WL_Label, y = factor(U), fill = Spearman_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Spearman_r)), size = 2) +
  scale_fill_viridis_c(option = "magma", name = "Spearman r") +
  theme_minimal(base_size = 10) +
  labs(title = "Spearman r Heatmap",
       x = "Candidates Ranked by Composite Score", y = "Target Length (U)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", size = 11))

fig3 <- ggarrange(h1, h2, ncol = 1, align = "v")
ggsave(file.path(OUT_DIR, "Fig3_Parameter_Landscape_Heatmap.pdf"), fig3, width = 8, height = 9)

log_msg("Module partC.12 finished successfully.")