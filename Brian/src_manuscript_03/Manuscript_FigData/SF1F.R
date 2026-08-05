#!/usr/bin/env Rscript
# derived from: partC.12.WL_compare.R
# 
# SF1F: select best W/L candidates based on composite score and visualize their performance on GT siRNA dataset

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. Configuration & Initialization =================
BASE_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partC"
EVAL_DIR <- file.path(BASE_DIR, "03.GT_siRNA_acc", "02.siRNA_matrix")
CAND_FILE <- file.path(BASE_DIR, "01.Self_WL_list", "03.selected_WL_Candidates.csv")
OUT_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/Manuscript_FigData/supple"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE_SIZE <- 6.5 # Base font size under extreme size constraint
GRAY_PALETTE <- c("Score" = "gray40", "Box1" = "gray70", "Box2" = "gray85", "Ref" = "black")

# Highlight colors for reference lines
COLOR_U8 <- "#E74C3C"  # Bright red
COLOR_U16 <- "#C0392B" # Dark red

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
min_max_scale <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# ================= 1. Data Loading & Factor Ordering =================
log_msg("Step 1: Loading GT data and pre-calculated WL rankings...")

cont_df <- fread(file.path(EVAL_DIR, "01.Evaluation_Continuous.csv"))[U >= 6 & U <= 16]
bin_df  <- fread(file.path(EVAL_DIR, "02.Evaluation_Binary.csv"))[U >= 6 & U <= 16]

eval_merged <- merge(
  cont_df[, .(W, L, U, Spearman_r)],
  bin_df[, .(W, L, U, AUC)],
  by = c("W", "L", "U")
)
eval_merged[, WL_Label := sprintf("W%d_L%d", W, L)]

# 【核心修正】：直接读取 partC.12 已经计算好的精确 Ranking 文件
RANKING_FILE <- file.path(BASE_DIR, "03.GT_siRNA_acc", "03.WL_selection", "01.WL_Final_Ranking.csv")
final_df <- fread(RANKING_FILE)

# 提取核心对照组用于画基准线
ref_U8  <- eval_merged[W == 80 & L == 40 & U == 8]
ref_U16 <- eval_merged[W == 80 & L == 40 & U == 16]

# 【排版逻辑】：强制把 Reference (W80_L40) 放在 X 轴最前面，其余候选点按真实得分降序排列
cand_labels <- final_df[WL_Label != "W80_L40"][order(-Composite_Score)]$WL_Label
ordered_labels <- c("W80_L40", cand_labels)

final_df$WL_Label    <- factor(final_df$WL_Label, levels = ordered_labels)
eval_merged$WL_Label <- factor(eval_merged$WL_Label, levels = ordered_labels)

# ================= 2. Compact Grayscale 3-Panel Figure =================
log_msg("Step 2: Rendering compact 3-panel grayscale plots with RED baselines...")

# Panel a (Top): Spearman boxplot with red baselines
p_spearman <- ggplot(eval_merged, aes(x = WL_Label, y = Spearman_r)) +
  geom_boxplot(fill = GRAY_PALETTE["Box1"], outlier.size = 0.5, outlier.shape = 1, linewidth = 0.3) +
  geom_hline(yintercept = ref_U8$Spearman_r,  linetype = "dashed", color = COLOR_U8,  linewidth = 0.5) +
  geom_hline(yintercept = ref_U16$Spearman_r, linetype = "dotted", color = COLOR_U16, linewidth = 0.7) +
  annotate("text", x = 2.5, y = ref_U8$Spearman_r + 0.006,  label = "W80_L40_U8",  color = COLOR_U8,  hjust = 0, fontface = "bold", size = 2) +
  annotate("text", x = 2.5, y = ref_U16$Spearman_r - 0.006, label = "W80_L40_U16", color = COLOR_U16, hjust = 0, fontface = "bold", size = 2) +
  labs(title = "a  Generalization on Ground Truth (Spearman r)", y = "Spearman", x = NULL) +
  theme_bw(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 2)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(2, 2, 0, 2)
  )

# Panel b (Middle): AUC boxplot with red baselines
p_auc <- ggplot(eval_merged, aes(x = WL_Label, y = AUC)) +
  geom_boxplot(fill = GRAY_PALETTE["Box2"], outlier.size = 0.5, outlier.shape = 1, linewidth = 0.3) +
  geom_hline(yintercept = ref_U8$AUC,  linetype = "dashed", color = COLOR_U8,  linewidth = 0.5) +
  geom_hline(yintercept = ref_U16$AUC, linetype = "dotted", color = COLOR_U16, linewidth = 0.7) +
  annotate("text", x = 2.5, y = ref_U8$AUC + 0.005,  label = "W80_L40_U8",  color = COLOR_U8,  hjust = 0, fontface = "bold", size = 2) +
  annotate("text", x = 2.5, y = ref_U16$AUC - 0.006, label = "W80_L40_U16", color = COLOR_U16, hjust = 0, fontface = "bold", size = 2) +
  labs(title = "b  Classification Power on Ground Truth (AUC)", y = "Binary AUC", x = NULL) +
  theme_bw(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", margin = margin(t = 2, b = 2)),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(0, 2, 0, 2)
  )

# Panel c (Bottom): Composite Score barplot
p_score <- ggplot(final_df, aes(x = WL_Label, y = Composite_Score)) +
  geom_bar(stat = "identity",
           fill = ifelse(final_df$WL_Label == "W80_L40", NA, GRAY_PALETTE["Score"]),
           color = "black", linewidth = 0.3, width = 0.7) +
  labs(title = "c  Global Parameter Ranking", y = "Composite Score", x = "Targeting Optimization Configurations") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_bw(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", margin = margin(t = 2, b = 2)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black",
                               face = ifelse(final_df$WL_Label == "W80_L40", "bold.italic", "plain")),
    axis.title.x = element_text(face = "bold", margin = margin(t = 2)),
    panel.grid = element_blank(),
    plot.margin = margin(0, 2, 2, 2)
  )

# Merge panels (top to bottom: Spearman → AUC → Composite)
fig_merge <- ggarrange(p_spearman, p_auc, p_score, ncol = 1, align = "v", heights = c(1.35, 1.35, 1.3))
out_pdf <- file.path(OUT_DIR, "Supple_WL_Selection_3Panel.pdf")
ggsave(out_pdf, plot = fig_merge, width = 4, height = 2.8, device = cairo_pdf)

# --- Save siRNA evaluation data ---
siRNA_summary <- eval_merged[, .(
  Spearman_Q1 = quantile(Spearman_r, 0.25, na.rm = TRUE),
  Spearman_Median = quantile(Spearman_r, 0.50, na.rm = TRUE),
  Spearman_Q3 = quantile(Spearman_r, 0.75, na.rm = TRUE),
  AUC_Q1 = quantile(AUC, 0.25, na.rm = TRUE),
  AUC_Median = quantile(AUC, 0.50, na.rm = TRUE),
  AUC_Q3 = quantile(AUC, 0.75, na.rm = TRUE)
), by = .(W, L, WL_Label)]
fwrite(siRNA_summary, file.path(OUT_DIR, "SF1F_siRNA_Distribution_Stats.csv"))
fwrite(final_df[, .(WL_Label, W, L, Composite_Score, Block_Score, Spearman_Q1, AUC_Q1)], 
       file.path(OUT_DIR, "SF1F_Composite_Scores.csv"))
fwrite(data.table(
  Reference = c("W80_L40_U8", "W80_L40_U16"),
  Spearman_r = c(ref_U8$Spearman_r, ref_U16$Spearman_r),
  AUC = c(ref_U8$AUC, ref_U16$AUC)
), file.path(OUT_DIR, "SF1F_Reference_Baselines.csv"))
log_msg(sprintf("Saved successfully to: %s", out_pdf))