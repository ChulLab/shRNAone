#!/usr/bin/env Rscript
# derived from: partF.04.example_shiny_track.R
#
# SF2G: Example Track Plot for RPLP0 with Top 50 Overlay (W35_L20_U11 vs W80_L40_U17)

library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
TRANSCRIPT_CSV <- file.path(PROJECT, "result/08.manuscript_03/partF/00.raw_data/transcript_list.csv")
MATRIX_DIR <- file.path(PROJECT, "result/08.manuscript_03/partF/02.target_matrix")

OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_GENE <- "RPLP0"

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# 统一颜色映射
REGION_COLORS <- c("5'UTR" = "#3498DB", "CDS" = "#F1C40F", "3'UTR" = "#E74C3C")
AREA_FILLS <- c("W35_L20_U11 (Optimal)" = "#A3E4D7", "W80_L40_U17 (Baseline)" = "#D5D8DC")
LINE_COLORS <- c("W35_L20_U11 (Optimal)" = "#117A65", "W80_L40_U17 (Baseline)" = "#5D6D7E")

# ================= 1. 数据定位 =================
log_msg(sprintf("Searching for gene: %s", TARGET_GENE))
tx_dt <- fread(TRANSCRIPT_CSV)

target_info <- tx_dt[gene_name == TARGET_GENE][1]
if (nrow(target_info) == 0 || is.na(target_info$transcript_accession)) {
  stop(sprintf("Gene %s not found in transcript_list.csv!", TARGET_GENE))
}

acc_id <- target_info$transcript_accession
csv_path <- file.path(MATRIX_DIR, paste0(acc_id, ".csv"))

# ================= 2. 数据处理与长格式转换 =================
df <- fread(csv_path)

df[, acc_W80_L40_U17 := as.numeric(acc_W80_L40_U17)]
df[, acc_W35_L20_U11 := as.numeric(acc_W35_L20_U11)]
df[, rank_W80_L40_U17 := as.numeric(rank_W80_L40_U17)]
df[, rank_W35_L20_U11 := as.numeric(rank_W35_L20_U11)]

df_long <- melt(
  df, 
  id.vars = c("transcript_version", "start_pos", "end_pos", "target_region", "target_sequence"),
  measure.vars = list(
    Score = c("acc_W35_L20_U11", "acc_W80_L40_U17"),
    Rank  = c("rank_W35_L20_U11", "rank_W80_L40_U17")
  ),
  variable.name = "Parameter",
  value.name = c("Score", "Rank")
)

df_long$Parameter <- factor(df_long$Parameter, levels = c(1, 2), labels = c("W35_L20_U11 (Optimal)", "W80_L40_U17 (Baseline)"))
df_long$target_region <- factor(df_long$target_region, levels = c("5'UTR", "CDS", "3'UTR"))

pseudo_count <- 1e-4
df_long$Plot_Score <- log10(df_long$Score + pseudo_count)

mean_dt <- df_long[!is.na(Score), .(Mean_Score = mean(Score)), by = Parameter]
mean_dt$Plot_Mean <- log10(mean_dt$Mean_Score + pseudo_count)

top50_dt <- df_long[!is.na(Rank) & Rank <= 50]

# ================= 3. 绘制左图：双轨 Track Plot =================
rug_ymin <- -4.25
rug_ymax <- -4.10

p_track <- ggplot(df_long, aes(x = end_pos, y = Plot_Score)) +
  geom_ribbon(aes(ymin = -4, ymax = Plot_Score, fill = Parameter), alpha = 0.5) +
  geom_line(aes(color = Parameter), linewidth = 0.3) +
  
  geom_segment(aes(xend = end_pos, y = rug_ymin, yend = rug_ymax, color = target_region), linewidth = 1) +
  
  geom_hline(data = mean_dt, aes(yintercept = Plot_Mean), linetype = "dashed", color = "#C0392B", linewidth = 0.4) +
  
  geom_text(data = mean_dt, aes(x = max(df_long$end_pos), y = Plot_Mean, label = sprintf("Mean: %.3f", Mean_Score)), 
            hjust = 1, vjust = -0.5, color = "#C0392B", size = 8 / .pt, family = "Arial") +
  
  geom_point(data = top50_dt, aes(fill = target_region), shape = 21, color = "black", size = 1.2, stroke = 0.3) +
  facet_wrap(~ Parameter, ncol = 1) +
  
  scale_fill_manual(values = c(AREA_FILLS, REGION_COLORS)) +
  scale_color_manual(values = c(LINE_COLORS, REGION_COLORS)) +
  
  # 【核心修改】：还原真实的负数对数标签
  scale_y_continuous(breaks = c(-4, -3, -2, -1, 0), labels = c("-4", "-3", "-2", "-1", "0")) +
  coord_cartesian(ylim = c(-4.25, 0.2)) +
  
  labs(title = sprintf("Accessibility Landscape (%s)", TARGET_GENE), x = "Transcript Position", y = "Predicted Score (Log10)") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "#EBEDEF", color = "black", linewidth = 0.3),
    strip.text = element_text(size = 8, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, margin = margin(b = 6))
  )

# ================= 4. 绘制右上图：共同选点散点图 =================
df$log_W80 <- log10(df$acc_W80_L40_U17 + pseudo_count)
df$log_W35 <- log10(df$acc_W35_L20_U11 + pseudo_count)

thresh_W35 <- min(df[rank_W35_L20_U11 <= 50]$acc_W35_L20_U11, na.rm = TRUE)
thresh_W80 <- min(df[rank_W80_L40_U17 <= 50]$acc_W80_L40_U17, na.rm = TRUE)
log_thresh_W35 <- log10(thresh_W35 + pseudo_count)
log_thresh_W80 <- log10(thresh_W80 + pseudo_count)

df[, Top_Category := "Others"]
df[rank_W35_L20_U11 <= 50 & rank_W80_L40_U17 <= 50, Top_Category := "Shared Top 50"]
df[rank_W35_L20_U11 <= 50 & rank_W80_L40_U17 > 50, Top_Category := "Optimal Top 50 Only"]
df[rank_W35_L20_U11 > 50 & rank_W80_L40_U17 <= 50, Top_Category := "Baseline Top 50 Only"]

df$Top_Category <- factor(df$Top_Category, levels = c("Shared Top 50", "Optimal Top 50 Only", "Baseline Top 50 Only", "Others"))
df$target_region <- factor(df$target_region, levels = c("5'UTR", "CDS", "3'UTR"))

df_others <- df[Top_Category == "Others"]
df_highlights <- df[Top_Category != "Others"]

p_scatter <- ggplot() +
  geom_point(data = df_others, aes(x = log_W80, y = log_W35), color = "gray80", size = 0.3) +
  geom_vline(xintercept = log_thresh_W80, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_hline(yintercept = log_thresh_W35, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black", linewidth = 0.3) +
  geom_point(data = df_highlights, aes(x = log_W80, y = log_W35, fill = target_region, shape = Top_Category), 
             size = 1.2, color = "black", stroke = 0.2) +
  
  scale_shape_manual(values = c("Shared Top 50" = 21, "Optimal Top 50 Only" = 24, "Baseline Top 50 Only" = 22)) +
  scale_fill_manual(values = REGION_COLORS) +
  
  # 【核心修改】：散点图坐标轴同步还原为真实的负数对数标签
  scale_x_continuous(breaks = c(-4, -2, 0), labels = c("-4", "-2", "0"), limits = c(-4.2, 0.2)) +
  scale_y_continuous(breaks = c(-4, -2, 0), labels = c("-4", "-2", "0"), limits = c(-4.2, 0.2)) +
  
  labs(title = "Top 50 Agreement", x = "Baseline Acc. (Log10)", y = "Optimal Acc. (Log10)", fill = "Region", shape = "Selection") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, margin = margin(b = 2)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = "gray60", linewidth = 0.2),
    legend.key.size = unit(2.5, "mm"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.margin = margin(1, 2, 1, 2),
    legend.spacing.y = unit(0.5, "mm")
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 1.5), order = 1),
         shape = guide_legend(override.aes = list(fill = "gray50", size = 1.5), order = 2))

# ================= 5. 绘制右下图：Top 50 区域分布柱状图 =================
bar_dt <- top50_dt[, .N, by = .(Parameter, target_region)]
bar_dt[, Param_Short := ifelse(grepl("W35", Parameter), "Optimal", "Baseline")]
bar_dt$Param_Short <- factor(bar_dt$Param_Short, levels = c("Optimal", "Baseline"))

p_bar <- ggplot(bar_dt, aes(x = Param_Short, y = N, fill = target_region)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.5) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 8 / .pt, family = "Arial") +
  scale_fill_manual(values = REGION_COLORS) +
  labs(title = "Top 50 Regional Dist.", x = NULL, y = "Count") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, margin = margin(b = 2)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

# ================= 6. 组装与保存 =================
right_col <- ggarrange(p_scatter, p_bar, ncol = 1, nrow = 2, align = "v", heights = c(1.3, 1))
merged_plot <- ggarrange(p_track, right_col, ncol = 2, nrow = 1, widths = c(2.5, 1))

out_pdf <- file.path(OUT_DIR, "Task6_Example_Track_Joint.pdf")
ggsave(out_pdf, plot = merged_plot, width = 7.5, height = 3, device = cairo_pdf)

# --- Save summary stats ---
mean_dt[, Param_Short := ifelse(grepl("W35", Parameter), "Optimal_W35L20", "Baseline_W80L40")]
fwrite(mean_dt[, .(Param_Short, Mean_Score)], file.path(OUT_DIR, "SF2G_RPLP0_Mean_Scores.csv"))
fwrite(bar_dt, file.path(OUT_DIR, "SF2G_RPLP0_Top50_Regional_Dist.csv"))
fwrite(data.table(
  Threshold = c("W35_Top50", "W80_Top50"),
  Raw_Value = c(thresh_W35, thresh_W80),
  Log10 = c(log_thresh_W35, log_thresh_W80)
), file.path(OUT_DIR, "SF2G_RPLP0_Top50_Thresholds.csv"))
log_msg(sprintf("Final plot assembled and saved to: %s (Width: 7.5, Height: 3)", out_pdf))