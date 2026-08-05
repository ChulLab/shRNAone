#!/usr/bin/env Rscript
# partF.04.example_shiny_track.R
#

library(data.table)
library(ggplot2)
library(dplyr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partF"
TRANSCRIPT_CSV <- file.path(PROJECT, "00.raw_data", "transcript_list.csv")
MATRIX_DIR <- file.path(PROJECT, "02.target_matrix")
OUT_DIR <- file.path(PROJECT, "03.example")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_GENE <- "RPLP0"

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. 自动定位 RPLP0 数据 =================
log_msg(sprintf("Searching for gene: %s", TARGET_GENE))
tx_dt <- fread(TRANSCRIPT_CSV)

target_info <- tx_dt[gene_name == TARGET_GENE][1]
if (nrow(target_info) == 0 || is.na(target_info$transcript_accession)) {
  stop(sprintf("Gene %s not found in transcript_list.csv!", TARGET_GENE))
}

acc_id <- target_info$transcript_accession
log_msg(sprintf("Found %s: Accession = %s (Length = %d nt)", TARGET_GENE, acc_id, target_info$transcript_length))

csv_path <- file.path(MATRIX_DIR, paste0(acc_id, ".csv"))
if (!file.exists(csv_path)) {
  stop(sprintf("Target matrix for %s not found in %s", acc_id, MATRIX_DIR))
}

# ================= 2. 数据处理与导出 =================
log_msg("Loading target matrix and separating Top 50 / Full data...")
df <- fread(csv_path)

df[, acc_W80_L40_U17 := as.numeric(acc_W80_L40_U17)]
df[, acc_W35_L20_U11 := as.numeric(acc_W35_L20_U11)]
df[, rank_W80_L40_U17 := as.numeric(rank_W80_L40_U17)]
df[, rank_W35_L20_U11 := as.numeric(rank_W35_L20_U11)]

out_csv_full <- file.path(OUT_DIR, paste0("Shiny_Demo_Data_", TARGET_GENE, "_full.csv"))
fwrite(df, out_csv_full)
log_msg(sprintf("Full dataset saved to: %s", out_csv_full))

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

mean_dt <- df_long[!is.na(Score), .(Mean_Score = mean(Score)), by = Parameter]
top50_dt <- df_long[!is.na(Rank) & Rank <= 50]

out_csv_top50 <- file.path(OUT_DIR, paste0("Shiny_Demo_Data_", TARGET_GENE, ".csv"))
fwrite(top50_dt, out_csv_top50)

# ================= 3. 为 Log 坐标系进行数据预转换 =================
# 添加 1e-4 防止 0 值在 log10 下变成 -Inf 导致图形崩溃
pseudo_count <- 1e-4

df_long$Plot_Score <- log10(df_long$Score + pseudo_count)
mean_dt$Plot_Mean <- log10(mean_dt$Mean_Score + pseudo_count)
top50_dt$Plot_Score <- log10(top50_dt$Score + pseudo_count)

# ================= 4. 绘制高级轨道图 (Log 转换版) =================
log_msg("Rendering Shiny App Demo Track Plot in Log10 Scale...")

region_colors <- c("5'UTR" = "#3498DB", "CDS" = "#F1C40F", "3'UTR" = "#E74C3C")
area_fills <- c("W35_L20_U11 (Optimal)" = "#A3E4D7", "W80_L40_U17 (Baseline)" = "#D5D8DC")
line_colors <- c("W35_L20_U11 (Optimal)" = "#117A65", "W80_L40_U17 (Baseline)" = "#5D6D7E")

# 设置地基：log10(1e-4) = -4，所以将 rug 放在 -4.25 到 -4.10 之间
rug_ymin <- -4.25
rug_ymax <- -4.10

p <- ggplot(df_long, aes(x = start_pos, y = Plot_Score)) +
  # 使用 geom_ribbon 固定底线 ymin = -4 (即原始值的 0)
  geom_ribbon(aes(ymin = -4, ymax = Plot_Score, fill = Parameter), alpha = 0.5) +
  geom_line(aes(color = Parameter), linewidth = 0.4) +
  
  geom_segment(aes(xend = start_pos, y = rug_ymin, yend = rug_ymax, color = target_region), linewidth = 1.2) +
  
  geom_hline(data = mean_dt, aes(yintercept = Plot_Mean), 
             linetype = "dashed", color = "#C0392B", linewidth = 0.6) +
  geom_text(data = mean_dt, aes(x = max(df_long$start_pos), y = Plot_Mean, 
                                label = sprintf("Mean: %.3f", Mean_Score)), 
            hjust = 1, vjust = -0.5, color = "#C0392B", size = 3, fontface = "italic") +
  
  geom_point(data = top50_dt, aes(fill = target_region), shape = 21, color = "black", size = 1.8, stroke = 0.4) +
  
  facet_wrap(~ Parameter, ncol = 1) +
  
  scale_fill_manual(values = c(area_fills, region_colors)) +
  scale_color_manual(values = c(line_colors, region_colors)) +
  
  # 手动将对数坐标还原为人类易读的标签
  scale_y_continuous(
    breaks = c(-4, -3, -2, -1, 0),
    labels = c("0", "0.001", "0.01", "0.1", "1.0")
  ) +
  # 限制 Y 轴的视野范围，留出上下边距
  coord_cartesian(ylim = c(-4.25, 0.2)) +
  
  labs(
    title = sprintf("PspCas13b Target Accessibility Landscape (%s, %s)", TARGET_GENE, acc_id),
    subtitle = "Points indicate Top 50 targets. Displayed on a pseudo-Log10 scale to reveal hidden patterns.",
    x = "Transcript Position (3' end of target)",
    y = "Predicted Accessibility Score (Log10 scaled)",
    color = "Legend",
    fill = "Legend"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30", size = 10, margin = margin(b = 10)),
    strip.background = element_rect(fill = "#EBEDEF", color = NA),
    strip.text = element_text(face = "bold", size = 11, color = "#2C3E50"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank()
  ) +
  
  guides(color = guide_legend(override.aes = list(linewidth = 0)),
         fill = guide_legend(override.aes = list(shape = 21, size = 4)))

out_pdf <- file.path(OUT_DIR, paste0("Shiny_Track_Example_", TARGET_GENE, "_Log.pdf"))
ggsave(out_pdf, plot = p, width = 10, height = 6.5)

log_msg(sprintf("Plot successfully generated and saved to: %s", out_pdf))
log_msg("Part F.04 execution complete. The hidden baseline patterns are now revealed!")