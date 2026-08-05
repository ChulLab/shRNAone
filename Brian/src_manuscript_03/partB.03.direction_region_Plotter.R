#!/usr/bin/env Rscript
# partB.03.direction_region_Plotter.R — Dedicated Visualization Script
#
# Generates highly polished, publication-ready supplementary figures 
# reading directly from the high-speed computational cache.

library(ggplot2)
library(dplyr)
library(data.table)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
CACHE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/03.direction_region/cache")
# OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/03.direction_region/cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")
REGION_COLORS <- c("3prime" = "#E41A1C", "5prime" = "#377EB8", "upstream" = "#4DAF4A", "downstream" = "#984EA3")

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# 验证缓存文件是否存在
req_files <- c("01_sig_data_cache.rds", "03_Variance_Partitioning_EtaSq.csv", "04_LMM_Contrasts.csv")
for (f in req_files) {
  if (!file.exists(file.path(CACHE_DIR, f))) stop(sprintf("Cache file missing: %s. Run compute core first.", f))
}

# ================= 1. ECDF 分布图 =================
log_msg("Rendering ECDF Distribution Plot...")
sig_data <- readRDS(file.path(CACHE_DIR, "01_sig_data_cache.rds"))

p_ecdf <- ggplot(sig_data, aes(x = spearman_r, color = direction)) +
  stat_ecdf(linewidth = 0.6, alpha = 0.8) +
  facet_grid(region ~ tool) +
  scale_color_manual(values = REGION_COLORS) +
  labs(title = "Empirical Cumulative Distribution of Spearman r",
       subtitle = "FDR < 0.05, Positive Correlations Only",
       x = "Spearman r", y = "Cumulative Probability") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.title = element_text(margin = margin(b = 2)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b = 6))
  )

out_p1 <- file.path(OUT_DIR, "Supple_01_ECDF_Distribution.pdf")
ggsave(out_p1, plot = p_ecdf, width = 6.5, height = 5, device = cairo_pdf)
log_msg(sprintf("Saved: %s", out_p1))


# ================= 2. 线性混合效应模型 (LMM) 森林图 =================
log_msg("Rendering LMM Forest Plot...")
lmm_dt <- fread(file.path(CACHE_DIR, "04_LMM_Contrasts.csv"))

# 提取 3prime 作为基准的对比项
forest_df <- lmm_dt[grepl("3prime -", contrast)]
forest_df[, Sig_Level := ifelse(p.value < 0.001, "P < 0.001", "ns")]

p_forest <- ggplot(forest_df, aes(x = estimate, y = contrast, color = Sig_Level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_point(size = 1.2) +
  geom_errorbarh(aes(xmin = estimate - 1.96*SE, xmax = estimate + 1.96*SE), height = 0.2, linewidth = 0.4) +
  facet_grid(region ~ Tool, scales = "free_y") +
  scale_color_manual(values = c("P < 0.001" = "#E41A1C", "ns" = "black"), name = "Significance") +
  labs(title = "LMM Fixed Effects: 3prime vs Other Directions",
       subtitle = "Positive estimate indicates higher correlation yield for 3prime",
       x = "Estimate Difference (3prime vs Others)", y = "Contrast") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.title = element_text(margin = margin(b = 2)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b = 6)),
    legend.position = "right"
  )

out_p2 <- file.path(OUT_DIR, "Supple_02_LMM_Forest.pdf")
ggsave(out_p2, plot = p_forest, width = 6.5, height = 4, device = cairo_pdf)
log_msg(sprintf("Saved: %s", out_p2))


# ================= 3. [NEW] 方差分解图 (Eta-Squared) =================
log_msg("Rendering Variance Partitioning Plot...")
eta_dt <- fread(file.path(CACHE_DIR, "03_Variance_Partitioning_EtaSq.csv"))

eta_long <- melt(eta_dt, id.vars = "Tool", variable.name = "Source", value.name = "Variance_Percentage")
eta_long[, Source := gsub("Variance_Explained_by_", "", Source)]
eta_long[, Source := factor(Source, levels = c("Interaction", "Region", "Direction"))]

p_eta <- ggplot(eta_long, aes(x = Tool, y = Variance_Percentage, fill = Source)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.6) +
  scale_fill_manual(values = c("Direction" = "#3498DB", "Region" = "#F39C12", "Interaction" = "#95A5A6")) +
  geom_text(aes(label = sprintf("%.1f%%", Variance_Percentage)), 
            position = position_stack(vjust = 0.5), size = 8 / .pt, family = "Arial", color = "white") +
  labs(title = "Variance Explained by Grid Parameters (Eta-squared)",
       subtitle = "Percentage of variance uniquely explained by targeting rules",
       x = "Tool", y = "Variance Explained (%)") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black", size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.title = element_text(margin = margin(b = 2)),
    plot.subtitle = element_text(color = "gray40", margin = margin(b = 6)),
    legend.position = "right",
    legend.key.size = unit(3, "mm")
  )

out_p3 <- file.path(OUT_DIR, "Supple_03_Variance_EtaSq.pdf")
ggsave(out_p3, plot = p_eta, width = 5.5, height = 3.5, device = cairo_pdf)
log_msg(sprintf("Saved: %s", out_p3))


# ================= 4. Top 5% 超几何分布富集柱状图 =================
log_msg("Rendering Top 5% Dominance Plot...")

top_counts <- sig_data %>%
  group_by(tool, region) %>%
  filter(spearman_r >= quantile(spearman_r, 0.95)) %>%
  group_by(tool, region, direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(tool, region) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  as.data.table()

p_enrich <- ggplot(top_counts, aes(x = region, y = Proportion, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  facet_wrap(~ tool) +
  scale_fill_manual(values = REGION_COLORS) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 1)), 
            position = position_stack(vjust = 0.5), size = 8 / .pt, family = "Arial", color = "white") +
  labs(title = "Dominance in Top 5% Performance Parameter Space",
       x = "Region", y = "Proportion of Top 5% Configurations") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", size = 8),
    axis.text = element_text(color = "black", size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.title = element_text(margin = margin(b = 6)),
    legend.position = "right"
  )

out_p4 <- file.path(OUT_DIR, "Supple_04_Top5_Dominance.pdf")
ggsave(out_p4, plot = p_enrich, width = 6.5, height = 3.5, device = cairo_pdf)
log_msg(sprintf("Saved: %s", out_p4))

log_msg("--- All Plotting Tasks Completed Successfully ---")