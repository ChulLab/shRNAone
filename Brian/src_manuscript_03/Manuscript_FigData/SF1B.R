#!/usr/bin/env Rscript
# derived from: partB.03.direction_region_Plotter.R
# 
# SF1B_modified: Direction distribution (Whole only) for shRNAone, CasRx, and PspCas13b

library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)
library(tidyr)

# ================= 0. Configuration & Initialization =================
PROJECT <- "/data/cai801/data/HKUcas"
CACHE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/03.direction_region/cache")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

REGION_COLORS <- c("Whole" = "gray70", "5'UTR" = "#3498DB", "CDS" = "#21908C", "3'UTR" = "#FDE725")
DIRECTION_COLORS_MUTED <- c("3prime" = "#D28E72", 
                            "5prime" = "#7A9E9F", 
                            "upstream" = "#B3BCA4", 
                            "downstream" = "#C5B4A3")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. Upper Panel: Direction Proportion (Whole only) =================
log_msg("Rendering Panel A: Direction Proportion in Top 5% (Whole region)...")
sig_data <- readRDS(file.path(CACHE_DIR, "01_sig_data_cache.rds"))

sig_data <- sig_data %>%
  mutate(tool = factor(tool, levels = c("shRNAone", "CasRx_day5", "PspCas13b")))

# Restrict to Whole region only
sig_whole <- sig_data %>% filter(region == "whole")

top_counts <- sig_whole %>%
  group_by(tool) %>%
  filter(spearman_r >= quantile(spearman_r, 0.95)) %>%
  group_by(tool, direction) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(tool) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  as.data.table()

# Ensure all 4 directions exist for every tool (fill missing with 0)
top_counts <- as.data.table(
  tidyr::complete(as.data.frame(top_counts), tool, direction,
                  fill = list(Count = 0L, Proportion = 0))
)
top_counts[, direction := factor(direction, levels = c("3prime", "5prime", "upstream", "downstream"))]
top_counts[, tool := factor(tool, levels = c("shRNAone", "CasRx_day5", "PspCas13b"))]

top_counts[, label_text := ifelse(Proportion > 0.03, 
                                   scales::percent(Proportion, accuracy = 1), "")]

p_enrich <- ggplot(top_counts, aes(x = direction, y = Proportion, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.6) +
  facet_wrap(~ tool, nrow = 1) +
  scale_fill_manual(values = DIRECTION_COLORS_MUTED) +
  geom_text(aes(label = label_text, y = Proportion + 0.02),
            size = 7 / .pt, family = "Arial", color = "black") +
  labs(x = NULL, y = "Top 5% Proportion", fill = "") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.25),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(t = 2, r = 2, b = 0, l = 2)
  )

# ================= 2. Lower Panel: LMM Fixed Effects (Whole only) =================
log_msg("Rendering Panel B: LMM Fixed Effects (Whole region only)...")
lmm_dt <- fread(file.path(CACHE_DIR, "04_LMM_Contrasts.csv"))

lmm_dt[, Tool := factor(Tool, levels = c("shRNAone", "CasRx_day5", "PspCas13b"))]

forest_df <- lmm_dt[grepl("3prime -", contrast)]
forest_df <- forest_df[region == "whole"]

forest_df[, Contrast_Short := gsub("3prime - ", "", contrast)]
forest_df[, Contrast_Short := factor(Contrast_Short, levels = c("5prime", "upstream", "downstream"))]

forest_df[, Sig_Label := ifelse(p.value < 0.001, "***", 
                                ifelse(p.value < 0.01, "**", 
                                       ifelse(p.value < 0.05, "*", "")))]
forest_df[, y_pos := estimate + 1.96*SE + 0.01]

p_lmm <- ggplot(forest_df, aes(x = Contrast_Short, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_bar(stat = "identity", fill = REGION_COLORS["Whole"],
           color = "black", linewidth = 0.2, width = 0.6) +
  geom_errorbar(aes(ymin = estimate - 1.96*SE, ymax = estimate + 1.96*SE),
                width = 0.3, linewidth = 0.2) +
  geom_text(aes(y = y_pos, label = Sig_Label),
            size = 7 / .pt, family = "Arial", vjust = 0, hjust = 0.5, angle = 0) +
  facet_wrap(~ Tool, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "Contrast (3prime vs Others)", y = "Est. Diff (LMM)", fill = "") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(t = 0, r = 2, b = 2, l = 2)
  )

# ================= 3. Legend Reorganization & Merged Output =================
log_msg("Merging plots and formatting legends...")

leg_top <- get_legend(p_enrich)

p_enrich_clean <- p_enrich + theme(legend.position = "none")
p_lmm_clean <- p_lmm   # already has no legend

main_plots <- ggarrange(p_enrich_clean, p_lmm_clean, 
                        ncol = 1, nrow = 2, 
                        align = "v", 
                        heights = c(1, 1.25))

merged_plot <- ggarrange(main_plots, leg_top, 
                         ncol = 1, nrow = 2, 
                         heights = c(10, 1))

# ================= 4. Global Title Addition =================
title_str <- sprintf(
  "WLU Scan Scope (W: 20-200, L: 0.25W-0.75W, U: 30) | Significance Filter: FDR < 0.05, r > 0\nValid Configurations (N): shRNAone = %d, CasRx_day5 = %d, PspCas13b = %d",
  sum(sig_data$tool == "shRNAone"),
  sum(sig_data$tool == "CasRx_day5"),
  sum(sig_data$tool == "PspCas13b")
)

final_plot_with_title <- annotate_figure(
  merged_plot,
  top = text_grob(title_str, face = "bold", size = 9, family = "Arial", lineheight = 1.2)
)

out_pdf <- file.path(OUT_DIR, "Supple_Direction_Region_Unified.pdf")
ggsave(out_pdf, plot = final_plot_with_title, width = 4.5, height = 4.0, device = cairo_pdf)

# --- Save direction and LMM stats ---
fwrite(top_counts, file.path(OUT_DIR, "SF1B_Top5_Direction_Proportions.csv"))
fwrite(forest_df, file.path(OUT_DIR, "SF1B_LMM_Whole_Contrasts.csv"))
log_msg(sprintf("Plot saved to: %s  |  CSV stats saved.", out_pdf))
