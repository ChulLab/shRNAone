#!/usr/bin/env Rscript
# derived from: partA.14.plot_wet_dry.R
#
# SF2D: Wet vs Dry Correlation for CasRx and PspCas13b across all models

library(data.table)
library(ggplot2)
library(ggpubr)

# -------------------- paths --------------------
PROJECT <- "/data/cai801/data/HKUcas"
FIG_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 去除 shRNAone
TOOL_PREFIX <- c(CasRx_day5 = "02", PspCas13b = "04")

# -------------------- model appearance --------------------
# 替换为高对比度的科学配色，并且全部统一为实线 (solid)
MODEL_META <- data.table(
  model = c("DeepCas13_23", "DeepCas13_18", "DeepCas13_30",
            "Hsu2023_30",
            "Sanjana2020_23", "Sanjana2020_30",
            "TIGER2023_23",
            "Fareh2024_30"),
  short = c("DeepCas13", "DeepCas13", "DeepCas13",
            "Hsu2023",
            "Sanjana", "Sanjana",
            "TIGER",
            "Fareh2024"),
  col   = c("#E74C3C", "#E74C3C", "#E74C3C", # 红色
            "#3498DB",                       # 蓝色
            "#2ECC71", "#2ECC71",            # 绿色
            "#9B59B6",                       # 紫色
            "#F39C12"),                      # 橙色
  lty   = c("solid", "solid", "solid",       # 全部使用实线
            "solid",
            "solid", "solid",
            "solid",
            "solid")
)

# -------------------- load --------------------
dt <- rbindlist(lapply(names(TOOL_PREFIX), function(tool) {
  f <- file.path(FIG_DIR, sprintf("%s.offset_wet_dry_%s.csv", TOOL_PREFIX[tool], tool))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d[, tool := tool]
  d
}), use.names = TRUE, fill = TRUE)

dt <- merge(dt, MODEL_META, by = "model", all.x = TRUE)
dt[is.na(short), `:=`(short = model, col = "grey30", lty = "solid")]

lev <- unique(MODEL_META$short)
dt[, short := factor(short, levels = lev)]
dt[, tool  := factor(tool,  levels = names(TOOL_PREFIX))]

best <- dt[, .SD[which.max(correlation_r)], by = .(tool, model)]
best[, short := factor(short, levels = lev)]
best[, tool  := factor(tool,  levels = names(TOOL_PREFIX))]

# -------------------- maps --------------------
# 【修复核心】：按组提取 unique 映射，避免单列 unique 导致的长度不一致问题
meta_unique <- unique(MODEL_META[, .(short, col, lty)])
col_map <- setNames(meta_unique$col, meta_unique$short)
lty_map <- setNames(meta_unique$lty, meta_unique$short)

# -------------------- panel function --------------------
make_panel <- function(tool_name, show_x = FALSE) {
  d <- dt[tool == tool_name]
  b <- best[tool == tool_name]

  # sort by r high → low, assign evenly spaced y on the right
  setorder(b, -correlation_r)
  y_rng <- range(d$correlation_r, na.rm = TRUE)
  y_top <- y_rng[2] + diff(y_rng) * 0.06
  y_bot <- y_rng[1] + diff(y_rng) * 0.12
  b[, lab_y := seq(y_top, y_bot, length.out = .N)]
  b[, label := sprintf("(%+d, %.2f)", offset, correlation_r)]

  x_seg_end <- 52
  x_text    <- 54

  ggplot(d, aes(x = offset, y = correlation_r,
                colour = short, linetype = short, group = model)) +
    geom_vline(xintercept = 0, colour = "grey60",
               linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
    # 细线风格
    geom_line(linewidth = 0.35) + 
    geom_point(data = b, aes(x = offset, y = correlation_r, colour = short),
               shape = 21, size = 1.5, fill = "white", stroke = 0.6,
               inherit.aes = FALSE, show.legend = FALSE) +
    geom_segment(data = b,
                 aes(x = offset, xend = x_seg_end,
                     y = correlation_r, yend = lab_y),
                 colour = "grey50", linewidth = 0.3,
                 inherit.aes = FALSE) +
    geom_text(data = b,
              aes(x = x_text, y = lab_y, label = label),
              hjust = 0, size = 1.6, colour = "black",
              inherit.aes = FALSE) +
    scale_colour_manual(values = col_map, name = NULL) +
    scale_linetype_manual(values = lty_map, name = NULL) +
    scale_x_continuous(breaks = c(-40, 0, 40), expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    coord_cartesian(xlim = c(-50, 50), clip = "off") +
    labs(title = tool_name,
         x = if (show_x) "Offset (nt)" else NULL,
         y = "Spearman r") +
    theme_bw(base_size = 7) +
    theme(
      plot.title       = element_text(size = 7, hjust = 0.5, margin = margin(b = 1)),
      axis.title       = element_text(size = 7),
      axis.text        = element_text(size = 6, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25),
      panel.grid       = element_blank(),
      panel.border     = element_rect(linewidth = 0.35, colour = "black"),
      legend.position  = "none",
      plot.margin      = margin(3, 20, 2, 4)
    )
}

# -------------------- build panels --------------------
p1 <- make_panel("CasRx_day5")
p2 <- make_panel("PspCas13b", show_x = TRUE)

# -------------------- shared legend --------------------
p_leg <- ggplot(dt, aes(offset, correlation_r, colour = short, linetype = short)) +
  geom_line(linewidth = 0.35) + 
  scale_colour_manual(values = col_map, name = NULL) +
  scale_linetype_manual(values = lty_map, name = NULL) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_bw(base_size = 6) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.key.width  = unit(0.6, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.text       = element_text(size = 5.5),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(-3, 0, 0, 0),
    legend.background = element_blank(),
    legend.key        = element_blank()
  )
leg <- get_legend(p_leg)

# -------------------- assemble --------------------
# 重新分配高度比例，适配 2 个面板
fig <- ggarrange(p1, p2, leg,
                 ncol = 1, nrow = 3,
                 heights = c(1, 1.08, 0.38),
                 align = "v")

out_pdf <- file.path(OUT_DIR, "Supple_Fig_offset_wet_dry.pdf")

# 整体宽度 2.2 不变，高度更改为 3
ggsave(out_pdf, fig, width = 2.2, height = 3,
       device = cairo_pdf, bg = "white")

# --- Save best offset stats ---
fwrite(best[, .(tool, model, short, offset, correlation_r)], 
       file.path(OUT_DIR, "SF2D_Best_Offset_Stats.csv"))
message("Saved: ", out_pdf)