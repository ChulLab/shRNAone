#!/usr/bin/env Rscript
# derived from: partD.01.optimal_U_selection.R
# 
# SF1G: Optimal U Selection for W=35, L=20 vs W=80, L=40 across all tools and regions

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("CasRx_day5", "PspCas13b", "shRNAone")
TARGET_REGIONS <- c("whole", "CDS", "3'UTR")
CAND_LABEL <- "W35_L20"
REF_LABEL  <- "W80_L40"

U_PLATEAU <- data.table(
  tool   = rep(TOOLS, each = length(TARGET_REGIONS)),
  region = rep(TARGET_REGIONS, times = length(TOOLS)),
  xmin   = c(8, 8, 8, 9, 9, 9, 6, 6, 6),
  xmax   = c(14, 14, 14, 23, 23, 23, 14, 14, 14),
  ymin   = -Inf,
  ymax   = Inf
)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. 数据载入与过滤 =================
log_msg("Step 1: Loading grid results...")
cols_to_keep <- c("W", "L", "U", "direction", "region", "spearman_r", "spearman_p")
dt_list <- lapply(TOOLS, function(t) {
  fpath <- file.path(IN_DIR, paste0(t, "_grid.csv"))
  if (file.exists(fpath)) {
    dt <- fread(fpath, select = cols_to_keep)
    dt[, tool := t]
    return(dt)
  }
})
all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

plot_dt <- all_dt[
  direction == "3prime" &
  region %in% TARGET_REGIONS &
  ((W == 35 & L == 20) | (W == 80 & L == 40))
]

plot_dt[, WL_Label := sprintf("W%d_L%d", W, L)]
plot_dt[, WL_Label := factor(WL_Label, levels = c(CAND_LABEL, REF_LABEL))]
plot_dt[, tool     := factor(tool,     levels = TOOLS)]
plot_dt[, region   := factor(region,   levels = TARGET_REGIONS)]

U_PLATEAU[, tool   := factor(tool,   levels = TOOLS)]
U_PLATEAU[, region := factor(region, levels = TARGET_REGIONS)]

# ================= 2. 锚定 whole 区域最佳 U 并映射 =================
log_msg("Step 2: Extracting optimal U...")
peak_u_info <- plot_dt[region == "whole",
                       .(Peak_U = U[which.max(spearman_r)]),
                       by = .(tool, WL_Label)]

peak_points <- merge(plot_dt, peak_u_info, by = c("tool", "WL_Label"))
peak_points <- peak_points[U == Peak_U]

peak_points[, label_txt := sprintf("U=%d, r=%.2f", Peak_U, spearman_r)]
peak_points[, v_align   := ifelse(WL_Label == CAND_LABEL, -0.6, 1.6)]

# ================= 3. 独立构建三个工具函数以彻底杜绝变量作用域污染 =================
make_tool_plot <- function(target_tool_name, show_y_axis_title = FALSE) {
  t_data <- plot_dt[tool == target_tool_name]
  t_peak <- peak_points[tool == target_tool_name]
  t_plat <- U_PLATEAU[tool == target_tool_name]   # now safe – U_PLATEAU is data.table

  p <- ggplot(t_data, aes(x = U, y = spearman_r, color = WL_Label, linetype = WL_Label)) +
    geom_rect(data = t_plat,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#D28E72", alpha = 0.15, color = NA, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(data = t_peak, aes(x = Peak_U, y = spearman_r),
               size = 0.8, stroke = 0.5) +
    geom_text(data = t_peak,
              aes(x = Peak_U, y = spearman_r, label = label_txt,
                  color = WL_Label, vjust = v_align),
              size = 1.8, fontface = "bold", show.legend = FALSE) +
    facet_grid(. ~ region) +
    coord_cartesian(xlim = c(1, 25)) +
    scale_x_continuous(breaks = c(5, 15, 25)) +
    scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.15, 0.4))) +
    scale_color_manual(values = setNames(c("#E74C3C", "gray40"),
                                         c(CAND_LABEL, REF_LABEL))) +
    scale_linetype_manual(values = setNames(c("solid", "dashed"),
                                            c(CAND_LABEL, REF_LABEL))) +
    labs(title = target_tool_name, x = "Target Length (U)") +
    theme_bw(base_size = 6, base_family = "Arial") +
    theme(
      plot.title      = element_text(face = "bold", size = 7, hjust = 0.5,
                                     margin = margin(b = 2)),
      strip.background = element_rect(fill = "grey90", color = "black",
                                      linewidth = 0.2),
      strip.text      = element_text(face = "bold", size = 5,
                                     margin = margin(t = 1, b = 1)),
      panel.spacing   = unit(0, "lines"),          # regions of the same tool glued together
      axis.title      = element_text(face = "bold", size = 6),
      axis.text       = element_text(color = "black", size = 5),
      axis.text.x     = element_text(margin = margin(t = 1)),
      panel.grid.minor = element_blank(),
      panel.border    = element_rect(color = "black", linewidth = 0.2),
      plot.margin     = margin(1, 3, 1, 3)
    )

  if (show_y_axis_title) {
    p <- p + ylab("Spearman (r)")
  } else {
    p <- p + ylab(NULL)
  }
  return(p)
}

log_msg("Step 3: Rendering individual tool panels...")
p_cas <- make_tool_plot("CasRx_day5", show_y_axis_title = TRUE)
p_psp <- make_tool_plot("PspCas13b",  show_y_axis_title = FALSE)
p_shr <- make_tool_plot("shRNAone",   show_y_axis_title = FALSE)

# three independent plots → each keeps its own y-scale;
# ggarrange places them side-by-side with a shared legend
fig_G <- ggarrange(
  p_cas, p_psp, p_shr,
  ncol = 3, nrow = 1,
  align = "h",
  common.legend = TRUE,
  legend = "bottom"
)

out_pdf <- file.path(OUT_DIR, "Supple_FigG_W35L20_vs_W80L40.pdf")
ggsave(out_pdf, plot = fig_G, width = 7.5, height = 1.5, device = cairo_pdf)

# --- Save peak U comparison stats ---
peak_out <- merge(peak_points, U_PLATEAU, by = c("tool", "region"))
peak_out <- peak_out[, .(tool, WL_Label, region, Peak_U, Peak_r = spearman_r, Plateau_Min = xmin, Plateau_Max = xmax)]
fwrite(peak_out, file.path(OUT_DIR, "SF1G_Peak_U_Comparison.csv"))
log_msg(sprintf("Successfully generated: %s", out_pdf))