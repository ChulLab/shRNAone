#!/usr/bin/env Rscript
# derived from: partC.03.WL_landscape.R
#
# SF1E: WL Landscape with Candidate Points Overlay

library(ggplot2)
library(data.table)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
WORK_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/01.Self_WL_list")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 核心锚点
REF_W <- 80
REF_L <- 40
OPT_W <- 35
OPT_L <- 20

# 经典热力图配色
heatmap_colors <- c("#313695", "#74ADD1", "#E0F3F8", "#FDAE61", "#F46D43", "#A50026")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. 加载数据 =================
log_msg("Loading performance matrix (whole region)...")
matrix_file <- file.path(WORK_DIR, "00.WL_Performance_Matrix_whole.csv")
if (!file.exists(matrix_file)) stop("Matrix file not found. Please run partC.02 first.")

matrix_dt <- fread(matrix_file)
# 强制物理约束：W >= L, 并限定数据加载范围
plot_data <- matrix_dt[W >= L & W <= 100 & W >= 20 & L <= 80 & L >= 10]


# 【新增】：加载代表点候选清单，用于在图上标出区块和红点
log_msg("Loading selected representative candidates...")
ranking_file <- file.path(PROJECT, "result/08.manuscript_03/partC/03.GT_siRNA_acc/03.WL_selection/01.WL_Final_Ranking.csv")
if (!file.exists(ranking_file)) stop("Final Ranking file not found.")

candidates_dt <- fread(ranking_file)
# 仅提取通过 NMS 算法筛选出的候选点 (排除 Reference 点，因为下面会单独用叉号标记)
spatial_points <- candidates_dt[Source != "Reference"]


# ================= 2. 提取核心锚点数据构造 Dashboard =================
log_msg("Extracting Reference and Optimum point data...")

# 提取指定 W, L 的各工具跑分
extract_scores <- function(dt, w, l) {
  row <- dt[W == w & L == l]
  if(nrow(row) > 0) {
    return(list(
      psp = row$Mean_r_PspCas13b,
      shr = row$Mean_r_shRNAone,
      cas = row$Mean_r_CasRx_day5,
      score = row$Weighted_Score_Mean
    ))
  } else {
    return(list(psp=NA, shr=NA, cas=NA, score=NA))
  }
}

ref_data <- extract_scores(matrix_dt, REF_W, REF_L)
opt_data <- extract_scores(matrix_dt, OPT_W, OPT_L)

# 构造左上角详细 Dashboard 文本 (取消前缀，保持最简洁)
dashboard_txt <- sprintf(
  "× Ref (%d, %d)  |  Weighted Score: %.3f\n    Psp: %.2f  |  shRNA: %.2f  |  CasRx: %.2f\n\n+ Opt (%d, %d)  |  Weighted Score: %.3f\n    Psp: %.2f  |  shRNA: %.2f  |  CasRx: %.2f",
  REF_W, REF_L, ref_data$score, ref_data$psp, ref_data$shr, ref_data$cas,
  OPT_W, OPT_L, opt_data$score, opt_data$psp, opt_data$shr, opt_data$cas
)

# ================= 3. 高密度极简出图 =================
log_msg("Rendering updated Figure E with spatial points overlay...")

p_land <- ggplot(plot_data, aes(x = W, y = L)) +
  # 1. 底层：栅格与等高线
  geom_raster(aes(fill = Weighted_Score_Mean), interpolate = TRUE) +
  geom_contour(aes(z = Weighted_Score_Mean), color = "white", alpha = 0.4, bins = 12, linewidth = 0.3) +
  
  # 2. 【新增】中层：候选区块标识 (极细 5x5 方框)
  geom_rect(data = spatial_points,
            aes(xmin = W - 2.5, xmax = W + 2.5, ymin = L - 2.5, ymax = L + 2.5),
            fill = NA, color = "gray20", linewidth = 0.25, alpha = 0.7, inherit.aes = FALSE) +
            
  # 3. 【新增】中层：候选点中心红点
  geom_point(data = spatial_points, aes(x = W, y = L),
             color = "#E74C3C", size = 1, shape = 16, inherit.aes = FALSE) +
  
  # 4. 顶层：Opt 最终优选点 (+) 和 Ref 基准点 (×) 覆盖在红点之上确保清晰
  geom_point(data = data.frame(W = OPT_W, L = OPT_L), aes(x = W, y = L), 
             color = "black", shape = 3, size = 2.2, stroke = 1, inherit.aes = FALSE) +
  geom_point(data = data.frame(W = REF_W, L = REF_L), aes(x = W, y = L), 
             color = "black", shape = 4, size = 2.2, stroke = 1, inherit.aes = FALSE) +
  
  # 左上角信息面板
  geom_text(data = data.frame(W = 22, L = 78, label = dashboard_txt),
            aes(x = W, y = L, label = label),
            color = "black", hjust = 0, vjust = 1, size = 2.2, fontface = "bold", family = "Arial", inherit.aes = FALSE) +
  
  scale_fill_gradientn(colors = heatmap_colors, name = "Weighted\nMean Score") +
  # W 轴 20-100，L 轴 10-80 截断
  coord_cartesian(xlim = c(20, 100), ylim = c(10, 80)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "WL combination weighted score landscape",
    x = "Macro Window Size (W)",
    y = "Max Folding Span (L)"
  ) +
  theme_bw(base_size = 8, base_family = "Arial") + 
  theme(
    plot.title = element_text(face = "bold", size = 9, margin = margin(b = 4)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 5.5),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

out_pdf <- file.path(OUT_DIR, "Supple_FigE_WL_Landscape_Compact.pdf")
ggsave(out_pdf, plot = p_land, width = 3.5, height = 2.8, device = cairo_pdf)

# --- Save Ref/Opt scores and spatial candidates ---
ref_opt_dt <- data.table(
  Point = c("Ref_W80_L40", "Opt_W35_L20"),
  W = c(REF_W, OPT_W), L = c(REF_L, OPT_L),
  PspCas13b = c(ref_data$psp, opt_data$psp),
  shRNAone  = c(ref_data$shr, opt_data$shr),
  CasRx     = c(ref_data$cas, opt_data$cas),
  Weighted_Score = c(ref_data$score, opt_data$score)
)
fwrite(ref_opt_dt, file.path(OUT_DIR, "SF1E_Ref_Opt_Scores.csv"))
fwrite(spatial_points[, .(W = W, L = L, Block_Score)], 
       file.path(OUT_DIR, "SF1E_Spatial_Candidates.csv"))
log_msg(sprintf("Figure E successfully rendered and saved to: %s", out_pdf))