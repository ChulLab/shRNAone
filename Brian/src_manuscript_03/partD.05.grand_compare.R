#!/usr/bin/env Rscript
# partD.05.grand_compare.R 
#

library(ggplot2)
library(data.table)
library(dplyr)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
PART_A_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
PART_B_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partD/02.best_combination")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COR_METHOD <- "spearman"

# 定义 Part D 最终敲定的“双参数标准”
TARGET_PARAMS <- list(
  "shRNAone"   = list(c(35, 20, 7),  c(80, 40, 7)),
  "PspCas13b"  = list(c(35, 20, 11), c(80, 40, 17)),
  "CasRx_day5" = list(c(35, 20, 11), c(80, 40, 9))
)

get_sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# ================= 1. 数据合并与重构 =================
log_msg("Step 1: Loading Part A external models baseline...")
part_a_file <- file.path(PART_A_DIR, "05.self_other_compare.csv")
if (!file.exists(part_a_file)) stop("Part A results not found!")
df_external <- fread(part_a_file)

log_msg("Step 2: Extracting targeted Self-Models from Part B grids...")
self_models_list <- list()

for (tool_name in names(TARGET_PARAMS)) {
  grid_file <- file.path(PART_B_DIR, paste0(tool_name, "_grid.csv"))
  if (!file.exists(grid_file)) next
  
  df_grid <- fread(grid_file)
  params <- TARGET_PARAMS[[tool_name]]
  
  for (p in params) {
    w_val <- p[1]; l_val <- p[2]; u_val <- p[3]
    
    # 精准提取：指定 W, L, U 且方向为 3prime 的行
    sub_df <- df_grid[W == w_val & L == l_val & U == u_val & direction == "3prime"]
    
    # 转换为与 Part A 完全一致的数据结构
    for (i in 1:nrow(sub_df)) {
      raw_region <- sub_df$region[i]
      # 抹平 Part B (whole) 与 Part A (Whole) 的命名差异
      clean_region <- ifelse(raw_region == "whole", "Whole", raw_region)
      
      # 过滤掉非核心区域（如非 CDS, 3'UTR, Whole 的片段）
      if (!clean_region %in% c("Whole", "CDS", "3'UTR")) next
      
      self_models_list[[length(self_models_list) + 1]] <- data.table(
        Tool = tool_name,
        Region = clean_region,
        Method = sprintf("Self_W%d_L%d", w_val, l_val),
        Correlation = sub_df$spearman_r[i],
        P_value = sub_df$spearman_p[i],
        Sig_Label = get_sig_label(sub_df$spearman_p[i]),
        N = NA # 网格数据未存 N，此处设为 NA 不影响绘图
      )
    }
  }
}

df_self <- rbindlist(self_models_list)
df_grand <- rbind(df_external, df_self, fill = TRUE)

# 保存大横评统计表
out_csv <- file.path(OUT_DIR, "01.Grand_Performance_Compare.csv")
fwrite(df_grand, out_csv)
log_msg(sprintf("Grand comparison stats saved to: %s", out_csv))

# ================= 2. 颜色映射与绘图 =================
log_msg("Step 3: Generating faceted bar plots...")

df_grand <- df_grand[!is.na(Correlation)]
df_grand$Tool <- factor(df_grand$Tool, levels = names(TARGET_PARAMS))
df_grand$Region <- factor(df_grand$Region, levels = c("Whole", "CDS", "3'UTR"))

# 设定 Method 顺序，让 Self_Models 排在最前面
methods <- unique(df_grand$Method)
self_methods <- grep("Self_", methods, value = TRUE)
other_methods <- setdiff(methods, self_methods)
df_grand$Method <- factor(df_grand$Method, levels = c(sort(self_methods), sort(other_methods)))

# 动态计算显著性星号的 Y 轴位置
df_grand$label_y <- ifelse(df_grand$Correlation >= 0, 
                           df_grand$Correlation + 0.04, 
                           df_grand$Correlation - 0.04)

# 添加高亮颜色分类组
df_grand <- df_grand %>%
  mutate(Color_Group = case_when(
    grepl("W35_L20", Method) ~ "Optimal (W35L20)",
    grepl("W80_L40", Method) ~ "Baseline (W80L40)",
    TRUE ~ "External Models"
  ))

# 匹配 Part A 的绘图风格并加入高亮配色
p_bar <- ggplot(df_grand, aes(x = Method, y = Correlation, fill = Color_Group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = Sig_Label, y = label_y), size = 2.0, 
            vjust = ifelse(df_grand$Correlation >= 0, 0, 1)) +
  # 统一调色盘：外部灰，W80蓝，W35红
  scale_fill_manual(values = c("Optimal (W35L20)" = "#E41A1C", 
                               "Baseline (W80L40)" = "#377EB8", 
                               "External Models" = "gray70")) +
  facet_grid(Region ~ Tool, scales = "free_x", space = "free_x") + 
  scale_y_continuous(limits = c(-0.1, 0.7), breaks = c(0, 0.25, 0.5)) +
  labs(
    title = sprintf("Grand Performance Comparison (%s r)", tools::toTitleCase(COR_METHOD)),
    y = sprintf("%s r", tools::toTitleCase(COR_METHOD)), x = NULL, fill = "Model Type"
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
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

pdf_path <- file.path(OUT_DIR, "02.Grand_Performance_Compare.pdf")
# 稍微增加高度以容纳底部的 Legend
pdf(pdf_path, width = 5.5, height = 3.2) 
print(p_bar)
dev.off()
log_msg(sprintf("Grand multi-region bar plot saved successfully to: %s", pdf_path))