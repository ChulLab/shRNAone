# partC.02.extract_WL_matrix.R 
# 
# 核心逻辑：
# 1. 锁定 3prime 方向，同时提取 whole 和 3'UTR 区域的数据。
# 2. 强制统一使用基于 whole 区域定义的黄金 U 平台范围（作为严格的物理约束）。
# 3. 扫描所有 (W, L) 组合，计算其在各自区域内的整体表现 (Mean_r) 和最佳表现 (Max_r)。
# 4. 采用自定义加权打分系统，分别输出 whole 和 3'UTR 两个高维性能矩阵。

library(data.table)
library(dplyr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/01.Self_WL_list")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")
REGIONS_TO_EXTRACT <- c("whole", "3'UTR")

# --- 核心约束 1: 黄金 U 平台 (基于 whole 区域严格定义) ---
U_PLATEAU <- list(
  shRNAone   = c(6, 14),
  PspCas13b  = c(9, 23),
  CasRx_day5 = c(8, 14)
)

# --- 核心约束 2: 加权打分系统 ---
WEIGHTS <- c(
  PspCas13b  = 0.50,
  shRNAone   = 0.30,
  CasRx_day5 = 0.20
)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 1. 极速加载并过滤数据 =================
log_msg("Loading grid results and applying Base constraints (3prime, whole & 3'UTR)...")

cols_to_keep <- c("W", "L", "U", "direction", "region", "spearman_r")
dt_list <- lapply(TOOLS, function(t) {
  fpath <- file.path(IN_DIR, paste0(t, "_grid.csv"))
  if (file.exists(fpath)) {
    dt <- fread(fpath, select = cols_to_keep)
    dt$tool <- t
    return(dt)
  }
})
all_dt <- rbindlist(dt_list)

# 仅保留锚定条件的子集：3prime 方向，包含我们需要的两个区域，并且要求正相关
anchor_dt <- all_dt[direction == "3prime" & region %in% REGIONS_TO_EXTRACT & spearman_r > 0]

# ================= 2. 应用 U 平台遮罩 (统一物理约束) =================
log_msg("Applying unified tool-specific U Plateau masks...")

plateau_dt <- data.table()

for (t in TOOLS) {
  u_min <- U_PLATEAU[[t]][1]
  u_max <- U_PLATEAU[[t]][2]
  
  # 对该工具提取指定 U 平台内的数据（此处的 U 范围约束同时作用于 whole 和 3'UTR）
  sub_dt <- anchor_dt[tool == t & U >= u_min & U <= u_max]
  plateau_dt <- rbind(plateau_dt, sub_dt)
}

# ================= 3. 聚合计算 (引入 Region 维度) =================
log_msg("Calculating Max_r and Mean_r grouped by (W, L), tool, and region...")

# 按 W, L, tool 以及 region 分组进行聚合计算
agg_dt <- plateau_dt[, .(
  Max_r = max(spearman_r, na.rm = TRUE),
  Mean_r = mean(spearman_r, na.rm = TRUE)
), by = .(W, L, tool, region)]

# ================= 4. 矩阵重构与独立加权输出 =================
log_msg("Constructing matrices and applying weighting system...")

for (tgt_region in REGIONS_TO_EXTRACT) {
  
  # 提取当前区域的数据
  sub_agg <- agg_dt[region == tgt_region]
  
  # 将长数据转化为宽数据矩阵
  wide_dt <- dcast(sub_agg, W + L ~ tool, value.var = c("Max_r", "Mean_r"))
  
  # 清理 NA 值 (用 0 填充缺失的组合)
  cols_to_fill <- setdiff(names(wide_dt), c("W", "L"))
  for (j in cols_to_fill) {
    if (j %in% names(wide_dt)) {
      set(wide_dt, which(is.na(wide_dt[[j]])), j, 0)
    }
  }
  
  # 安全地应用加权打分系统 (防止某些列因为完全没有数据而缺失)
  calc_weighted_score <- function(dt, metric, weights) {
    score <- rep(0, nrow(dt))
    for (t in names(weights)) {
      col_name <- paste0(metric, "_", t)
      if (col_name %in% names(dt)) {
        score <- score + (dt[[col_name]] * weights[t])
      }
    }
    return(score)
  }
  
  wide_dt[, Weighted_Score_Max := calc_weighted_score(wide_dt, "Max_r", WEIGHTS)]
  wide_dt[, Weighted_Score_Mean := calc_weighted_score(wide_dt, "Mean_r", WEIGHTS)]
  
  # 根据加权平均分降序排序
  setorder(wide_dt, -Weighted_Score_Mean)
  
  # 针对当前区域安全地重命名特定的列，避免文件名包含非法字符
  file_region_name <- gsub("'", "", tgt_region) # 将 3'UTR 转换为 3UTR
  
  # 输出结果
  out_file <- file.path(OUT_DIR, sprintf("00.WL_Performance_Matrix_%s.csv", file_region_name))
  fwrite(wide_dt, out_file)
  
  log_msg(sprintf("Extraction complete for [%s]. Matrix saved to: %s", tgt_region, out_file))
}

log_msg("All matrix extraction tasks completed successfully.")