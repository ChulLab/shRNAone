# partC.04.select_representative_WL.R 
# 
# 核心逻辑：
# 1. 导入02.Block_Stats_Strict_whole_Weighted_Score_Mean.csv
# 2. 选取代表性参数：在 5x5 区块内，精确提取 W 和 L 均为 5 的倍数的坐标作为代表。
# 3. 应用空间非极大值抑制 (Spatial NMS) 算法：
#    - 要求候选区块之间保持足够的物理距离（W 至少相隔 5，L 至少相隔 5）。
# 4. 提取 Top 20 个空间独立的稳健极值点。
# 5. 混入经典对照 (W=80, L=40)。

library(data.table)
library(dplyr)

# ================= 0. 配置与初始化 =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partC/01.Self_WL_list")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partC/01.Self_WL_list")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 目标：从 whole 区域的 Mean 打分中提取最稳健的代表
TARGET_FILE <- "02.Block_Stats_Strict_whole_Weighted_Score_Mean.csv"

# 空间隔离阈值 (微调为 5)
MIN_DIST_W <- 5  
MIN_DIST_L <- 5   
TARGET_COUNT <- 20 # 提取 20 个代表

REF_W <- 80
REF_L <- 40

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 1. 加载并初始化 =================
log_msg(sprintf("Loading strict block stats: %s", TARGET_FILE))
block_stats <- fread(file.path(IN_DIR, TARGET_FILE))

# 取该 5x5 区块内唯一的 5 的倍数
block_stats[, Rep_W := ceiling(Block_W / 5) * 5]
block_stats[, Rep_L := ceiling(Block_L / 5) * 5]

# ================= 2. 空间非极大值抑制 (Spatial NMS) =================
log_msg("Applying Spatial NMS to select independent representative clusters...")

selected_candidates <- data.table()

for (i in 1:nrow(block_stats)) {
  current_row <- block_stats[i]
  
  if (nrow(selected_candidates) == 0) {
    selected_candidates <- rbind(selected_candidates, current_row)
    next
  }
  
  # 计算当前代表点与所有已入选代表点的距离
  is_far_enough <- sapply(1:nrow(selected_candidates), function(j) {
    dist_w <- abs(current_row$Rep_W - selected_candidates$Rep_W[j])
    dist_l <- abs(current_row$Rep_L - selected_candidates$Rep_L[j])
    return(dist_w >= MIN_DIST_W | dist_l >= MIN_DIST_L)
  })
  
  if (all(is_far_enough)) {
    selected_candidates <- rbind(selected_candidates, current_row)
  }
  
  if (nrow(selected_candidates) >= TARGET_COUNT) {
    break
  }
}

# ================= 3. 整理输出清单 =================
candidate_list <- selected_candidates[, .(
  Source = sprintf("Cluster_Rank_%d", 1:.N),
  W = Rep_W,
  L = Rep_L,
  Block_Score = round(Block_Score, 4)
)]

baseline <- data.table(
  Source = "Reference",
  W = REF_W,
  L = REF_L,
  Block_Score = NA_real_
)

final_list <- rbind(candidate_list, baseline)

out_file <- file.path(OUT_DIR, "03.selected_WL_Candidates.csv")
fwrite(final_list, out_file)

log_msg("Selection complete. Candidate list generated:")
print(final_list)
log_msg(sprintf("Saved to: %s", out_file))