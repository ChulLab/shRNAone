#!/usr/bin/env Rscript
# derived from: partE.02.clusters_acc_compare.R
#
# Fig6D: Accessibility comparison for TTR and PCSK9 clusters (window=21)

library(data.table)
library(ggplot2)
library(pROC)
library(patchwork)

# ================= 0. Configuration =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_CLUSTER <- file.path(PROJECT, "result/08.manuscript_03/partE/01.spatial_clusters/01.Spatial_Clusters_Details.csv")

OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/main")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MAT_PATHS <- list(
  TTR_W35 = file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix/TTR_shRNAone_TTR_based_W35_L20_U7.csv"),
  TTR_W80 = file.path(PROJECT, "result/08.manuscript_03/partD/00.TTR_pred/03.matrix/TTR_shRNAone_RNAxs_based_W80_L40_U7.csv"),
  PCSK9_W35 = file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix/PCSK9_shRNAone_TTR_based_W35_L20_U7.csv"),
  PCSK9_W80 = file.path(PROJECT, "result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix/PCSK9_shRNAone_RNAxs_based_W80_L40_U7.csv")
)

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. Load Data =================
log_msg("Loading cluster labels and accessibility matrices...")
cluster_dt <- fread(IN_CLUSTER)

load_acc <- function(path, param_name, gene_name) {
  dt <- fread(path)[, .(start_pos, Accessibility)]
  dt[, Parameter := param_name]
  dt[, Gene := gene_name]
  return(dt)
}

acc_all <- rbindlist(list(
  load_acc(MAT_PATHS$TTR_W35, "W35_L20_U7", "TTR"),
  load_acc(MAT_PATHS$TTR_W80, "W80_L40_U7", "TTR"),
  load_acc(MAT_PATHS$PCSK9_W35, "W35_L20_U7", "PCSK9"),
  load_acc(MAT_PATHS$PCSK9_W80, "W80_L40_U7", "PCSK9")
))

merged_dt <- merge(acc_all, cluster_dt[, .(Gene, start_pos, Region_Type, Window_NT)], 
                   by = c("Gene", "start_pos"), all.x = FALSE, allow.cartesian = TRUE)

# ================= 2. Data Processing =================
log_msg("Preparing data, calculating ROC thresholds and significance brackets...")

window_dt <- merged_dt[Window_NT == 21]

# Unify Y-axis to log scale
window_dt$Accessibility_log <- log10(window_dt$Accessibility + 0.01)

# 客观标注（与 partE.01 的计算逻辑完全对应）
# Optimal   → Rolling_log2FC > 99.5% CI  → >CI99.5
# Background → within CI
# Poor      → Rolling_log2FC < 0.5% CI   → <CI0.5
window_dt[, Region_Type_Clean := factor(
  fcase(
    Region_Type == "Optimal",    ">CI99.5",
    Region_Type == "Background", "within CI",
    Region_Type == "Poor",       "<CI0.5"
  ),
  levels = c(">CI99.5", "within CI", "<CI0.5")
)]

window_dt[Parameter == "W35_L20_U7", Parameter_Full := "W35 L20 U7"]
window_dt[Parameter == "W80_L40_U7", Parameter_Full := "W80 L40 U7"]
window_dt$Parameter_Full <- factor(window_dt$Parameter_Full, levels = c("W35 L20 U7", "W80 L40 U7"))
window_dt$Gene <- factor(window_dt$Gene, levels = c("TTR", "PCSK9"))

# n 按 Parameter + Gene + Clean 计算
window_dt[, n := .N, by = .(Parameter_Full, Gene, Region_Type_Clean)]


# ================= 3. Plotting Functions =================
build_plot_data <- function(sub_dt) {
  sig_list <- list()
  thresh_list <- list()
  blank_list <- list()
  
  param_max <- max(sub_dt$Accessibility_log, na.rm = TRUE)
  param_min <- min(sub_dt$Accessibility_log, na.rm = TRUE)
  param_span <- param_max - param_min
  if (param_span == 0) param_span <- 1
  
  y_limit <- param_max + param_span * 0.55
  tip_len <- param_span * 0.03
  
  for (g in levels(sub_dt$Gene)) {
    sub_g <- sub_dt[Gene == g]
    if (nrow(sub_g) == 0) next
    
    gene_max <- max(sub_g$Accessibility_log, na.rm = TRUE)
    
    y_pos_1 <- gene_max + param_span * 0.15
    y_pos_2 <- gene_max + param_span * 0.30
    y_pos_3 <- gene_max + param_span * 0.45
    
    v_opt  <- sub_g[Region_Type_Clean == ">CI99.5", Accessibility_log]
    v_bg   <- sub_g[Region_Type_Clean == "within CI", Accessibility_log]
    v_poor <- sub_g[Region_Type_Clean == "<CI0.5", Accessibility_log]
    
    get_sig <- function(v1, v2) {
      if (length(v1) >= 2 && length(v2) >= 2) {
        pval <- wilcox.test(v1, v2)$p.value
        return(ifelse(pval < 0.0001, "****",
                      ifelse(pval < 0.001, "***",
                             ifelse(pval < 0.01, "**",
                                    ifelse(pval < 0.05, "*", "ns")))))
      }
      return(NA)
    }
    
    s1 <- get_sig(v_opt, v_bg)
    if (!is.na(s1)) sig_list[[length(sig_list) + 1]] <- data.table(Gene = g, x1 = 1, x2 = 2, y = y_pos_1, label = s1, tip = tip_len)
    
    s2 <- get_sig(v_bg, v_poor)
    if (!is.na(s2)) sig_list[[length(sig_list) + 1]] <- data.table(Gene = g, x1 = 2, x2 = 3, y = y_pos_2, label = s2, tip = tip_len)
    
    s3 <- get_sig(v_opt, v_poor)
    if (!is.na(s3)) sig_list[[length(sig_list) + 1]] <- data.table(Gene = g, x1 = 1, x2 = 3, y = y_pos_3, label = s3, tip = tip_len)
    
    # ROC Cutoff
    response <- ifelse(sub_g$Region_Type_Clean == ">CI99.5", 1, 0)
    predictor <- sub_g$Accessibility
    if (length(unique(response)) == 2) {
      roc_obj <- suppressMessages(roc(response, predictor, direction = "<", quiet = TRUE))
      auc_val <- as.numeric(auc(roc_obj))
      best_coords <- coords(roc_obj, "best", ret = c("threshold"), best.method = "youden")
      best_thresh <- best_coords$threshold[1]
      
      thresh_list[[length(thresh_list) + 1]] <- data.table(
        Gene = g,
        Threshold_Log = log10(best_thresh + 0.01),
        Label = sprintf("Cutoff: %.2f (AUC: %.2f)", best_thresh, auc_val)
      )
    }
    
    blank_list[[length(blank_list) + 1]] <- data.table(
      Gene = g,
      Y_Lim = y_limit,
      Region_Type_Clean = ">CI99.5"
    )
  }
  
  return(list(sig = rbindlist(sig_list), thresh = rbindlist(thresh_list), blank = rbindlist(blank_list)))
}

plot_param <- function(param_str, is_left = TRUE) {
  sub_dt <- window_dt[Parameter_Full == param_str]
  res <- build_plot_data(sub_dt)
  
  p <- ggplot(sub_dt, aes(x = Region_Type_Clean, y = Accessibility_log)) +
    geom_blank(data = res$blank, aes(y = Y_Lim)) +
    geom_boxplot(outlier.shape = NA, fill = "gray80", color = "black", linewidth = 0.2, alpha = 0.85) +
    geom_hline(data = res$thresh, aes(yintercept = Threshold_Log), linetype = "dashed", color = "#C0392B", linewidth = 0.4) +
    geom_text(data = res$thresh, aes(x = 2, y = Threshold_Log, label = Label), 
              hjust = 0.5, vjust = -0.5, size = 6.5 / .pt, family = "Arial", color = "#C0392B", inherit.aes = FALSE) +
    geom_segment(data = res$sig, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = res$sig, aes(x = x1, xend = x1, y = y, yend = y - tip), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = res$sig, aes(x = x2, xend = x2, y = y, yend = y - tip), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = res$sig, aes(x = (x1 + x2)/2, y = y + (tip * 0.8), label = label), 
              inherit.aes = FALSE, size = 8 / .pt, family = "Arial", vjust = 0) +
    
    facet_grid(. ~ Gene) +
    
    # 同一行显示名称 + (n=*)，不再分行
    scale_x_discrete(labels = function(x) {
      sapply(x, function(lv) {
        n_val <- sub_dt[Region_Type_Clean == lv, n][1]
        sprintf("%s (n=%d)", lv, n_val)
      })
    }) +
    
    labs(x = NULL, y = expression(log[10](Acc+0.01))) +
    ggtitle(param_str) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(b = 4)),
      strip.background = element_rect(fill = "#E5E7E9", color = "black", linewidth = 0.3),
      strip.text = element_text(size = 8, color = "black"),
      axis.text.x = element_text(color = "black", size = 6.5, angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black", size = 7),
      plot.margin = margin(5, 2, 2, 2, "pt"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.3),
      panel.spacing = unit(0.2, "lines")
    )
  
  if (!is_left) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  return(p)
}

log_msg("Generating final plots with Patchwork layout...")

p1 <- plot_param("W35 L20 U7", is_left = TRUE)
p2 <- plot_param("W80 L40 U7", is_left = FALSE)

final_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1)) & theme(legend.position = "none")

out_pdf <- file.path(OUT_DIR, "Task4_Clusters_Acc_Window21.pdf")
ggsave(out_pdf, plot = final_plot, width = 4.5, height = 2, device = cairo_pdf)

# --- Save cluster comparison stats ---
cluster_stats <- list()
for (param_str in c("W35 L20 U7", "W80 L40 U7")) {
  sub_dt <- window_dt[Parameter_Full == param_str]
  # n per group
  for (g in levels(sub_dt$Gene)) {
    sub_g <- sub_dt[Gene == g]
    for (ct in c(">CI99.5", "within CI", "<CI0.5")) {
      vals <- sub_g[Region_Type_Clean == ct, Accessibility]
      cluster_stats[[length(cluster_stats) + 1]] <- data.table(
        Parameter = param_str, Gene = g, Region_Type = ct,
        N = length(vals), Median_Acc = median(vals, na.rm = TRUE),
        IQR_Acc = IQR(vals, na.rm = TRUE)
      )
    }
    # ROC
    response <- ifelse(sub_g$Region_Type_Clean == ">CI99.5", 1, 0)
    predictor <- sub_g$Accessibility
    if (length(unique(response)) == 2) {
      roc_obj <- suppressMessages(roc(response, predictor, direction = "<", quiet = TRUE))
      auc_val <- as.numeric(auc(roc_obj))
      best_coords <- coords(roc_obj, "best", ret = c("threshold"), best.method = "youden")
      best_thresh <- best_coords$threshold[1]
      cluster_stats[[length(cluster_stats) + 1]] <- data.table(
        Parameter = param_str, Gene = g, Region_Type = "ROC",
        N = length(response), Median_Acc = auc_val, IQR_Acc = best_thresh
      )
    }
    # Wilcoxon pairwise
    for (pair in list(c(">CI99.5", "within CI"), c("within CI", "<CI0.5"), c(">CI99.5", "<CI0.5"))) {
      v1 <- sub_g[Region_Type_Clean == pair[1], Accessibility]
      v2 <- sub_g[Region_Type_Clean == pair[2], Accessibility]
      if (length(v1) >= 2 && length(v2) >= 2) {
        wt <- wilcox.test(v1, v2)
        cluster_stats[[length(cluster_stats) + 1]] <- data.table(
          Parameter = param_str, Gene = g,
          Region_Type = paste(pair[1], "vs", pair[2]),
          N = length(v1) + length(v2),
          Median_Acc = NA_real_, IQR_Acc = wt$p.value
        )
      }
    }
  }
}
fwrite(rbindlist(cluster_stats, fill = TRUE), file.path(OUT_DIR, "Task4_Cluster_Comparison_Stats.csv"))
log_msg("Saved Task4 CSV statistics.")

log_msg(sprintf("Saved Task 4 PDF: %s", out_pdf))