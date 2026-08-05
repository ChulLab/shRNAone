#!/usr/bin/env Rscript
# partD.02.best_combination_display.R
#

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pROC)

# ================= 0. 核心配置区 =================
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

PARAMS <- list(
  TTR_based = list(
    shRNAone   = c(W=35, L=20, U=7),
    PspCas13b  = c(W=35, L=20, U=11),
    CasRx_day5 = c(W=35, L=20, U=11)
  ),
  RNAxs_based = list(
    shRNAone   = c(W=80, L=40, U=7),
    PspCas13b  = c(W=80, L=40, U=17),
    CasRx_day5 = c(W=80, L=40, U=9)
  )
)

WET_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partA/01.clean_wet"
PCSK9_FASTA <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.PCSK9_pred/01.fasta/NM_174936.3.fasta"
PCSK9_LUNP  <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.PCSK9_pred/02.lunp_files"
PCSK9_MAT   <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix"
TTR_FASTA <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.TTR_pred/01.fasta/NM_000371.4.fasta"
TTR_LUNP  <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.TTR_pred/02.lunp_files"
TTR_MAT   <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.TTR_pred/03.matrix"
SIRNA_EVAL_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partC/03.GT_siRNA_acc/02.siRNA_matrix"
OUT_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/02.best_combination"

invisible(lapply(c(PCSK9_LUNP, PCSK9_MAT, TTR_LUNP, TTR_MAT, OUT_DIR), dir.create, showWarnings = FALSE, recursive = TRUE))

REGION_COLORS <- c("Whole" = "gray70", "5'UTR" = "#3498DB", "CDS" = "#21908C", "3'UTR" = "#FDE725")
MINI_MARGIN <- margin(2, 2, 2, 2, "pt")

# ================= 1. 核心数据与绘图构建函数 =================

run_rnaplfold <- function(fasta_path, W, L, out_dir) {
  lines <- readLines(fasta_path, n = 1)
  seq_name <- gsub("^>| .*$", "", lines[1])
  target_lunp_name <- file.path(out_dir, sprintf("lunp_W%03d_L%03d.txt", W, L))
  if (!file.exists(target_lunp_name)) {
    cwd <- getwd()
    setwd(out_dir) 
    cmd <- sprintf("RNAplfold -W %d -L %d -u 30 < %s", W, L, fasta_path)
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    file.rename(sprintf("%s_lunp", seq_name), target_lunp_name)
    setwd(cwd)
  }
  return(target_lunp_name)
}

read_wet_data <- function(tool_name, gene_name="TTR") {
  if (gene_name == "TTR") {
    file_path <- file.path(WET_DIR, sprintf("TTR_%s_clean.csv", tool_name))
  } else {
    file_path <- file.path(WET_DIR, "PCSK9_shRNAone_3nt_clean.csv") 
  }
  if(!file.exists(file_path)) stop(sprintf("Wet data not found: %s", file_path))
  dt <- fread(file_path)
  setnames(dt, "mean_log2FC", "log2FC")
  return(dt)
}

build_3prime_matrix <- function(lunp_file, wet_dt, target_U) {
  lunp_dt <- fread(lunp_file, skip = 2, fill = TRUE, na.strings = "NA")
  setnames(lunp_dt, 1, "end_pos") 
  acc_col <- target_U + 1 
  if (ncol(lunp_dt) < acc_col) return(NULL)
  
  acc_dt <- lunp_dt[, c(1, acc_col), with = FALSE]
  setnames(acc_dt, 2, "Accessibility")
  merged <- merge(wet_dt, acc_dt, by = "end_pos")
  merged <- merged[!is.na(Accessibility) & !is.na(log2FC)]
  return(merged)
}

get_track_plot <- function(df, title_str) {
  df <- df[order(df$start_pos), ] 
  df$region_col <- ifelse(df$region %in% names(REGION_COLORS), REGION_COLORS[df$region], "#999999") 
  df$acc_log <- log10(df$Accessibility + 0.01) 
  df$acc_log <- ifelse(df$acc_log > 0, 0, df$acc_log) 
  
  y_min <- min(df$log2FC, na.rm = TRUE) 
  y_max <- max(df$log2FC, na.rm = TRUE) 
  if (y_min == y_max) { y_min <- y_min - 1; y_max <- y_max + 1 } 
  y_range <- y_max - y_min 
  df$acc_mapped <- y_min + ((df$acc_log + 2) / 2) * y_range 
  
  ggplot(df, aes(x = start_pos)) +
    geom_ribbon(aes(ymin = y_min, ymax = acc_mapped), fill = "#CCCCCC", alpha = 0.5) + 
    geom_line(aes(y = log2FC, color = region_col, group = 1), linewidth = 0.3) + 
    scale_color_identity() + 
    scale_y_continuous(
      name = expression(log[2][FC]),
      sec.axis = sec_axis(~ (. - y_min) / y_range * 2 - 2, name = expression(log[10](Acc+0.01))) 
    ) +
    labs(title = title_str, x = "Position") +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", margin = margin(b = 2)),
      axis.title.y = element_text(size = 5),
      axis.title.y.right = element_text(size = 5),
      axis.text = element_text(size = 5),
      plot.margin = MINI_MARGIN,
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    )
}

get_bar_plot <- function(df, regions_to_plot) {
  compute_corr <- function(sub_df, reg_filter = NULL) {
    if (!is.null(reg_filter)) sub_df <- sub_df[sub_df$region == reg_filter, ]
    if (nrow(sub_df) < 3) return(NULL)
    ct <- cor.test(sub_df$Accessibility, sub_df$log2FC, method = "spearman", exact = FALSE) 
    list(r = ct$estimate, p_sig = ifelse(ct$p.value < 0.001, "***", ifelse(ct$p.value < 0.01, "**", ifelse(ct$p.value < 0.05, "*", "ns"))))
  }
  
  bar_data_list <- list()
  for (reg in regions_to_plot) {
    rc <- compute_corr(df, if(reg == "Whole") NULL else reg)
    if (!is.null(rc)) bar_data_list[[reg]] <- data.frame(Region = reg, r = rc$r, p_sig = rc$p_sig, stringsAsFactors = FALSE)
  }
  df_bar <- do.call(rbind, bar_data_list)
  df_bar$Region <- factor(df_bar$Region, levels = regions_to_plot)
  
  df_bar$label_text <- sprintf("%.2f%s", df_bar$r, df_bar$p_sig)
  df_bar$label_y <- ifelse(df_bar$r >= 0, df_bar$r + 0.06, df_bar$r - 0.06) 
  
  ggplot(df_bar, aes(x = Region, y = r, fill = Region)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.7) + 
    geom_text(aes(label = label_text, y = label_y), size = 1.8, vjust = ifelse(df_bar$r >= 0, 0, 1)) + 
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) + 
    scale_fill_manual(values = REGION_COLORS) + 
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
    labs(x = NULL, y = "Spearman r") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 5),
      axis.text.y = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      plot.margin = MINI_MARGIN,
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3), 
      legend.position = "none" 
    ) +
    coord_cartesian(clip = "off")
}

get_box_plot <- function(dt, regions_to_plot) {
  get_extremes <- function(sub_dt, reg_name) {
    if(nrow(sub_dt) < 10) return(NULL)
    q20 <- quantile(sub_dt$log2FC, 0.2, na.rm = TRUE)
    q80 <- quantile(sub_dt$log2FC, 0.8, na.rm = TRUE)
    
    top <- sub_dt[log2FC >= q80]
    bot <- sub_dt[log2FC <= q20]
    
    if(nrow(top)>0) top$Group <- "Top 20%"
    if(nrow(bot)>0) bot$Group <- "Bot 20%"
    
    res <- rbind(top, bot)
    if(nrow(res)>0) res$Plot_Region <- reg_name
    return(res)
  }
  
  box_list <- list()
  if ("Whole" %in% regions_to_plot) box_list[["Whole"]] <- get_extremes(dt, "Whole")
  for (r in unique(dt$region)) {
    if (r %in% regions_to_plot) box_list[[r]] <- get_extremes(dt[region == r], r)
  }
  
  box_dt <- rbindlist(box_list)
  box_dt$Plot_Region <- factor(box_dt$Plot_Region, levels = regions_to_plot)
  box_dt$Group <- factor(box_dt$Group, levels = c("Bot 20%", "Top 20%"))
  
  # 自动计算 upper whisker 以自适应 Y 轴 (避免因 outlier.shape=NA 而导致的坐标轴畸变)
  upper_lim <- max(sapply(split(box_dt, list(box_dt$Plot_Region, box_dt$Group)), function(sub) {
    if(nrow(sub) < 5) return(0)
    boxplot.stats(sub$Accessibility)$stats[5]
  }), na.rm = TRUE)
  if(is.na(upper_lim) || upper_lim == 0) upper_lim <- max(box_dt$Accessibility, na.rm = TRUE)
  
  # 留出顶部空间专门用于画显著性横线
  y_max <- upper_lim * 1.35 
  y_bracket <- upper_lim * 1.15
  y_text <- upper_lim * 1.18
  tip_len <- upper_lim * 0.03
  
  # 纯手工计算 Wilcoxon p-value，彻底规避 stat_compare_means 的绘图缺陷
  sig_list <- list()
  for (reg in regions_to_plot) {
    sub_bot <- box_dt[Plot_Region == reg & Group == "Bot 20%", Accessibility]
    sub_top <- box_dt[Plot_Region == reg & Group == "Top 20%", Accessibility]
    if(length(sub_bot) >= 3 && length(sub_top) >= 3) {
      p_val <- wilcox.test(sub_bot, sub_top)$p.value
      p_sig <- ifelse(p_val < 0.0001, "****", ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns"))))
      sig_list[[reg]] <- data.table(Plot_Region = reg, p_sig = p_sig)
    }
  }
  sig_dt <- rbindlist(sig_list)
  sig_dt$Plot_Region <- factor(sig_dt$Plot_Region, levels = regions_to_plot)

  ggplot(box_dt, aes(x = Group, y = Accessibility, fill = Plot_Region)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, linewidth = 0.2) +
    # 底层绘制：精准控制横线与竖线的位置，保证各区域绝对对齐且无重叠
    geom_segment(data = sig_dt, aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = sig_dt, aes(x = 1, xend = 1, y = y_bracket, yend = y_bracket - tip_len), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = sig_dt, aes(x = 2, xend = 2, y = y_bracket, yend = y_bracket - tip_len), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = sig_dt, aes(x = 1.5, y = y_text, label = p_sig), inherit.aes = FALSE, size = 2.2, fontface = "bold", vjust = 0) +
    facet_wrap(~ Plot_Region, nrow = 1) +
    scale_fill_manual(values = REGION_COLORS) +
    theme_bw(base_size = 7) +
    labs(title = NULL, x = "", y = "Acc.") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      plot.margin = MINI_MARGIN,
      strip.background = element_blank(),
      strip.text = element_text(size = 6, face = "bold")
    ) +
    coord_cartesian(ylim = c(0, y_max), clip = "off")
}

# ================= 2. 任务一：TTR 全景组合图 =================
log_msg("--- Task 1: TTR Combined Matrix (6 rows) ---")
ttr_regions <- c("Whole", "CDS", "3'UTR")
task1_plots <- list()

for (tool in names(PARAMS$TTR_based)) {
  wet_dt <- read_wet_data(tool, "TTR")
  for (param_type in c("TTR_based", "RNAxs_based")) {
    p <- PARAMS[[param_type]][[tool]]
    W <- p["W"]; L <- p["L"]; U <- p["U"]
    
    lunp_file <- run_rnaplfold(TTR_FASTA, W, L, TTR_LUNP)
    merged_dt <- build_3prime_matrix(lunp_file, wet_dt, U)
    
    mat_out <- file.path(TTR_MAT, sprintf("TTR_%s_%s_W%d_L%d_U%d.csv", tool, param_type, W, L, U))
    fwrite(merged_dt, mat_out)
    
    title_str <- sprintf("%s | %s (W%d L%d U%d)", tool, param_type, W, L, U)
    
    task1_plots <- c(task1_plots, list(
      get_track_plot(merged_dt, title_str),
      get_bar_plot(merged_dt, ttr_regions),
      get_box_plot(merged_dt, ttr_regions)
    ))
  }
}

merged_task1 <- ggarrange(
  plotlist = task1_plots, 
  ncol = 3, nrow = 6, 
  widths = c(4, 1, 2), 
  align = "h"
)
ggsave(file.path(OUT_DIR, "Task1_TTR_Combined_View.pdf"), plot = merged_task1, width = 7, height = 8)
log_msg("Task 1 PDF Generated: 7x8 inches layout.")

# ================= 3. 任务二：PCSK9 全景组合图 =================
log_msg("--- Task 2: PCSK9 Combined Matrix (2 rows) ---")
tool <- "shRNAone"
wet_dt <- read_wet_data(tool, "PCSK9")
pcsk9_regions <- c("Whole", "5'UTR", "CDS", "3'UTR")
task2_plots <- list()

for (param_type in c("TTR_based", "RNAxs_based")) {
  p <- PARAMS[[param_type]][[tool]]
  W <- p["W"]; L <- p["L"]; U <- p["U"]
  
  lunp_file <- run_rnaplfold(PCSK9_FASTA, W, L, PCSK9_LUNP)
  merged_dt <- build_3prime_matrix(lunp_file, wet_dt, U)
  
  mat_out <- file.path(PCSK9_MAT, sprintf("PCSK9_%s_%s_W%d_L%d_U%d.csv", tool, param_type, W, L, U))
  fwrite(merged_dt, mat_out)
  
  title_str <- sprintf("PCSK9 %s | %s (W%d L%d U%d)", tool, param_type, W, L, U)
  
  task2_plots <- c(task2_plots, list(
    get_track_plot(merged_dt, title_str),
    get_bar_plot(merged_dt, pcsk9_regions),
    get_box_plot(merged_dt, pcsk9_regions)
  ))
}

merged_task2 <- ggarrange(
  plotlist = task2_plots, 
  ncol = 3, nrow = 2, 
  widths = c(4, 1, 2), 
  align = "h"
)
ggsave(file.path(OUT_DIR, "Task2_PCSK9_Combined_View.pdf"), plot = merged_task2, width = 7, height = 3)
log_msg("Task 2 PDF Generated: 7x3 inches layout.")

# ================= 4. 任务三：siRNA 验证 =================
log_msg("--- Task 3: siRNA Matrix Evaluation ---")
cont_eval <- fread(file.path(SIRNA_EVAL_DIR, "01.Evaluation_Continuous.csv"))
bin_eval  <- fread(file.path(SIRNA_EVAL_DIR, "02.Evaluation_Binary.csv"))

siRNA_targets <- list(
  "shRNAone (TTR Best)"   = "W_35_L_20_U_7",
  "shRNAone (RNAxs Best)" = "W_80_L_40_U_7",
  "Baseline (U8)"         = "W_80_L_40_U_8",
  "Baseline (U16)"        = "W_80_L_40_U_16"
)

results <- list()
for (name in names(siRNA_targets)) {
  param_col <- siRNA_targets[[name]]
  r_val <- cont_eval[Parameter == param_col, Spearman_r]
  auc_val <- bin_eval[Parameter == param_col, AUC]
  if (length(r_val) > 0 && length(auc_val) > 0) {
    results[[name]] <- data.table(Group = name, Parameter = param_col, Spearman = r_val[1], AUC = auc_val[1])
  }
}

res_dt <- rbindlist(results)
fwrite(res_dt, file.path(OUT_DIR, "Task3_siRNA_Evaluation_Summary.csv"))

if (nrow(res_dt) > 0) {
  res_dt$Group <- factor(res_dt$Group, levels = names(siRNA_targets))
  melt_dt <- melt(res_dt, id.vars = c("Group", "Parameter"), measure.vars = c("Spearman", "AUC"), variable.name = "Metric", value.name = "Performance")
  
  p_sirna <- ggplot(melt_dt, aes(x = Group, y = Performance)) +
    geom_bar(stat = "identity", position = "dodge", fill = "gray70", color = "black", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Performance)), vjust = -0.6, size = 2.2, fontface = "bold") +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
    theme_bw(base_size = 8) +
    labs(title = "siRNA Independent Verification Comparison", x = "", y = "Score") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 6),
      axis.text.y = element_text(size = 6),
      legend.position = "none",
      strip.background = element_rect(fill = "#F0F3F4"),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(size = 9, face = "bold")
    )
  ggsave(file.path(OUT_DIR, "Task3_siRNA_Performance_Bar.pdf"), p_sirna, width = 4, height = 2)
}

log_msg("--- All Tasks Completed Successfully ---")