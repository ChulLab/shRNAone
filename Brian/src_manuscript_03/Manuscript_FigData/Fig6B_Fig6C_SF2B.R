#!/usr/bin/env Rscript
# derived from: partD.02.best_combination_display.R
#
# Fig6B: siRNA Matrix Evaluation
# Fig6C: PCSK9 Combined Matrix (TTR_based)
# SF2B: PCSK9 Combined Matrix (RNAxs_based)


library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. Core Configuration =================
log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

PARAMS <- list(
  TTR_based = list(shRNAone = c(W=35, L=20, U=7)),
  RNAxs_based = list(shRNAone = c(W=80, L=40, U=7))
)

# Path configuration
PCSK9_MAT <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.PCSK9_pred/03.matrix"
SIRNA_EVAL_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partC/03.GT_siRNA_acc/02.siRNA_matrix"
OUT_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/Manuscript_FigData/main"

invisible(lapply(c(OUT_DIR), dir.create, showWarnings = FALSE, recursive = TRUE))

REGION_COLORS <- c("Whole" = "gray70", "5'UTR" = "#3498DB", "CDS" = "#21908C", "3'UTR" = "#FDE725")
MINI_MARGIN <- margin(2, 2, 2, 2, "pt")

# ================= 1. Core Plotting Functions (Task 2) =================

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
    geom_text(aes(label = label_text, y = label_y), size = 8 / .pt, family = "Arial", vjust = ifelse(df_bar$r >= 0, 0, 1)) + 
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) + 
    scale_fill_manual(values = REGION_COLORS) + 
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
    labs(x = NULL, y = "Spearman r") +
    theme_minimal(base_size = 8, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
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
  
  # Convert Accessibility to log10(Acc + 0.01)
  box_dt$Accessibility <- log10(box_dt$Accessibility + 0.01)
  
  upper_lim <- max(sapply(split(box_dt, list(box_dt$Plot_Region, box_dt$Group)), function(sub) {
    if(nrow(sub) < 5) return(-Inf)
    boxplot.stats(sub$Accessibility)$stats[5]
  }), na.rm = TRUE)
  if(is.na(upper_lim) || is.infinite(upper_lim)) upper_lim <- max(box_dt$Accessibility, na.rm = TRUE)
  
  # Calculate dynamic margins based on log scale range because limits can be negative
  y_min_val <- min(box_dt$Accessibility, na.rm = TRUE)
  y_span <- upper_lim - y_min_val
  if (y_span == 0) y_span <- 1
  
  y_max <- upper_lim + y_span * 0.35 
  y_bracket <- upper_lim + y_span * 0.15
  y_text <- upper_lim + y_span * 0.18
  tip_len <- y_span * 0.03
  
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
    geom_segment(data = sig_dt, aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = sig_dt, aes(x = 1, xend = 1, y = y_bracket, yend = y_bracket - tip_len), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = sig_dt, aes(x = 2, xend = 2, y = y_bracket, yend = y_bracket - tip_len), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = sig_dt, aes(x = 1.5, y = y_text, label = p_sig), inherit.aes = FALSE, size = 8 / .pt, family = "Arial", vjust = 0) +
    facet_wrap(~ Plot_Region, nrow = 1) +
    scale_fill_manual(values = REGION_COLORS) +
    theme_bw(base_size = 8, base_family = "Arial") +
    labs(title = NULL, x = "", y = expression(log[10](Acc+0.01))) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = MINI_MARGIN,
      strip.background = element_blank()
    ) +
    coord_cartesian(ylim = c(y_min_val, y_max), clip = "off")
}


# ================= 2. Task 2: PCSK9 Combined Matrix (No Track) =================
log_msg("--- Task 2: PCSK9 Combined Matrix (No Track) ---")
tool <- "shRNAone"
pcsk9_regions <- c("Whole", "5'UTR", "CDS", "3'UTR")
task2_plots <- list()

for (param_type in names(PARAMS)) {
  p <- PARAMS[[param_type]][[tool]]
  W <- p["W"]; L <- p["L"]; U <- p["U"]
  
  mat_in <- file.path(PCSK9_MAT, sprintf("PCSK9_%s_%s_W%d_L%d_U%d.csv", tool, param_type, W, L, U))
  if(!file.exists(mat_in)) {
      log_msg(sprintf("Warning: Cannot find %s", mat_in))
      next
  }
  
  merged_dt <- fread(mat_in)
  
  task2_plots <- c(task2_plots, list(
    get_bar_plot(merged_dt, pcsk9_regions),
    get_box_plot(merged_dt, pcsk9_regions)
  ))
}

if(length(task2_plots) > 0) {
    merged_task2 <- ggarrange(
      plotlist = task2_plots, 
      ncol = 2, nrow = 2, 
      widths = c(1, 2.5), 
      align = "h"
    )
    
    out_task2 <- file.path(OUT_DIR, "Task2_PCSK9_Combined_View.pdf")
    # Height modified to 2.4 as required
    ggsave(out_task2, plot = merged_task2, width = 4, height = 3, device = cairo_pdf)
    log_msg("Task 2 PDF Generated: Width 3.75, Height 2.4 layout.")
}

# --- Save PCSK9 correlation + box stats ---
pcs_corr <- list()
pcs_box  <- list()
for (param_type in names(PARAMS)) {
  p <- PARAMS[[param_type]][[tool]]
  W <- p["W"]; L <- p["L"]; U <- p["U"]
  mat_in <- file.path(PCSK9_MAT, sprintf("PCSK9_%s_%s_W%d_L%d_U%d.csv", tool, param_type, W, L, U))
  if (!file.exists(mat_in)) next
  dt <- fread(mat_in)
  for (reg in pcsk9_regions) {
    sub <- if (reg == "Whole") dt else dt[region == reg]
    if (nrow(sub) >= 3) {
      ct <- cor.test(sub$Accessibility, sub$log2FC, method = "spearman", exact = FALSE)
      pcs_corr[[length(pcs_corr) + 1]] <- data.table(
        Tool = tool, Param_Type = param_type, Region = reg,
        W = W, L = L, U = U, Spearman_r = ct$estimate, P_value = ct$p.value, N = nrow(sub)
      )
    }
  }
  q20 <- quantile(dt$log2FC, 0.2, na.rm = TRUE)
  q80 <- quantile(dt$log2FC, 0.8, na.rm = TRUE)
  for (reg in pcsk9_regions) {
    sub <- if (reg == "Whole") dt else dt[region == reg]
    top <- sub[log2FC >= q80, Accessibility]
    bot <- sub[log2FC <= q20, Accessibility]
    if (length(top) >= 3 && length(bot) >= 3) {
      wt <- wilcox.test(top, bot)
      pcs_box[[length(pcs_box) + 1]] <- data.table(
        Tool = tool, Param_Type = param_type, Region = reg,
        W = W, L = L, U = U,
        Median_Acc_Top20 = median(top, na.rm = TRUE),
        Median_Acc_Bot20 = median(bot, na.rm = TRUE),
        N_Top = length(top), N_Bot = length(bot),
        Wilcoxon_P = wt$p.value
      )
    }
  }
}
fwrite(rbindlist(pcs_corr), file.path(OUT_DIR, "Task2_PCSK9_Correlation_Stats.csv"))
fwrite(rbindlist(pcs_box),  file.path(OUT_DIR, "Task2_PCSK9_Box_Wilcoxon_Stats.csv"))
log_msg("Saved Task2 CSV statistics.")


# ================= 3. Task 3: siRNA Matrix Evaluation (n labeled in title, no main title) =================
log_msg("--- Task 3: siRNA Matrix Evaluation ---")

cont_file <- file.path(SIRNA_EVAL_DIR, "01.Evaluation_Continuous.csv")
bin_file  <- file.path(SIRNA_EVAL_DIR, "02.Evaluation_Binary.csv")

if(file.exists(cont_file) && file.exists(bin_file)) {
    cont_eval <- fread(cont_file)
    bin_eval  <- fread(bin_file)
    
    siRNA_targets <- list(
      "W35 L20 U7"   = "W_35_L_20_U_7",
      "W80 L40 U7"   = "W_80_L_40_U_7",
      "W80 L40 U8"   = "W_80_L_40_U_8",
      "W80 L40 U16"  = "W_80_L_40_U_16"
    )
    
    results <- list()
    for (name in names(siRNA_targets)) {
      param_col <- siRNA_targets[[name]]
      
      r_val <- cont_eval[Parameter == param_col, Spearman_r]
      n_cont <- cont_eval[Parameter == param_col, Valid_N]
      
      auc_val <- bin_eval[Parameter == param_col, AUC]
      n_func <- bin_eval[Parameter == param_col, N_Func]
      n_nonfunc <- bin_eval[Parameter == param_col, N_NonFunc]
      
      if (length(r_val) > 0 && length(auc_val) > 0) {
        results[[name]] <- data.table(
          Group = name, 
          Spearman = r_val[1], 
          Spearman_N = n_cont[1],
          AUC = auc_val[1],
          AUC_N = (n_func[1] + n_nonfunc[1])
        )
      }
    }
    
    res_dt <- rbindlist(results)
    
    if (nrow(res_dt) > 0) {
      # Dynamically retrieve n and label in the title
      spearman_n_val <- max(res_dt$Spearman_N, na.rm = TRUE)
      auc_n_val <- max(res_dt$AUC_N, na.rm = TRUE)
      
      spearman_label <- sprintf("Spearman (n=%d)", spearman_n_val)
      auc_label <- sprintf("AUC (n=%d)", auc_n_val)
      
      melt_dt <- data.table(
         Group = rep(names(siRNA_targets), 2),
         Metric = factor(rep(c(spearman_label, auc_label), each = length(siRNA_targets)), levels = c(spearman_label, auc_label)),
         Performance = c(res_dt$Spearman, res_dt$AUC)
      )
      
      melt_dt$Group <- factor(melt_dt$Group, levels = names(siRNA_targets))
      melt_dt$Group_num <- as.numeric(melt_dt$Group)
      
      melt_dt$Y_min <- ifelse(grepl("Spearman", melt_dt$Metric), 0.25, 0.70)
      
      p_sirna <- ggplot(melt_dt) +
        geom_hline(aes(yintercept = Y_min), linewidth = 0.3, color = "black") +
        geom_rect(aes(xmin = Group_num - 0.35, xmax = Group_num + 0.35, 
                      ymin = Y_min, ymax = Performance), 
                  fill = "gray70", color = "black", linewidth = 0.3) +
        # Remove n above the bar, keeping only the value
        geom_text(aes(x = Group_num, y = Performance, 
                      label = sprintf("%.2f", Performance)), 
                  vjust = -0.4, size = 8 / .pt, family = "Arial") +
        facet_wrap(~ Metric, scales = "free_y") +
        scale_x_continuous(breaks = 1:length(names(siRNA_targets)), labels = names(siRNA_targets)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + 
        theme_bw(base_size = 8, base_family = "Arial") +
        # Remove the main title
        labs(title = NULL, x = "", y = "Score") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          strip.background = element_blank()
        )
      
      out_task3 <- file.path(OUT_DIR, "Task3_siRNA_Performance_Bar.pdf")
      # Height compressed to 1.25
      ggsave(out_task3, plot = p_sirna, width = 4, height = 2, device = cairo_pdf)
      log_msg("Task 3 PDF Generated: Width 4, Height 2 layout.")

      fwrite(res_dt, file.path(OUT_DIR, "Task3_siRNA_Evaluation_Summary.csv"))
      log_msg("Saved Task3 CSV statistics.")
    }
} else {
    log_msg("Warning: siRNA evaluation files not found.")
}

log_msg("--- All Plotting Tasks Completed Successfully ---")