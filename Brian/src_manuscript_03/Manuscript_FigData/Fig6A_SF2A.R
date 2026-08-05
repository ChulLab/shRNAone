#!/usr/bin/env Rscript
# derived from: partD.02.best_combination_display.R
#
# Fig6A: TTR-based comparison for shRNAone, CasRx, and PspCas13b (full)
# SF2A: RNAxs-based comparison for shRNAone, CasRx, and PspCas13b (without track plot)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ================= 0. Core Configuration =================
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

# Input path for pre-calculated matrices
TTR_MAT <- "/data/cai801/data/HKUcas/result/08.manuscript_03/partD/00.TTR_pred/03.matrix"
# Output path for new manuscript figures
OUT_DIR <- "/data/cai801/data/HKUcas/result/08.manuscript_03/Manuscript_FigData/main"

invisible(lapply(c(OUT_DIR), dir.create, showWarnings = FALSE, recursive = TRUE))

REGION_COLORS <- c("Whole" = "gray70", "5'UTR" = "#3498DB", "CDS" = "#21908C", "3'UTR" = "#FDE725")
MINI_MARGIN <- margin(2, 2, 2, 2, "pt")

# ================= 1. Core Plotting Functions =================

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
    theme_minimal(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(margin = margin(b = 2)),
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
  
  # Convert Accessibility to log10(Acc + 0.01) for the rightmost plot display
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

# ================= 2. Task 1: Read direct data and plot separately =================
log_msg("--- Task 1: TTR Separate Matrices (TTR_based vs RNAxs_based) ---")
ttr_regions <- c("Whole", "CDS", "3'UTR")

# Unified order: shRNAone on top, CasRx in the middle, and PspCas13b at the bottom
tools <- c("shRNAone", "CasRx_day5", "PspCas13b")

corr_stats <- list()
box_stats  <- list()

for (param_type in names(PARAMS)) {
  log_msg(sprintf("Processing group: %s", param_type))
  current_plots <- list()
  
  for (tool in tools) {
    p <- PARAMS[[param_type]][[tool]]
    W <- p["W"]; L <- p["L"]; U <- p["U"]
    
    # Construct input matrix file path, directly read previously calculated results
    mat_in <- file.path(TTR_MAT, sprintf("TTR_%s_%s_W%d_L%d_U%d.csv", tool, param_type, W, L, U))
    
    if (!file.exists(mat_in)) {
      log_msg(sprintf("Warning: Data matrix not found for %s (%s). Expected path: %s", tool, param_type, mat_in))
      next
    }
    
    merged_dt <- fread(mat_in)
    
    # --- Collect correlation stats ---
    for (reg in ttr_regions) {
      sub <- if (reg == "Whole") merged_dt else merged_dt[region == reg]
      if (nrow(sub) >= 3) {
        ct <- cor.test(sub$Accessibility, sub$log2FC, method = "spearman", exact = FALSE)
        corr_stats[[length(corr_stats) + 1]] <- data.table(
          Tool = tool, Param_Type = param_type, Region = reg,
          W = W, L = L, U = U,
          Spearman_r = ct$estimate, P_value = ct$p.value, N = nrow(sub)
        )
      }
    }
    # --- Collect box (Top vs Bot 20%) Wilcoxon stats ---
    q20 <- quantile(merged_dt$log2FC, 0.2, na.rm = TRUE)
    q80 <- quantile(merged_dt$log2FC, 0.8, na.rm = TRUE)
    for (reg in ttr_regions) {
      sub <- if (reg == "Whole") merged_dt else merged_dt[region == reg]
      top <- sub[log2FC >= q80, Accessibility]
      bot <- sub[log2FC <= q20, Accessibility]
      med_top <- median(log10(top + 0.01), na.rm = TRUE)
      med_bot <- median(log10(bot + 0.01), na.rm = TRUE)
      if (length(top) >= 3 && length(bot) >= 3) {
        wt <- wilcox.test(top, bot)
        box_stats[[length(box_stats) + 1]] <- data.table(
          Tool = tool, Param_Type = param_type, Region = reg,
          W = W, L = L, U = U,
          Median_Acc_Top20_log10 = med_top, Median_Acc_Bot20_log10 = med_bot,
          N_Top = length(top), N_Bot = length(bot),
          Wilcoxon_P = wt$p.value
        )
      }
    }
    
    # Formatting display names to accurately match the image's layout
    display_tool <- ifelse(tool == "CasRx_day5", "CasRx", tool)
    display_param <- sub("_based", "", param_type)
    title_str <- sprintf("%s | %s (W%d L%d U%d)", display_tool, display_param, W, L, U)
    
    # Store the three plots of the tool into the current group's list
    current_plots <- c(current_plots, list(
      get_track_plot(merged_dt, title_str),
      get_bar_plot(merged_dt, ttr_regions),
      get_box_plot(merged_dt, ttr_regions)
    ))
  }
  
  # If data is successfully read, arrange and save
  if (length(current_plots) > 0) {
    # Each PDF contains only 3 tools under the current parameter group (3 rows in total)
    merged_fig <- ggarrange(
      plotlist = current_plots, 
      ncol = 3, nrow = length(tools), 
      widths = c(4, 1, 2), 
      align = "h"
    )
    
    # Output file is an independent PDF, width and height set to 7.5 x 3.5
    out_pdf <- file.path(OUT_DIR, sprintf("Task1_TTR_%s_View.pdf", param_type))
    # ! Added device = cairo_pdf to ggsave here to solve system font mapping error
    ggsave(out_pdf, plot = merged_fig, width = 7.5, height = 3.5, device = cairo_pdf)
    log_msg(sprintf("Saved PDF: %s (Width: 7.5, Height: 3.5)", out_pdf))
  }
}

# --- Save collected stats ---
fwrite(rbindlist(corr_stats), file.path(OUT_DIR, "Task1_TTR_Correlation_Stats.csv"))
fwrite(rbindlist(box_stats),  file.path(OUT_DIR, "Task1_TTR_Box_Wilcoxon_Stats.csv"))
log_msg("Saved Task1 CSV statistics.")

log_msg("--- All Plotting Tasks Completed Successfully ---")