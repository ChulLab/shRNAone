# 15.offset_wet_innate.R — Pairwise offset scanning between wet-lab methods.
#
# Compares 3 wet tools pairwise directly in memory:
#   - shRNAone (Base) vs CasRx_day5 (Shift)
#   - shRNAone (Base) vs PspCas13b (Shift)
#   - CasRx_day5 (Base) vs PspCas13b (Shift)
#
# Logic: Aligns by ALIGN_BY (e.g., end_pos or start_pos), scans offsets from -50 to +50 nt.
# Positive offset means the Shift tool's position is downstream of the Base tool.
#
# Generates:
#   1. 01.offset_wet_innate.csv (correlation statistics for all offsets)
#   2. 01.offset_wet_innate.pdf (line chart with best offset highlighted)
#
# conda env: system R

library(ggplot2)
library(grid)
library(ggrepel) # Added for non-overlapping text annotations

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush=TRUE)
}

# ================= Configuration =================
PROJECT <- "/data/cai801/data/HKUcas"
WET_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
FIGURE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# Analytical settings
COR_METHOD <- "pearson"  # Options: "spearman" or "pearson"
ALIGN_BY <- "start_pos"     # Options: "end_pos" or "start_pos" (ensure this column exists in wet-lab csvs)
OFFSET_RANGE <- 50

# File names for the final wet datasets
WET_TOOLS <- list(
  "shRNAone" = "TTR_shRNAone_clean.csv",
  "PspCas13b" = "TTR_PspCas13b_clean.csv",
  "CasRx_day5" = "TTR_CasRx_day5_clean.csv"
)

# Pair combinations to test (Base vs Shift)
PAIRS <- list(
  c("shRNAone", "CasRx_day5"),
  c("shRNAone", "PspCas13b"),
  c("CasRx_day5", "PspCas13b")
)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-10) return("< 1e-10")
  if (p < 1e-5) return("< 1e-5")
  if (p < 1e-3) return("< 1e-3")
  sprintf("%.4f", p)
}

# Dynamic labels based on selected method
method_label <- tools::toTitleCase(COR_METHOD)
log_msg(sprintf("01. Offset scanning: wet vs wet (innate pairwise, align by: %s, method: %s)", ALIGN_BY, method_label))

all_offset_results <- data.frame()
best_stats <- data.frame()

for (pair in PAIRS) {
  tool_a <- pair[1] # Base
  tool_b <- pair[2] # Shift
  # Clarify Base and Shift in the legend
  pair_name <- sprintf("%s (Base) vs %s (Shift)", tool_a, tool_b)
  log_msg(sprintf("\n  === %s ===", pair_name))
  
  # Load wet data
  file_a <- file.path(WET_DIR, WET_TOOLS[[tool_a]])
  file_b <- file.path(WET_DIR, WET_TOOLS[[tool_b]])
  
  if (!file.exists(file_a) || !file.exists(file_b)) {
    log_msg(sprintf("    [SKIP] Missing file for %s", pair_name))
    next
  }
  
  df_a <- read.csv(file_a)
  df_b <- read.csv(file_b)
  
  # Validate column existence
  if (!(ALIGN_BY %in% colnames(df_a)) || !(ALIGN_BY %in% colnames(df_b))) {
    log_msg(sprintf("    [ERROR] Column '%s' not found in one of the dataframes. Skipping pair.", ALIGN_BY))
    next
  }
  
  # Keep only the alignment column and mean_log2FC to avoid clutter
  df_a <- df_a[, c(ALIGN_BY, "mean_log2FC")]
  df_b <- df_b[, c(ALIGN_BY, "mean_log2FC")]
  colnames(df_a) <- c("pos_a", "log2FC_a")
  colnames(df_b) <- c("pos_b", "log2FC_b")
  
  pair_results <- data.frame()
  best_r <- -Inf
  best_offset <- NA
  best_p <- NA
  
  for (offset in -OFFSET_RANGE:OFFSET_RANGE) {
    # Alignment rule: pos_a + offset = pos_b 
    # Therefore we shift B back to A's coordinate space for joining
    df_b$target_pos_a <- df_b$pos_b - offset
    
    # Inner merge to only keep overlapping targets
    merged <- merge(df_a, df_b, by.x = "pos_a", by.y = "target_pos_a")
    merged <- na.omit(merged)
    
    if (nrow(merged) >= 3) {
      test <- cor.test(merged$log2FC_a, merged$log2FC_b, method = COR_METHOD, exact = FALSE)
      r <- test$estimate
      p <- test$p.value
      
      pair_results <- rbind(pair_results, data.frame(
        pair = pair_name,
        offset = offset,
        correlation_r = r,
        correlation_p = p,
        n_points = nrow(merged)
      ))
      
      # STRICT RULE: Must be the highest real correlation (positive peak)
      if (!is.na(r) && r > best_r) {
        best_r <- r
        best_offset <- offset
        best_p <- p
      }
    }
  }
  
  all_offset_results <- rbind(all_offset_results, pair_results)
  
  best_stats <- rbind(best_stats, data.frame(
    pair = pair_name,
    best_offset = best_offset,
    cor_r = best_r,
    p_value = best_p,
    p_fmt = format_p(best_p)
  ))
  
  log_msg(sprintf("    Best offset: %d nt, r = %.4f", best_offset, best_r))
}

# ================= 1. Save Data =================
out_csv <- file.path(FIGURE_DIR, "01.offset_wet_innate.csv")
write.csv(all_offset_results, out_csv, row.names = FALSE)
log_msg(sprintf("  Saved offset scan data: %s", basename(out_csv)))

# ================= 2. Plotting =================
# Prepare label text for best points dynamically
best_stats$label_text <- sprintf("%d nt\nr=%.3f, p=%s", 
                                 best_stats$best_offset, 
                                 best_stats$cor_r, 
                                 best_stats$p_fmt)

p_offset <- ggplot(all_offset_results, aes(x = offset, y = correlation_r, color = pair, group = pair)) +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed", alpha = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2, alpha = 0.6) +
  # Highlight BEST offset points
  geom_point(data = best_stats, 
             aes(x = best_offset, y = cor_r, color = pair), 
             inherit.aes = FALSE, size = 3.5, shape = 21, fill = "white", stroke = 1.2, 
             show.legend = FALSE) +
  # Repel overlapping text annotations
  geom_text_repel(data = best_stats, 
                  aes(x = best_offset, y = cor_r, label = label_text), 
                  inherit.aes = FALSE, color = "black", size = 3, 
                  box.padding = 0.6, point.padding = 0.3, 
                  min.segment.length = 0, segment.color = "grey50",
                  show.legend = FALSE) +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(clip = "off") + 
  labs(
    title = sprintf("Innate Wet-Lab Comparison: Offset Scan (%s)", method_label),
    x = sprintf("Offset of Shift (nt) relative to Base (Aligned by %s)", ALIGN_BY),
    y = sprintf("%s r", method_label),
    color = "Pairwise Comparison"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "vertical", # Stack legend items vertically if too long
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10)
  )

pdf_path <- file.path(FIGURE_DIR, "01.offset_wet_innate.pdf")
pdf(pdf_path, width = 8, height = 5.5) # Slightly taller to fit stacked legends
print(p_offset)
dev.off()
log_msg(sprintf("  Saved plot: %s", basename(pdf_path)))

log_msg("15. Done.")