# 16.offset_dry_innate.R — Pairwise offset scanning between ALL dry prediction models (Heatmap).
#
# Scans offsets from -50 to +50 nt for all pairs of dry models.
# Positive offset means the Shift model's position is downstream of the Base model.
#
# Generates:
#   1. 16.offset_dry_innate.csv (correlation statistics for the BEST offsets of all pairs)
#   2. 16.offset_dry_innate.pdf (Heatmap with row=Base, col=Shift, annotating offset, r, and p)
#
# conda env: system R

library(ggplot2)
library(grid)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush=TRUE)
}

# ================= Configuration =================
PROJECT <- "/data/cai801/data/HKUcas"
MODEL_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
FIGURE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")

dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# Analytical settings
COR_METHOD <- "pearson"  # Options: "spearman" or "pearson"
ALIGN_BY <- "start_pos"     # Options: "end_pos" or "start_pos"
OFFSET_RANGE <- 50

# All dry models for the All-vs-All matrix
ALL_MODELS <- c(
  "Sanjana2020_23", "TIGER2023_23", "Hsu2023_30",
  "DeepCas13_18", "DeepCas13_23", "DeepCas13_30",
  "Fareh2024_30", "Sanjana2020_30"
)

format_p_short <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-10) return("<1e-10")
  if (p < 1e-5) return("<1e-5")
  if (p < 1e-3) return("<1e-3")
  sprintf("%.3f", p)
}

method_label <- tools::toTitleCase(COR_METHOD)
log_msg(sprintf("02. Offset scanning: All-vs-All dry models (align by: %s, method: %s)", ALIGN_BY, method_label))

# Pre-load all models to save I/O time
loaded_models <- list()
for (m in ALL_MODELS) {
  file_path <- file.path(MODEL_DIR, paste0(m, "_aligned.csv"))
  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    if (ALIGN_BY %in% colnames(df)) {
      df <- df[, c(ALIGN_BY, "prediction_score")]
      colnames(df) <- c("pos", "score")
      df <- aggregate(score ~ pos, data = df, FUN = mean)
      loaded_models[[m]] <- df
    }
  }
}

best_stats <- data.frame()

# ================= Matrix Scanning =================
for (base_name in names(loaded_models)) {
  df_base <- loaded_models[[base_name]]
  colnames(df_base) <- c("pos_base", "score_base")
  
  for (shift_name in names(loaded_models)) {
    log_msg(sprintf("  Scanning: [Base] %s vs [Shift] %s", base_name, shift_name))
    
    # Self-comparison optimization
    if (base_name == shift_name) {
      best_stats <- rbind(best_stats, data.frame(
        Base = base_name, Shift = shift_name,
        best_offset = 0, cor_r = 1.0, p_value = 0,
        p_fmt = "<1e-10"
      ))
      next
    }
    
    df_shift <- loaded_models[[shift_name]]
    colnames(df_shift) <- c("pos_shift", "score_shift")
    
    best_r <- -Inf
    best_offset <- NA
    best_p <- NA
    
    for (offset in -OFFSET_RANGE:OFFSET_RANGE) {
      # Alignment rule: pos_base + offset = pos_shift 
      df_shift$target_pos_base <- df_shift$pos_shift - offset
      merged <- merge(df_base, df_shift, by.x = "pos_base", by.y = "target_pos_base")
      merged <- na.omit(merged)
      
      if (nrow(merged) >= 3) {
        test <- cor.test(merged$score_base, merged$score_shift, method = COR_METHOD, exact = FALSE)
        r <- test$estimate
        p <- test$p.value
        
        # Keep highest real correlation (positive peak)
        if (!is.na(r) && r > best_r) {
          best_r <- r
          best_offset <- offset
          best_p <- p
        }
      }
    }
    
    if (!is.na(best_offset)) {
      best_stats <- rbind(best_stats, data.frame(
        Base = base_name, Shift = shift_name,
        best_offset = best_offset, cor_r = best_r, p_value = best_p,
        p_fmt = format_p_short(best_p)
      ))
    }
  }
}

# ================= 1. Save Data =================
out_csv <- file.path(FIGURE_DIR, "02.offset_dry_innate.csv")
write.csv(best_stats, out_csv, row.names = FALSE)
log_msg(sprintf("  Saved heatmap stats data: %s", basename(out_csv)))

# ================= 2. Plotting (Heatmap) =================
# Format Base and Shift as factors to lock order
best_stats$Base <- factor(best_stats$Base, levels = rev(ALL_MODELS)) # Reverse so first model is at the top
best_stats$Shift <- factor(best_stats$Shift, levels = ALL_MODELS)

# Prepare cell label text
best_stats$cell_label <- sprintf("Off: %d\nr: %.2f\np: %s", 
                                 best_stats$best_offset, 
                                 best_stats$cor_r, 
                                 best_stats$p_fmt)

p_heat <- ggplot(best_stats, aes(x = Shift, y = Base, fill = cor_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = cell_label), color = ifelse(best_stats$cor_r > 0.6, "white", "black"), size = 2.5, lineheight = 0.9) +
  scale_fill_viridis_c(name = sprintf("%s r", method_label), option = "mako", direction = -1, limits = c(min(best_stats$cor_r), 1)) +
  labs(
    title = "Innate Dry-Lab Comparison: Pairwise Best Offset Heatmap",
    subtitle = sprintf("Correlation: %s | Alignment anchor: %s | Scanning Range: ±%d nt", method_label, ALIGN_BY, OFFSET_RANGE),
    x = "Shift Model",
    y = "Base Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

pdf_path <- file.path(FIGURE_DIR, "02.offset_dry_innate.pdf")
pdf(pdf_path, width = 11, height = 9) 
print(p_heat)
dev.off()
log_msg(sprintf("  Saved heatmap plot: %s", basename(pdf_path)))

log_msg("16. Done.")