#!/usr/bin/env Rscript
# derived from: partB.04.parameter_impact_Plotter.R
#
# SF1C: Random Forest (Non-linear Importance) for W, L, and U parameters

library(ggplot2)
library(data.table)

PROJECT <- "/data/cai801/data/HKUcas"
CACHE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/04.parameter_impact/cache")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Reverting to high-grade grayscale to avoid cognitive confusion with Region/Direction colors
PARAM_COLORS <- c("W" = "gray40", "L" = "gray65", "U" = "gray85")
BASE_SIZE <- 7 

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. Data Loading & Preprocessing =================
rf_res <- fread(file.path(CACHE_DIR, "03_RF_res.csv"))

# Calculate Relative Importance
rf_res[, Rel_Imp := IncMSE / sum(IncMSE) * 100, by = Tool]

# Enforce tool order from left to right
rf_res[, Tool := factor(Tool, levels = c("shRNAone", "CasRx_day5", "PspCas13b"))]

# Enforce parameter order from left to right
rf_res[, Parameter := factor(Parameter, levels = c("W", "L", "U"))]


# ================= 2. Plotting Design =================
log_msg("Rendering Random Forest (Non-linear Importance) Plot...")

p_rf <- ggplot(rf_res, aes(x = Parameter, y = Rel_Imp, fill = Parameter)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.7) +
  facet_wrap(~ Tool, nrow = 1) +
  scale_fill_manual(values = PARAM_COLORS) +
  labs(title = "Non-linear Global Importance (RF)", y = "Rel. Imp. (%)", x = "Thermodynamic Parameter") +
  theme_bw(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = BASE_SIZE + 1, margin = margin(t = 2, b = 6), hjust = 0.5),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", face = "bold", margin = margin(t = 2, b = 2)),
    axis.text.x = element_text(color = "black"),      
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.margin = margin(2, 2, 2, 2) 
  )

# ================= 3. Output =================
out_pdf <- file.path(OUT_DIR, "Supple_Parameter_Impact_RF.pdf")

# Adjusted height to 1.5 since it is now a single-panel plot
ggsave(out_pdf, plot = p_rf, width = 2.7, height = 1.2, device = cairo_pdf)

# --- Save RF importance ---
fwrite(rf_res, file.path(OUT_DIR, "SF1C_RF_Importance.csv"))
log_msg(sprintf("Successfully generated compact dense plot: %s", out_pdf))