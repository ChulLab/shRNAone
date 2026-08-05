#!/usr/bin/env Rscript
# derived from: partC.01.define_refU_plateau.R
# 
# SF1D: define plateau region for across all tools based on W=80, L=40, whole region, 3prime direction

library(ggplot2)
library(dplyr)
library(data.table)

# ================= 0. Configuration & Initialization =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/supple")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Enforce left-to-right presentation order
TOOLS <- c("shRNAone", "CasRx_day5", "PspCas13b")
TARGET_DIR <- "3prime"   # Inherit Part B conclusion, only analyze 3prime
TARGET_REGION <- "whole" # Based on requirements, only analyze whole region

PLATEAU_THRESHOLD <- 0.75
COLOR_3P <- "#D28E72"    # Use the established Morandi color for 3prime
BASE_SIZE <- 7

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. Rapid Loading & Filtering =================
log_msg("Loading and filtering data...")
cols_to_keep <- c("W", "L", "U", "direction", "region", "spearman_r", "spearman_p")

all_data_list <- lapply(TOOLS, function(t) {
  # Wait to use original filenames even if plotting order is different
  fpath <- file.path(IN_DIR, paste0(t, "_grid.csv")) 
  if (file.exists(fpath)) {
    dt <- fread(fpath, select = cols_to_keep)
    dt$tool <- t
    return(dt)
  }
})
ref_data <- rbindlist(all_data_list)

# Core filter: W=80, L=40, whole region, 3prime direction
# Define significance label for conditional plotting
sub_data <- ref_data %>%
  filter(W == 80, L == 40, region == TARGET_REGION, direction == TARGET_DIR) %>%
  mutate(tool = factor(tool, levels = TOOLS)) %>%
  mutate(Sig_Label = ifelse(spearman_p < 0.05 & spearman_r > 0, "p < 0.05", "ns"))

if(nrow(sub_data) == 0) stop("No data found after filtering. Please check input files.")

# ================= 2. Plateau Calculation =================
log_msg("Calculating Peak and Plateau statistics...")
plateau_stats <- sub_data %>%
  group_by(tool) %>%
  summarise(
    Max_Cor = max(spearman_r, na.rm = TRUE),
    Peak_U  = U[which.max(spearman_r)],
    Plateau_Min = min(U[spearman_r >= PLATEAU_THRESHOLD * max(spearman_r, na.rm = TRUE) & spearman_r > 0], na.rm = TRUE),
    Plateau_Max = max(U[spearman_r >= PLATEAU_THRESHOLD * max(spearman_r, na.rm = TRUE) & spearman_r > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Break label into two lines to fit the ultra-compact 2.7-inch width
    Label_Text = sprintf("Peak U=%d\n(Plat: %d-%d)", Peak_U, Plateau_Min, Plateau_Max),
    Label_Y = Max_Cor + 0.15 # Added slight padding for two-line text
  )

# ================= 3. Compact High-Density Plotting =================
log_msg("Rendering compact plot (2.7 x 1.8 inches)...")

y_max_limit <- max(plateau_stats$Label_Y, na.rm = TRUE) + 0.05
y_min_limit <- min(sub_data$spearman_r, na.rm = TRUE) - 0.05

p_traj <- ggplot(sub_data, aes(x = U, y = spearman_r)) +
  # Plateau shaded background
  geom_rect(data = plateau_stats %>% filter(!is.na(Plateau_Min)), 
            aes(xmin = Plateau_Min, xmax = Plateau_Max, ymin = -Inf, ymax = Inf), 
            fill = COLOR_3P, alpha = 0.2, inherit.aes = FALSE, color = NA) +
  
  # Y=0 reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  
  # Trajectory line
  geom_line(color = COLOR_3P, alpha = 0.9, linewidth = 0.8) +
  
  # Conditional significance points mapped to shape 
  geom_point(aes(shape = Sig_Label), color = COLOR_3P, size = 1.0) +
  scale_shape_manual(values = c("p < 0.05" = 16, "ns" = NA), name = "Significance") +
  
  # Dual-line label (centered in each panel's upper bound)
  geom_text(data = plateau_stats,
            aes(x = mean(range(sub_data$U)), y = Label_Y, label = Label_Text),
            color = "black", size = 5 / .pt, family = "Arial", fontface = "bold", lineheight = 0.9) +
  
  facet_wrap(~ tool, nrow = 1) +
  scale_y_continuous(limits = c(y_min_limit, y_max_limit)) +
  labs(x = "Unpaired Length Constraint (U)", y = "Spearman (r)") +
  theme_bw(base_size = BASE_SIZE, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(t = -8, b = 0, l = 0, r = 0),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", color = "black", margin = margin(t=2, b=2)),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.margin = margin(2, 2, 2, 2)
  )

out_pdf <- file.path(OUT_DIR, "Supple_Ref_U_Plateau_Compact.pdf")
ggsave(out_pdf, plot = p_traj, width = 2.7, height = 1.3, device = cairo_pdf)

# --- Save plateau stats ---
fwrite(plateau_stats, file.path(OUT_DIR, "SF1D_U_Plateau_Stats.csv"))
log_msg(sprintf("Plot successfully rendered and saved to: %s", out_pdf))