#!/usr/bin/env Rscript
# partB.04.parameter_impact_Plotter.R
#
# Dedicated script for rendering highly condensed, publication-quality figures 
# based strictly on cached computational results.

library(ggplot2)
library(data.table)

# ================= 0. Configuration =================
PROJECT <- "/data/cai801/data/HKUcas"
CACHE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/04.parameter_impact/cache")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/04.parameter_impact/cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Muted & Minimalist Palette (Aligned with earlier Manuscript Figures)
PARAM_COLORS <- c("U" = "#D28E72", "W" = "#7A9E9F", "L" = "#B3BCA4")
BASE_FAMILY <- "Arial"
BASE_SIZE <- 9 # Optimal size for compact Supplementary Figures

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# Shared clean theme for bar plots
clean_bar_theme <- theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
    strip.text = element_text(color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.3),
    plot.title = element_text(size = BASE_SIZE + 1, face = "bold")
  )

# ================= 1. Plot ANOVA (Linear Variance) =================
log_msg("Plotting Fig 1: ANOVA Variance Decomposition...")
anova_res <- fread(file.path(CACHE_DIR, "01_ANOVA_res.csv"))

p_anova <- ggplot(anova_res, aes(x = Parameter, y = Variance_Explained_Pct, fill = Parameter)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  facet_wrap(~ Tool) +
  scale_fill_manual(values = PARAM_COLORS) +
  labs(title = "Linear Model: Variance Explained", y = "Explained Variance (%)", x = NULL) +
  clean_bar_theme

ggsave(file.path(OUT_DIR, "01.ANOVA_Variance_Plot.pdf"), plot = p_anova, 
       width = 6, height = 3.5, device = cairo_pdf)


# ================= 2. Plot Standardized LM (Linear Betas) =================
log_msg("Plotting Fig 2: Standardized Linear Regression...")
lm_res <- fread(file.path(CACHE_DIR, "02_LM_res.csv"))

p_lm <- ggplot(lm_res, aes(x = Parameter, y = Beta, fill = Parameter)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  facet_wrap(~ Tool) +
  scale_fill_manual(values = PARAM_COLORS) +
  labs(title = "Linear Model: Standardized Effect", y = "Standardized Beta Coefficient", x = NULL) +
  clean_bar_theme

ggsave(file.path(OUT_DIR, "02.Regression_Coefficients.pdf"), plot = p_lm, 
       width = 6, height = 3.5, device = cairo_pdf)


# ================= 3. Plot Random Forest (Non-linear Importance) =================
log_msg("Plotting Fig 3: Random Forest Feature Importance...")
rf_res <- fread(file.path(CACHE_DIR, "03_RF_res.csv"))

p_rf <- ggplot(rf_res, aes(x = reorder(Parameter, IncMSE), y = IncMSE, fill = Parameter)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  coord_flip() + # Flip coordinates for easier reading of importance ranking
  facet_wrap(~ Tool, scales = "free_x") +
  scale_fill_manual(values = PARAM_COLORS) +
  labs(title = "Non-Linear Model: RF Feature Importance", y = "Permutation Importance (IncMSE)", x = NULL) +
  clean_bar_theme

ggsave(file.path(OUT_DIR, "03.RF_Importance.pdf"), plot = p_rf, 
       width = 6.5, height = 3, device = cairo_pdf)


# ================= 4. Plot GAM (Non-linear Significance) =================
# Maps the F-statistic of the non-linear smooth splines. Higher F-value strongly
# supports the existence of complex, non-linear physical interactions.
log_msg("Plotting Fig 4: Generalized Additive Model (GAM) F-statistics...")
gam_res <- fread(file.path(CACHE_DIR, "05_GAM_res.csv"))

p_gam <- ggplot(gam_res, aes(x = reorder(Parameter, F_value), y = F_value, fill = Parameter)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  coord_flip() + 
  facet_wrap(~ Tool, scales = "free_x") +
  scale_fill_manual(values = PARAM_COLORS) +
  labs(title = "Non-Linear Model: GAM Smooth Term Significance", 
       y = "F-statistic of Non-linear Spline", x = NULL) +
  clean_bar_theme

ggsave(file.path(OUT_DIR, "05.GAM_F_statistic.pdf"), plot = p_gam, 
       width = 6.5, height = 3, device = cairo_pdf)


# ================= 5. Plot Marginal Trajectory =================
# log_msg("Plotting Fig 5: Marginal Trajectory...")
# study_data <- readRDS(file.path(CACHE_DIR, "04_Trajectory_data.rds"))

# p_traj <- ggplot(study_data, aes(x = U, y = spearman_r)) +
#   # Background individual trajectories (extremely transparent)
#   geom_line(aes(group = interaction(W, L)), alpha = 0.03, color = "black") +
#   # Main structural trend mapped using robust LOESS
#   geom_smooth(method = "loess", color = "#E41A1C", linewidth = 1.2, se = FALSE, span = 0.4) +
#   facet_grid(region ~ tool, scales = "free_y") +
#   labs(title = "Global Target Accessibility Landscape",
#        x = "Unpaired Length Constraint (U)", y = "Spearman Correlation (r)") +
#   theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
#   theme(
#     strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.3),
#     strip.text = element_text(color = "black", face = "bold"),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "black", linewidth = 0.3)
#   )

# ggsave(file.path(OUT_DIR, "04.Marginal_Effect_of_U.pdf"), plot = p_traj, 
#        width = 7, height = 5, device = cairo_pdf)

log_msg("--- All High-Density Plots Generated and Saved! ---")