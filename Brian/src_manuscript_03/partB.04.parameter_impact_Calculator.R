#!/usr/bin/env Rscript
# partB.04.parameter_impact_Calc.R
#
# Core objective: Quantify the relative influence of thermodynamic parameters 
# (W, L, U) using both linear and non-linear frameworks.
#
# Methods included:
# 1. ANOVA Variance Partitioning (Linear sum-of-squares)
# 2. Standardized Multiple Linear Regression (Linear marginal effect)
# 3. Random Forest via 'ranger' (Non-linear Permutation Importance, Multi-threaded)
# 4. Generalized Additive Models (GAM) (Statistical significance of non-linear smooth terms)

library(dplyr)
library(broom)
library(data.table)
library(ranger) # Ultra-fast multi-threaded Random Forest
library(mgcv)   # Generalized Additive Models

# ================= 0. Configuration & Environment =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
CACHE_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/04.parameter_impact/cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)

# ================= 1. Data Loading & Preprocessing =================
log_msg("Step 1: Loading and aggregating grid search data...")
all_data <- rbindlist(lapply(TOOLS, function(tool) {
  file_path <- file.path(IN_DIR, paste0(tool, "_grid.csv"))
  if (file.exists(file_path)) {
    df <- fread(file_path)
    df$tool <- tool
    return(df)
  }
}))

# Filter for physically meaningful performance (r > 0) and the optimal anchor (3prime)
study_data <- all_data[spearman_r > 0 & direction == "3prime"]
log_msg(sprintf("Global parameter space loaded successfully: %d rows retained.", nrow(study_data)))

# ================= 2. Core Statistical Calculation Loop =================
anova_res_list <- list()
lm_res_list <- list()
rf_res_list <- list()
gam_res_list <- list()

for (t in TOOLS) {
  log_msg(sprintf(">>> Processing Tool: %s", t))
  sub_df <- study_data[tool == t]
  
  if(nrow(sub_df) < 50) {
    log_msg(sprintf("    Skipping %s due to insufficient data points.", t))
    next
  }
  
  # ---------- 2.1 ANOVA Variance Decomposition ----------
  log_msg("    Running ANOVA Type I...")
  fit_aov <- aov(spearman_r ~ W + L + U + region, data = sub_df)
  ss <- summary(fit_aov)[[1]]$`Sum Sq`
  vars <- trimws(rownames(summary(fit_aov)[[1]]))
  valid_idx <- which(vars %in% c("W", "L", "U"))
  pct_ss <- ss[valid_idx] / sum(ss[valid_idx]) * 100
  
  anova_res_list[[t]] <- data.table(
    Tool = t, Parameter = vars[valid_idx], Variance_Explained_Pct = pct_ss
  )
  
  # ---------- 2.2 Standardized Linear Regression ----------
  log_msg("    Running Standardized Linear Regression...")
  sub_df_scaled <- copy(sub_df)
  sub_df_scaled[, c("W_sc", "L_sc", "U_sc") := .(scale(W), scale(L), scale(U))]
  
  fit_lm <- lm(spearman_r ~ W_sc + L_sc + U_sc + region, data = sub_df_scaled)
  lm_tidy <- setDT(tidy(fit_lm))[term %in% c("W_sc", "L_sc", "U_sc")]
  lm_tidy[, term := gsub("_sc", "", term)] # Clean up names back to W, L, U
  
  lm_res_list[[t]] <- data.table(
    Tool = t, Parameter = lm_tidy$term, Beta = lm_tidy$estimate, P_value = lm_tidy$p.value
  )
  
  # ---------- 2.3 Non-linear Feature Importance (Ranger) ----------
  log_msg("    Running Fast Random Forest (Permutation Importance)...")
  num_cores <- max(1, parallel::detectCores() - 1) # Leave one core free
  rf_model <- ranger(spearman_r ~ W + L + U, data = sub_df, 
                     num.trees = 500, 
                     importance = "permutation", 
                     scale.permutation.importance = TRUE,
                     num.threads = num_cores)
  
  imp <- rf_model$variable.importance
  rf_res_list[[t]] <- data.table(
    Tool = t, Parameter = names(imp), IncMSE = imp
  )
  
  # ---------- 2.4 Generalized Additive Model (GAM) ----------
  # Fits non-linear smooth splines s() to W, L, U to formally test their non-linear significance.
  log_msg("    Running Generalized Additive Model (GAM)...")
  fit_gam <- gam(spearman_r ~ s(W) + s(L) + s(U) + region, data = sub_df, method = "REML")
  gam_sum <- summary(fit_gam)
  
  # Extract F-statistic and p-value for the smooth terms
  s_table <- as.data.frame(gam_sum$s.table)
  gam_res_list[[t]] <- data.table(
    Tool = t, 
    Parameter = gsub("s\\(|\\)", "", rownames(s_table)), # Extract "W" from "s(W)"
    F_value = s_table$F, 
    P_value = s_table$`p-value`
  )
}

# ================= 3. Data Persistence (Caching) =================
log_msg("Step 3: Saving computed statistical results to cache...")
fwrite(rbindlist(anova_res_list), file.path(CACHE_DIR, "01_ANOVA_res.csv"))
fwrite(rbindlist(lm_res_list), file.path(CACHE_DIR, "02_LM_res.csv"))
fwrite(rbindlist(rf_res_list), file.path(CACHE_DIR, "03_RF_res.csv"))
fwrite(rbindlist(gam_res_list), file.path(CACHE_DIR, "05_GAM_res.csv"))
saveRDS(study_data, file.path(CACHE_DIR, "04_Trajectory_data.rds"))

log_msg("--- All Calculations Completed Successfully! ---")