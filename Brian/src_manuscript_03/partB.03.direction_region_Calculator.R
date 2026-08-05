#!/usr/bin/env Rscript
# partB.03.direction_region_ComputeCore.R — High-Performance Statistical Engine
# 
# Statistical Methods & Additions:
# 1. FDR Correction (BH): Controls false discovery rate.
# 2. ECDF & K-S Test: Compares the stochastic dominance of distributions.
# 3. [NEW] Effect Size (Cohen's d): Quantifies the exact magnitude of performance 
#    differences, immune to the p-value inflation common in massive grid datasets.
# 4. [NEW] Variance Partitioning (Eta-squared): Evaluates the relative importance 
#    by measuring the percentage of variance uniquely explained by 'direction'.
# 5. LMM (Optimized): Isolates fixed effects using nloptwrap for fast convergence.
# 6. Top 5% Hypergeometric Enrichment: Orthogonal representation testing.
#
# Execution Architecture:
# - Strictly headless computation.
# - Leverages data.table for C-level memory efficiency.
# - Leverages foreach/doParallel for 80-core symmetric multiprocessing.
# - Saves all intermediate checkpoints as compressed .rds and .csv files.

library(data.table)
library(lme4)
library(lmerTest)
library(emmeans)
library(parallel)
library(doParallel)

# ================= 0. Configuration & Caching Setup =================
PROJECT <- "/data/cai801/data/HKUcas"
IN_DIR  <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/03.direction_region/cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOOLS <- c("shRNAone", "PspCas13b", "CasRx_day5")
NUM_CORES <- 80

# Configure multithreading for data.table and parallel loops
setDTthreads(NUM_CORES)
registerDoParallel(cores = NUM_CORES)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 1. High-Speed Data Loading & FDR Correction =================
log_msg("Step 1: Loading Grid Data via data.table...")

dt_list <- lapply(TOOLS, function(t) {
  file_path <- file.path(IN_DIR, paste0(t, "_grid.csv"))
  if (file.exists(file_path)) {
    dt <- fread(file_path, nThread = NUM_CORES)
    dt[, tool := t]
    return(dt)
  }
})
all_data <- rbindlist(dt_list, fill = TRUE)

# Filter positive space and calculate FDR
all_data <- all_data[spearman_r > 0]
all_data[, fdr_spearman := p.adjust(spearman_p, method = "BH")]

# Filter significant data
sig_data <- all_data[fdr_spearman < 0.05]

# Factorize to ensure 3prime is the contrast base
sig_data[, direction := factor(direction, levels = c("3prime", "5prime", "upstream", "downstream"))]
sig_data[, region := factor(region, levels = c("whole", "CDS", "3'UTR"))]

# Save checkpoint
saveRDS(sig_data, file.path(OUT_DIR, "01_sig_data_cache.rds"))
log_msg(sprintf("Checkpoint 1 Saved. Valid significant points: %d", nrow(sig_data)))

# ================= 2. Parallelized K-S Tests & Effect Size (Cohen's d) =================
log_msg("Step 2: Executing Parallel K-S Tests and computing Cohen's d...")

# Create task grid for parallelization
task_grid <- expand.grid(Tool = TOOLS, Region = unique(sig_data$region), stringsAsFactors = FALSE)

ks_results <- foreach(i = 1:nrow(task_grid), .combine = rbind, .packages = c("data.table")) %dopar% {
  t <- task_grid$Tool[i]
  r <- task_grid$Region[i]
  
  sub_dt <- sig_data[tool == t & region == r]
  if (nrow(sub_dt) < 50) return(NULL)
  
  base_dist <- sub_dt[direction == "3prime", spearman_r]
  res_list <- list()
  
  for (d in c("5prime", "upstream", "downstream")) {
    comp_dist <- sub_dt[direction == d, spearman_r]
    
    if (length(base_dist) > 10 && length(comp_dist) > 10) {
      # 1. K-S Test
      ks_test <- ks.test(base_dist, comp_dist, alternative = "less")
      
      # 2. Cohen's d Effect Size (Magnitude of difference)
      # d = (mean(base) - mean(comp)) / pooled_sd
      n1 <- length(base_dist); n2 <- length(comp_dist)
      var1 <- var(base_dist);  var2 <- var(comp_dist)
      pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
      cohen_d <- (mean(base_dist) - mean(comp_dist)) / pooled_sd
      
      res_list[[length(res_list) + 1]] <- data.table(
        Tool = t, Region = r, Base = "3prime", Contrast = d,
        D_statistic = ks_test$statistic, P_value = ks_test$p.value,
        Effect_Size_d = cohen_d
      )
    }
  }
  rbindlist(res_list)
}

fwrite(ks_results, file.path(OUT_DIR, "02_KS_EffectSize_Results.csv"))
log_msg("Checkpoint 2 Saved. K-S Tests & Effect Sizes computed.")

# ================= 3. [NEW] Variance Partitioning (Eta-Squared) =================
# Measures the proportion of total variance explained by direction vs region
log_msg("Step 3: Executing Variance Partitioning (ANOVA Eta-Squared)...")

eta_results <- foreach(t = TOOLS, .combine = rbind, .packages = c("data.table")) %dopar% {
  sub_dt <- sig_data[tool == t]
  if (nrow(sub_dt) < 100) return(NULL)
  
  # Standard ANOVA
  aov_model <- aov(spearman_r ~ direction * region, data = sub_dt)
  
  #  Sum Sq: 1=direction, 2=region, 3=direction:region, 4=Residuals
  ss <- summary(aov_model)[[1]][, "Sum Sq"]
  ss_total <- sum(ss)
  eta_sq <- ss / ss_total
  
  data.table(
    Tool = t,
    Variance_Explained_by_Direction = eta_sq[1] * 100,
    Variance_Explained_by_Region    = eta_sq[2] * 100,
    Variance_Explained_by_Interaction = eta_sq[3] * 100
  )
}

fwrite(eta_results, file.path(OUT_DIR, "03_Variance_Partitioning_EtaSq.csv"))
log_msg("Checkpoint 3 Saved. Variance Partitioning complete.")

# ================= 4. Optimized Linear Mixed-Effects Model (LMM) =================
log_msg("Step 4: Constructing LMM with nloptwrap optimization...")

# LMM is notoriously slow on massive datasets. We use 'nloptwrap' optimizer 
# and disable derivative calculations to massively accelerate convergence.
lmer_ctrl <- lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE)

lmm_results <- foreach(t = TOOLS, .combine = rbind, .packages = c("data.table", "lme4", "lmerTest", "emmeans")) %dopar% {
  sub_dt <- sig_data[tool == t]
  if (nrow(sub_dt) < 100) return(NULL)
  
  model <- lmer(spearman_r ~ direction * region + (1|W) + (1|L) + (1|U), 
                data = sub_dt, control = lmer_ctrl)
                
  emm_dir <- emmeans(model, pairwise ~ direction | region)
  contrasts <- as.data.table(emm_dir$contrasts)
  contrasts[, Tool := t]
  
  return(contrasts)
}

fwrite(lmm_results, file.path(OUT_DIR, "04_LMM_Contrasts.csv"))
log_msg("Checkpoint 4 Saved. LMM estimates extracted.")

# ================= 5. Parallelized Top 5% Hypergeometric Enrichment =================
log_msg("Step 5: Calculating Hypergeometric Enrichment...")

enrich_results <- foreach(i = 1:nrow(task_grid), .combine = rbind, .packages = c("data.table")) %dopar% {
  t <- task_grid$Tool[i]
  r <- task_grid$Region[i]
  
  sub_dt <- sig_data[tool == t & region == r]
  if (nrow(sub_dt) == 0) return(NULL)
  
  threshold <- quantile(sub_dt$spearman_r, 0.95)
  top_dt <- sub_dt[spearman_r >= threshold]
  
  res_list <- list()
  for (d in unique(sub_dt$direction)) {
    q <- nrow(top_dt[direction == d]) 
    m <- nrow(sub_dt[direction == d]) 
    n <- nrow(sub_dt) - m                        
    k <- nrow(top_dt)                            
    
    p_val <- phyper(q - 1, m, n, k, lower.tail = FALSE)
    
    res_list[[length(res_list) + 1]] <- data.table(
      Tool = t, Region = r, Direction = d,
      Top_Count = q, Total_in_Top = k,
      Enrichment_Ratio = (q / k) / (m / (m + n)),
      P_value = p_val
    )
  }
  rbindlist(res_list)
}

fwrite(enrich_results, file.path(OUT_DIR, "05_Top5_Enrichment_Stats.csv"))
log_msg("Checkpoint 5 Saved. Hypergeometric enrichment complete.")

# ================= Cleanup =================
stopImplicitCluster()
log_msg("=== Compute Core Execution Complete. Data cached for rapid plotting. ===")