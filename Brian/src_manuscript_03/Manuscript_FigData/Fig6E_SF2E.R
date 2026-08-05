#!/usr/bin/env Rscript
# derived from: partD.05.grand_compare.R 
#
# Fig6E: Grand comparison for CasRx, and PspCas13b (whole)
# SF3E: Grand comparison for CasRx, and PspCas13b (3'UTR and CDS)

library(ggplot2)
library(data.table)
library(dplyr)

log_msg <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg), flush = TRUE)
}

# ================= 0. Configuration and Initialization =================
PROJECT <- "/data/cai801/data/HKUcas"
PART_A_DIR <- file.path(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
PART_B_DIR <- file.path(PROJECT, "result/08.manuscript_03/partB/02.grid_results")

# Unified output to the main manuscript figure directory
OUT_DIR <- file.path(PROJECT, "result/08.manuscript_03/Manuscript_FigData/main")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COR_METHOD <- "spearman"

# Define the finalized "two-parameter standards" from Part D
TARGET_PARAMS <- list(
  "CasRx_day5" = list(c(35, 20, 11), c(80, 40, 9)),
  "PspCas13b"  = list(c(35, 20, 11), c(80, 40, 17))
)

get_sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# ================= 1. Data Merging and Reconstruction =================
log_msg("Step 1: Loading Part A external models baseline...")
part_a_file <- file.path(PART_A_DIR, "05.self_other_compare.csv")
if (!file.exists(part_a_file)) stop("Part A results not found!")
df_external <- fread(part_a_file)

log_msg("Step 2: Extracting targeted Self-Models from Part B grids...")
self_models_list <- list()

for (tool_name in names(TARGET_PARAMS)) {
  grid_file <- file.path(PART_B_DIR, paste0(tool_name, "_grid.csv"))
  if (!file.exists(grid_file)) next
  
  df_grid <- fread(grid_file)
  params <- TARGET_PARAMS[[tool_name]]
  
  for (p in params) {
    w_val <- p[1]; l_val <- p[2]; u_val <- p[3]
    
    sub_df <- df_grid[W == w_val & L == l_val & U == u_val & direction == "3prime"]
    
    for (i in 1:nrow(sub_df)) {
      raw_region <- sub_df$region[i]
      clean_region <- ifelse(raw_region == "whole", "Whole", raw_region)
      if (!clean_region %in% c("Whole", "CDS", "3'UTR")) next
      
      self_models_list[[length(self_models_list) + 1]] <- data.table(
        Tool = tool_name,
        Region = clean_region,
        Method = sprintf("Self_W%d_L%d", w_val, l_val),
        Correlation = sub_df$spearman_r[i],
        P_value = sub_df$spearman_p[i],
        Sig_Label = get_sig_label(sub_df$spearman_p[i]),
        N = NA 
      )
    }
  }
}

df_self <- rbindlist(self_models_list)
df_grand <- rbind(df_external, df_self, fill = TRUE)

# Filter out tools that are no longer in TARGET_PARAMS (e.g., shRNAone from external baseline)
df_grand <- df_grand[Tool %in% names(TARGET_PARAMS)]

out_csv <- file.path(OUT_DIR, "Task5_Grand_Performance_Compare.csv")
fwrite(df_grand, out_csv)
log_msg(sprintf("Grand comparison stats saved to: %s", out_csv))

# ================= 2. Data Formatting and Color Mapping =================
log_msg("Step 3: Formatting data for plots...")

df_grand <- df_grand[!is.na(Correlation)]
df_grand$Tool <- factor(df_grand$Tool, levels = names(TARGET_PARAMS))
df_grand$Region <- factor(df_grand$Region, levels = c("Whole", "CDS", "3'UTR"))

methods <- unique(df_grand$Method)
self_methods <- grep("Self_", methods, value = TRUE)
other_methods <- setdiff(methods, self_methods)
df_grand$Method <- factor(df_grand$Method, levels = c(sort(self_methods), sort(other_methods)))

df_grand$label_y <- ifelse(df_grand$Correlation >= 0, 
                           df_grand$Correlation + 0.04, 
                           df_grand$Correlation - 0.04)

df_grand$vjust_val <- ifelse(df_grand$Correlation >= 0, 0, 1)

df_grand <- df_grand %>%
  mutate(Color_Group = case_when(
    grepl("W35_L20", Method) ~ "W35L20",
    grepl("W80_L40", Method) ~ "W80L40",
    TRUE ~ "External Models"
  ))

# ================= 3. Core Plotting Functions =================
generate_plot <- function(dt) {
  ggplot(dt, aes(x = Method, y = Correlation, fill = Color_Group)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.7) +
    # Apply global font rules: Arial, 8pt
    geom_text(aes(label = Sig_Label, y = label_y, vjust = vjust_val), 
              size = 8 / .pt, family = "Arial") +
    scale_fill_manual(values = c("W35L20" = "#E41A1C", 
                                 "W80L40" = "#377EB8", 
                                 "External Models" = "gray70")) +
    facet_grid(Region ~ Tool, scales = "free_x", space = "free_x") + 
    scale_y_continuous(limits = c(-0.15, 0.7), breaks = c(0, 0.25, 0.5)) +
    labs(y = sprintf("%s r", tools::toTitleCase(COR_METHOD)), x = NULL, fill = "Model Type") +
    # Completely remove face="bold", clean layout
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
      strip.text = element_text(size = 8, color = "black"),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.3),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8, margin = margin(r = 2)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(3, "mm"),
      legend.margin = margin(t = -5),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
}

# ================= 4. Generating Grouped Plots =================
log_msg("Step 4: Generating Split PDFs...")

# Figure 1: Whole Region
df_whole <- df_grand[Region == "Whole"]
p_whole <- generate_plot(df_whole)
pdf_whole_path <- file.path(OUT_DIR, "Task5_Grand_Compare_Whole.pdf")
ggsave(pdf_whole_path, plot = p_whole, width = 3, height = 2, device = cairo_pdf)
log_msg(sprintf("Saved Whole Region PDF: %s (Width: 3, Height: 2)", pdf_whole_path))

# Figure 2: CDS & 3'UTR Regions
df_sub <- df_grand[Region %in% c("CDS", "3'UTR")]
p_sub <- generate_plot(df_sub)
pdf_sub_path <- file.path(OUT_DIR, "Task5_Grand_Compare_Subregions.pdf")
ggsave(pdf_sub_path, plot = p_sub, width = 3, height = 3, device = cairo_pdf)
log_msg(sprintf("Saved CDS & 3'UTR Region PDF: %s (Width: 3, Height: 3)", pdf_sub_path))

log_msg("--- All Tasks Completed Successfully ---")