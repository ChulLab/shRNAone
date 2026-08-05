#!/usr/bin/env python3
"""
13.offset_wet_dry.py — Offset scanning and data preparation.

For each wet tool (CasRx_day5, shRNAone, PspCas13b):
  1. Scan offsets from -50 to +50 nt using ALIGN_BY coordinate.
  2. Compute correlation (COR_METHOD) at each offset.
  3. Find the optimal offset based on MAXIMUM REAL CORRELATION (r > best_r, ignoring negative peaks).
  4. Export 2 CSVs per tool:
     - offset_wet_dry_...csv: Correlations across all offsets.
     - compare_wet_dry_...csv: Summary stats at the BEST offset.
"""

import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime
from scipy import stats

PROJECT = "/data/cai801/data/HKUcas"
WET_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
MODEL_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
FIGURE_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/04.offset_cor")
os.makedirs(FIGURE_DIR, exist_ok=True)

def log(msg):
    print(f"[{datetime.now():'%Y-%m-%d %H:%M:%S']} {msg}", flush=True)

def compute_correlation(x, y, method="spearman"):
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean, y_clean = x[mask], y[mask]
    if len(x_clean) < 3:
        return np.nan, np.nan
    if method == "spearman":
        r, p = stats.spearmanr(x_clean, y_clean)
    else:
        r, p = stats.pearsonr(x_clean, y_clean)
    return r, p

# ================= Configuration =================
COR_METHOD = "spearman"  # Options: "spearman" or "pearson"
ALIGN_BY = "end_pos"     # Options: "end_pos" or "start_pos"
OFFSET_RANGE = 50

TOOL_CONFIGS = [
    {
        "name": "CasRx_day5",
        "wet_file": "TTR_CasRx_day5_clean.csv",
        "log2fc_col": "mean_log2FC",
        "lg": 23,
        "dry_models": ["DeepCas13_23", "Hsu2023_30", "Sanjana2020_23", "TIGER2023_23"],
        "prefix": "02"
    },
    {
        "name": "shRNAone",
        "wet_file": "TTR_shRNAone_clean.csv",
        "log2fc_col": "mean_log2FC",
        "lg": 18,
        "dry_models": ["DeepCas13_18", "Hsu2023_30", "Sanjana2020_23", "TIGER2023_23"],
        "prefix": "03"
    },
    {
        "name": "PspCas13b",
        "wet_file": "TTR_PspCas13b_clean.csv",
        "log2fc_col": "mean_log2FC",
        "lg": 30,
        "dry_models": ["Fareh2024_30", "DeepCas13_30", "Hsu2023_30", "Sanjana2020_30", "TIGER2023_23"],
        "prefix": "04"
    }
]

log(f"04. Offset scanning: wet vs dry (method: {COR_METHOD}, align_by: {ALIGN_BY})")

for tc in TOOL_CONFIGS:
    tool_name = tc["name"]
    prefix = tc["prefix"]
    log(f"\n  === {tool_name} ===")

    # Load wet data
    wet_path = os.path.join(WET_DIR, tc["wet_file"])
    if not os.path.exists(wet_path):
        log(f"    [SKIP] Wet file not found: {wet_path}")
        continue
        
    df_wet = pd.read_csv(wet_path)
    if ALIGN_BY not in df_wet.columns:
        log(f"    [ERROR] {ALIGN_BY} not found in wet data.")
        continue
        
    log(f"    Wet data: {len(df_wet)} targets")

    wet_lookup = {}
    for _, row in df_wet.iterrows():
        pos = int(row[ALIGN_BY])
        wet_lookup[pos] = float(row[tc["log2fc_col"]])

    offset_results = {}
    best_stats = {}

    for model_name in tc["dry_models"]:
        model_file = os.path.join(MODEL_DIR, f"{model_name}_aligned.csv")
        if not os.path.exists(model_file):
            continue

        df_model = pd.read_csv(model_file)
        if ALIGN_BY not in df_model.columns:
            continue
            
        model_lookup = {}
        for _, row in df_model.iterrows():
            pos = int(row[ALIGN_BY])
            if pos not in model_lookup:
                model_lookup[pos] = []
            model_lookup[pos].append(float(row["prediction_score"]))

        for p_idx in model_lookup:
            model_lookup[p_idx] = np.mean(model_lookup[p_idx])

        # Scan offsets
        offset_scores = {}
        best_offset = None
        best_r = -np.inf  # Track highest positive correlation
        best_p = np.nan
        best_n = 0

        for offset in range(-OFFSET_RANGE, OFFSET_RANGE + 1):
            wet_vals, dry_vals = [], []
            for p_idx, log2fc in wet_lookup.items():
                shifted_p = p_idx + offset
                if shifted_p in model_lookup:
                    wet_vals.append(log2fc)
                    dry_vals.append(model_lookup[shifted_p])

            if len(wet_vals) >= 3:
                r, p = compute_correlation(np.array(wet_vals), np.array(dry_vals), COR_METHOD)
                offset_scores[offset] = (r, p, len(wet_vals))

                # STRICT RULE: Must be the absolute highest real value (r > best_r)
                if not np.isnan(r) and r > best_r:
                    best_r = r
                    best_offset = offset
                    best_p = p
                    best_n = len(wet_vals)

        offset_results[model_name] = offset_scores
        best_stats[model_name] = {"offset": best_offset, "r": best_r, "p": best_p, "n": best_n}
        log(f"    {model_name}: best offset = {best_offset}, r = {best_r:.4f}")

    # Save A: Offset scan data
    offset_csv_rows = []
    for model_name, scores in offset_results.items():
        for offset, (r, p, n) in sorted(scores.items()):
            offset_csv_rows.append({"model": model_name, "offset": offset, "correlation_r": r, "correlation_p": p, "n_points": n})
    pd.DataFrame(offset_csv_rows).to_csv(os.path.join(FIGURE_DIR, f"{prefix}.offset_wet_dry_{tool_name}.csv"), index=False)

    # Save B: Compare Stats (Changed spearman_r to generic cor_r)
    compare_rows = []
    for model_name, m_stats in best_stats.items():
        compare_rows.append({
            "tool": tool_name, "dry_model": model_name, "best_offset": m_stats["offset"],
            "cor_r": m_stats["r"], "p_value": m_stats["p"], "n_targets": m_stats["n"]
        })
    pd.DataFrame(compare_rows).to_csv(os.path.join(FIGURE_DIR, f"{prefix}.compare_wet_dry_{tool_name}.csv"), index=False)
    log(f"    Saved stat CSVs for {tool_name}")

log("\n04. Done.")