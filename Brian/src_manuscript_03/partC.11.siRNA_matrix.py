#!/usr/bin/env python3
# partC.11.siRNA_matrix.py
#
# 功能：
# 1. 提取所有候选 (W, L) 及 U 范围 (6-23) 下的 3'prime (target_end) accessibility。
# 2. 与 siRNA_continuous.csv 合并，生成大宽表 00.siRNA_raw_matrix.csv。
# 3. 计算 Continuous 模式下的 Spearman 相关性矩阵 (与 i_score)。
# 4. 计算 Binary 模式下的分类区分度 (AUC, Wilcoxon P-value, Cliff's Delta)。

import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.metrics import roc_auc_score
from datetime import datetime

# ================= 0. 目录与初始化 =================
BASE_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partC"
CLEAN_DIR = os.path.join(BASE_DIR, "03.GT_siRNA_acc", "00.clean_siRNA")
LUNP_DIR = os.path.join(BASE_DIR, "03.GT_siRNA_acc", "01.lunp_files")
MATRIX_DIR = os.path.join(BASE_DIR, "03.GT_siRNA_acc", "02.siRNA_matrix")
CAND_FILE = os.path.join(BASE_DIR, "01.Self_WL_list", "03.selected_WL_Candidates.csv")

os.makedirs(MATRIX_DIR, exist_ok=True)

U_RANGE = range(6, 17)  # U = 6 to 16 (inclusive)

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

# ================= 1. 构建 siRNA 原始特征矩阵 =================
log("Step 1: Extracting accessibility matrix for WLU combinations...")

# 加载基础信息
cands = pd.read_csv(CAND_FILE)
wl_pairs = cands[["W", "L"]].drop_duplicates()
sirna_cont = pd.read_csv(os.path.join(CLEAN_DIR, "siRNA_continuous.csv"))
sirna_bin = pd.read_csv(os.path.join(CLEAN_DIR, "siRNA_binary.csv"))

genes = sirna_cont["Accession"].unique()
acc_dict = {row["ID"]: {} for _, row in sirna_cont.iterrows()}

# 遍历 W, L 参数组合
for _, row in wl_pairs.iterrows():
    W = int(row["W"])
    L = int(row["L"])
    log(f"  Processing W={W}, L={L}...")
    
    for gene in genes:
        job_id = f"{gene}_W{W}_L{L}"
        lunp_file = os.path.join(LUNP_DIR, f"{job_id}_lunp")
        
        if not os.path.exists(lunp_file):
            continue
            
        # 跳过文件头，解析 plfold 输出
        # lunp 文件格式：第一列为碱基位置，后续列为不同 U 值下的概率
        lp = pd.read_csv(lunp_file, skiprows=2, sep=r"\s+", header=None)
        if lp.shape[1] < max(U_RANGE) + 1:
            continue
            
        pos_array = lp.iloc[:, 0].values.astype(int)
        pos2idx = {int(p): i for i, p in enumerate(pos_array)}
        
        # 找到属于当前 gene 的所有 siRNA
        subset = sirna_cont[sirna_cont["Accession"] == gene]
        
        for _, sr in subset.iterrows():
            target_end = int(sr["target_end"])
            row_idx = pos2idx.get(target_end, -1)
            
            if row_idx < 0:
                continue
                
            sid = sr["ID"]
            # 提取 U=6 到 23 的可及性值
            for U in U_RANGE:
                val = lp.iloc[row_idx, U]
                if not pd.isna(val):
                    col_name = f"W_{W}_L_{L}_U_{U}"
                    acc_dict[sid][col_name] = float(val)

# 转换为 DataFrame 并合并
acc_df = pd.DataFrame.from_dict(acc_dict, orient="index")
acc_df.index.name = "ID"
acc_df = acc_df.reset_index()

raw_matrix = pd.merge(sirna_cont, acc_df, on="ID", how="left")
out_raw = os.path.join(MATRIX_DIR, "00.siRNA_raw_matrix.csv")
raw_matrix.to_csv(out_raw, index=False)
log(f"Raw matrix saved to {out_raw} with shape {raw_matrix.shape}")


# ================= 2. 计算 Continuous 相关性矩阵 =================
log("Step 2: Calculating Continuous Evaluation (Spearman correlation)...")

param_cols = [c for c in raw_matrix.columns if c.startswith("W_")]
cont_results = []

for col in param_cols:
    df_clean = raw_matrix[["i_score", col]].dropna()
    if len(df_clean) < 20:
        continue
        
    r, p = spearmanr(df_clean[col], df_clean["i_score"])
    
    parts = col.split("_")
    cont_results.append({
        "Parameter": col,
        "W": int(parts[1]),
        "L": int(parts[3]),
        "U": int(parts[5]),
        "Spearman_r": r,
        "P_value": p,
        "Valid_N": len(df_clean)
    })

cont_df = pd.DataFrame(cont_results).sort_values(by="Spearman_r", ascending=False)
out_cont = os.path.join(MATRIX_DIR, "01.Evaluation_Continuous.csv")
cont_df.to_csv(out_cont, index=False)
log(f"Continuous evaluation matrix saved to {out_cont}")


# ================= 3. 计算 Binary 区分度矩阵 =================
log("Step 3: Calculating Binary Evaluation (AUC & Statistics)...")

# 保留在 Binary 数据集中 (即去除了 intermediate 的极端组)
bin_matrix = pd.merge(sirna_bin[["ID", "group"]], raw_matrix, on="ID", how="left")
bin_results = []

for col in param_cols:
    df_clean = bin_matrix[["group", col]].dropna()
    
    func_vals = df_clean[df_clean["group"] == "functional"][col].values
    nonfunc_vals = df_clean[df_clean["group"] == "nonfunctional"][col].values
    
    if len(func_vals) < 5 or len(nonfunc_vals) < 5:
        continue
        
    # Wilcoxon Rank-Sum Test
    stat, p_val = mannwhitneyu(func_vals, nonfunc_vals, alternative='two-sided')
    
    # ROC AUC 计算 (functional 为正类)
    y_true = (df_clean["group"] == "functional").astype(int)
    y_score = df_clean[col]
    try:
        auc = roc_auc_score(y_true, y_score)
    except ValueError:
        auc = np.nan
        
    # Cliff's Delta
    # Delta = (P(x > y) - P(x < y))
    nx, ny = len(func_vals), len(nonfunc_vals)
    greater = sum(sum(x > y for y in nonfunc_vals) for x in func_vals)
    less = sum(sum(x < y for y in nonfunc_vals) for x in func_vals)
    cliff_delta = (greater - less) / (nx * ny) if (nx * ny) > 0 else np.nan
    
    parts = col.split("_")
    bin_results.append({
        "Parameter": col,
        "W": int(parts[1]),
        "L": int(parts[3]),
        "U": int(parts[5]),
        "AUC": auc,
        "Cliff_Delta": cliff_delta,
        "Wilcoxon_P": p_val,
        "N_Func": nx,
        "N_NonFunc": ny
    })

bin_df = pd.DataFrame(bin_results).sort_values(by="AUC", ascending=False)
out_bin = os.path.join(MATRIX_DIR, "02.Evaluation_Binary.csv")
bin_df.to_csv(out_bin, index=False)
log(f"Binary evaluation matrix saved to {out_bin}")

log("Module partC.11 finished successfully. Matrices are ready for visualization.")