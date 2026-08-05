#!/usr/bin/env python3
"""
02.grid_search.py — Ultra-fast parallelized grid search over (W, L, U, direction, region).

Features:
  - Multiprocessing execution across W/L files using 40 CPU cores.
  - NumPy matrix vectorization for blazingly fast Pearson & Spearman computations.
  - Direct alignment by start_pos / end_pos coordinates.
  - Computes both Spearman and Pearson (r, p) across regions: [whole, CDS, 3'UTR].

Output: /data/cai801/data/HKUcas/result/08.manuscript_03/partB/02.grid_results/
"""

import os, sys, glob, pandas as pd, numpy as np
from datetime import datetime
from scipy.stats import spearmanr, pearsonr
from multiprocessing import Pool, cpu_count

# ================= Configuration =================
PROJECT = "/data/cai801/data/HKUcas"
LUNP_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partB/01.lunp_files")
OUT_DIR  = os.path.join(PROJECT, "result/08.manuscript_03/partB/02.grid_results")
WET_DIR  = os.path.join(PROJECT, "result/08.manuscript_03/partA/01.clean_wet")
os.makedirs(OUT_DIR, exist_ok=True)

NUM_CORES = 40

# 提取 TTR 的 3 个湿实验工具及 guide 长度
TOOLS = {
    "shRNAone":   {"wet_file": "TTR_shRNAone_clean.csv",   "L_guide": 18},
    "PspCas13b":  {"wet_file": "TTR_PspCas13b_clean.csv",  "L_guide": 30},
    "CasRx_day5": {"wet_file": "TTR_CasRx_day5_clean.csv", "L_guide": 23},
}

DIRECTIONS = ["3prime", "5prime", "upstream", "downstream"]
REGIONS    = ["whole", "CDS", "3'UTR"]

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def load_wet_data(path):
    """Load wet-lab data and rename columns for unified lookup."""
    df = pd.read_csv(path)
    # 统一位点: 以 start_pos 为准 (或 end_pos)
    df = df.rename(columns={"start_pos": "pos", "mean_log2FC": "log2FC"})
    df = df[["pos", "log2FC", "region"]].dropna(subset=["log2FC"])
    df["pos"] = df["pos"].astype(int)
    return df

def process_single_file(args):
    """Worker task processing a single lunp_Wxxx_Lxxx.txt file against all tools."""
    fpath, W, L_val, wet_data_dict = args
    results = []

    try:
        # 显式指定 na_values="NA"
        lunp = pd.read_csv(
            fpath, 
            skiprows=2, 
            sep=r"\s+", 
            header=None, 
            na_values="NA", 
            engine='c'
        )
    except Exception:
        return results

    n_cols = lunp.shape[1]
    if n_cols < 3:
        return results

    # 第 0 列为 1-based 转录本 3' 端位置 i$，第 1~end 列为 u=1..30 的可及性概率
    pos_arr = lunp.iloc[:, 0].values.astype(int)
    acc_matrix = lunp.iloc[:, 1:].values.astype(np.float64) # [N_positions, Max_U]
    
    # 建立 1-based 快速哈希定位表
    pos_max = int(pos_arr.max()) + 200
    pos2idx = np.full(pos_max, -1, dtype=int)
    valid_mask = pos_arr < pos_max
    pos2idx[pos_arr[valid_mask]] = np.where(valid_mask)[0]

    max_u_file = acc_matrix.shape[1]

    for tool_name, tool_info in wet_data_dict.items():
        wet_df = tool_info["df"]
        Lg = tool_info["L_guide"]
        max_U = min(Lg, max_u_file, 30)
        if max_U < 2:
            continue

        pos_w = wet_df["pos"].values     # 1-based start_pos
        y_w = wet_df["log2FC"].values
        reg_w = wet_df["region"].values

        for U in range(2, max_U + 1):
            u_col = U - 1  # 0-indexed column in acc_matrix (1-based u=U)

            for direction in DIRECTIONS:
                if direction == "3prime":
                    tgt = pos_w + Lg - 1
                elif direction == "5prime":
                    tgt = pos_w + U - 1
                elif direction == "upstream":
                    tgt = pos_w - 1
                elif direction == "downstream":
                    tgt = pos_w + Lg - 1 + U

                # 快速安全映射
                in_bounds = (tgt >= 1) & (tgt < pos_max)
                idx = np.full(len(tgt), -1, dtype=int)
                idx[in_bounds] = pos2idx[tgt[in_bounds]]

                valid = idx >= 0
                if valid.sum() < 10:
                    continue

                acc_vals = acc_matrix[idx[valid], u_col]
                y_vals = y_w[valid]
                reg_vals = reg_w[valid]

                # 过滤 NA 和 NaN（因为开头 1~30 行会有 NA 变为 np.nan）
                finite_mask = np.isfinite(acc_vals) & np.isfinite(y_vals)
                if finite_mask.sum() < 10:
                    continue

                acc_clean = acc_vals[finite_mask]
                y_clean = y_vals[finite_mask]
                reg_clean = reg_vals[finite_mask]

                for r in REGIONS:
                    if r == "whole":
                        rmask = np.ones(len(reg_clean), dtype=bool)
                    else:
                        rmask = (reg_clean == r)

                    n_pts = rmask.sum()
                    if n_pts < 10:
                        continue

                    sub_acc = acc_clean[rmask]
                    sub_y = y_clean[rmask]

                    # 计算 Spearman 与 Pearson
                    sr, sp = spearmanr(sub_acc, sub_y)
                    pr, pp = pearsonr(sub_acc, sub_y)

                    results.append({
                        "tool": tool_name,
                        "W": W,
                        "L": L_val,
                        "U": U,
                        "direction": direction,
                        "region": r,
                        "spearman_r": float(sr),
                        "spearman_p": float(sp),
                        "pearson_r": float(pr),
                        "pearson_p": float(pp),
                        "n_pts": int(n_pts)
                    })

    return results

def main():
    log("=== 02.grid_search (Ultra-Fast Parallel Version) ===")

    # 预加载 3 个湿实验工具的数据
    wet_data_dict = {}
    for tool_name, cfg in TOOLS.items():
        wpath = os.path.join(WET_DIR, cfg["wet_file"])
        if not os.path.exists(wpath):
            log(f"  [SKIP] Wet file not found: {wpath}")
            continue
        wet_df = load_wet_data(wpath)
        wet_data_dict[tool_name] = {
            "df": wet_df,
            "L_guide": cfg["L_guide"]
        }
        log(f"  Loaded {tool_name}: {len(wet_df)} targets")

    if not wet_data_dict:
        sys.exit("[ERROR] No wet data loaded.")

    # 收集所有的 lunp 文件路径与参数
    tasks = []
    w_dirs = sorted([d for d in os.listdir(LUNP_DIR) if d.startswith("W") and os.path.isdir(os.path.join(LUNP_DIR, d))])
    
    for wd in w_dirs:
        W = int(wd[1:])
        w_path = os.path.join(LUNP_DIR, wd)
        for lf in os.listdir(w_path):
            if lf.startswith("lunp_") and lf.endswith(".txt"):
                # 解析 L 参数: lunp_W030_L015.txt
                parts = lf.replace(".txt", "").split("_")
                L_val = int(parts[2][1:])
                fpath = os.path.join(w_path, lf)
                tasks.append((fpath, W, L_val, wet_data_dict))

    total_files = len(tasks)
    log(f"Total lunp files to search: {total_files} across {len(wet_data_dict)} tools")
    log(f"Utilizing {NUM_CORES} CPU cores for parallel grid search...")

    t0 = datetime.now()
    results_by_tool = {t: [] for t in wet_data_dict}

    # 多进程并行极速计算
    with Pool(processes=NUM_CORES) as pool:
        for i, res_list in enumerate(pool.imap_unordered(process_single_file, tasks, chunksize=20), 1):
            for item in res_list:
                results_by_tool[item["tool"]].append(item)
            if i % 1000 == 0 or i == total_files:
                log(f"  Progress: [{i}/{total_files}] ({i/total_files*100:.1f}%)")

    # 分别保存每个 Tool 的结果 CSV
    for tool_name, rows in results_by_tool.items():
        if not rows:
            continue
        df_out = pd.DataFrame(rows)
        # 按照 (W, L, U) 排序
        df_out = df_out.sort_values(by=["region", "direction", "W", "L", "U"]).reset_index(drop=True)
        csv_path = os.path.join(OUT_DIR, f"{tool_name}_grid.csv")
        df_out.to_csv(csv_path, index=False)
        log(f"Saved: {csv_path} ({len(df_out)} rows)")

        # 输出各方向与区域的 Top 结果
        for d in DIRECTIONS:
            for r in REGIONS:
                sub = df_out[(df_out["direction"] == d) & (df_out["region"] == r)]
                if len(sub) > 0:
                    best = sub.loc[sub["spearman_r"].idxmax()]
                    log(f"  [{tool_name}] {d:10s} {r:7s}: W={int(best.W)} L={int(best.L)} U={int(best.U)} | Spearman r={best.spearman_r:.4f} (p={best.spearman_p:.2e}) | Pearson r={best.pearson_r:.4f}")

    elapsed = (datetime.now() - t0).total_seconds()
    log(f"=== Done in {elapsed:.1f}s ===")

if __name__ == "__main__":
    main()