#!/usr/bin/env python3
"""
Extract clean wet-lab data.
Input:  /data/cai801/data/HKUcas/result/08.manuscript_03/partA/00.original_wet
Output: /data/cai801/data/HKUcas/result/08.manuscript_03/partA/01.clean_wet
"""

import os, sys, pandas as pd
from datetime import datetime

# 修改后的绝对路径
DATA_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partA/00.original_wet"
OUT_DIR  = "/data/cai801/data/HKUcas/result/08.manuscript_03/partA/01.clean_wet"
os.makedirs(OUT_DIR, exist_ok=True)

# 配置要处理的文件及对应的位点列名（NM ID 不同）
FILES = {
    "TTR_PspCas13b": {
        "fname": "Log2FC_TTR_PspCas13b_1nt.csv",
        "fc_col": "Log2FC_mean",
        "start_col": "start_nt_NM_000371.4",
        "end_col": "end_nt_NM_000371.4"
    },
    "TTR_shRNAone": {
        "fname": "Log2FC_TTR_shRNAone_1nt_260702_3reps.csv",
        "fc_col": "Log2FC_mean",
        "start_col": "start_nt_NM_000371.4",
        "end_col": "end_nt_NM_000371.4"
    },
    "TTR_CasRx_day5": {
        "fname": "Log2FC_TTR_RfxCas13d_U1-crRNA_1nt.csv",
        "fc_col": "rep3_day5_Log2FC",
        "start_col": "start_nt_NM_000371.4",
        "end_col": "end_nt_NM_000371.4"
    },
    "PCSK9_shRNAone_3nt": {
        "fname": "Log2FC_PCSK9_shRNAone_3nt_260706.csv",
        "fc_col": "Log2FC_mean",
        "start_col": "start_nt_NM_174936.3",
        "end_col": "end_nt_NM_174936.3"
    }
}

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def main():
    log("=== 01.3 Clean Final Wet Data ===")

    for label, cfg in FILES.items():
        fp = os.path.join(DATA_DIR, cfg["fname"])
        if not os.path.exists(fp):
            log(f"  [SKIP] {fp} not found")
            continue

        df = pd.read_csv(fp)

        # 提取指定的字段
        src_cols = [
            "name", 
            "target_seq", 
            cfg["start_col"], 
            cfg["end_col"], 
            "region", 
            cfg["fc_col"]
        ]
        
        # 检查列是否存在
        missing_cols = [c for c in src_cols if c not in df.columns]
        if missing_cols:
            log(f"  [ERROR] {cfg['fname']} 缺少以下列: {missing_cols}")
            continue

        df = df[src_cols].copy()

        # 过滤有效 region
        df = df[df["region"].isin(["5'UTR", "CDS", "3'UTR"])]

        # 列重命名映射规则：（name, target_seq, start_pos, end_pos, region, mean_log2FC）
        rename_dict = {
            cfg["start_col"]: "start_pos",
            cfg["end_col"]: "end_pos",
            cfg["fc_col"]: "mean_log2FC"
        }
        df = df.rename(columns=rename_dict)
        
        # 过滤 FC 缺失值
        df = df.dropna(subset=["mean_log2FC"])

        # 保证输出列的顺序符合规范
        target_order = ["name", "target_seq", "start_pos", "end_pos", "region", "mean_log2FC"]
        df = df[target_order]

        out_path = os.path.join(OUT_DIR, f"{label}_clean.csv")
        df.to_csv(out_path, index=False)
        log(f"  {label}: {len(df)} rows → {out_path}")
        log(f"    Regions: {df['region'].value_counts().to_dict()}")

    log("=== Done ===")

if __name__ == "__main__":
    main()