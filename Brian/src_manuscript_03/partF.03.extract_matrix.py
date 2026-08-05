#!/usr/bin/env python3
"""
partF.03.extract_matrix.py
专为 PspCas13b 设计：滑动 30nt 窗口，融合 W80L40_U17 和 W35L20_U11 的 Accessibility，
并为每个转录本内的靶点分别计算降序排名 (Rank)，最终按物理坐标顺序输出。
"""

import os
import sys
import csv
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# ================= 配置区 =================
BASE_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partF"
CSV_FILE = os.path.join(BASE_DIR, "00.raw_data", "transcript_list.csv")
CACHE_DIR = os.path.join(BASE_DIR, "01.lunp_files")
OUT_DIR = os.path.join(BASE_DIR, "02.target_matrix")

GUIDE_LEN = 30
MAX_WORKERS = max(1, multiprocessing.cpu_count() - 2)

def log_msg(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def get_region(pos, cds_start, cds_end):
    if pos < cds_start: return "5'UTR"
    if pos <= cds_end: return "CDS"
    return "3'UTR"

def parse_lunp(filepath, target_U):
    """解析单个 lunp 文件，提取特定 U 的值。返回 {end_pos: accessibility}"""
    acc_dict = {}
    if not os.path.exists(filepath): return acc_dict
    
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", ">")): continue
            parts = line.split()
            if len(parts) > target_U and parts[target_U] != "NA":
                acc_dict[int(parts[0])] = float(parts[target_U])
    return acc_dict

def process_transcript(row):
    acc = row["transcript_accession"]
    version = row["transcript_version"]
    seq = row["transcript_sequence"]
    
    try:
        cds_start, cds_end = int(row["CDS_start"]), int(row["CDS_end"])
    except ValueError:
        return ('ERROR', acc, "Invalid CDS")

    out_csv = os.path.join(OUT_DIR, f"{acc}.csv")
    if os.path.exists(out_csv):
        return ('SKIP', acc, 0)

    lunp_w80 = os.path.join(CACHE_DIR, f"{acc}_W80_L40_lunp")
    lunp_w35 = os.path.join(CACHE_DIR, f"{acc}_W35_L20_lunp")
    
    if not (os.path.exists(lunp_w80) and os.path.exists(lunp_w35)):
        return ('MISSING_LUNP', acc, 0)

    # 提取 W80L40_U17 和 W35L20_U11 的 3' 端可及性
    acc_w80_u17 = parse_lunp(lunp_w80, 17)
    acc_w35_u11 = parse_lunp(lunp_w35, 11)

    results = []
    max_start = len(seq) - GUIDE_LEN + 1 
    
    for start_pos in range(1, max_start + 1):
        end_pos = start_pos + GUIDE_LEN - 1 # 3' 末端坐标
        
        val_w80 = acc_w80_u17.get(end_pos, "NA")
        val_w35 = acc_w35_u11.get(end_pos, "NA")
        
        # 只保留至少在一个参数下有打分的靶点
        if val_w80 == "NA" and val_w35 == "NA":
            continue
            
        results.append({
            "transcript_version": version,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "target_sequence": seq[start_pos - 1 : end_pos],
            "target_region": get_region(end_pos, cds_start, cds_end),
            "acc_W80_L40_U17": val_w80,
            "acc_W35_L20_U11": val_w35
        })
        
    if results:
        # --- 计算 Rank 排行 ---
        
        # 1. 按照 W80_L40_U17 降序排名 (分数越高越好，NA 视为 -1.0 垫底)
        results.sort(key=lambda x: float(x["acc_W80_L40_U17"]) if x["acc_W80_L40_U17"] != "NA" else -1.0, reverse=True)
        for i, r in enumerate(results, 1):
            r["rank_W80_L40_U17"] = i if r["acc_W80_L40_U17"] != "NA" else "NA"
            
        # 2. 按照 W35_L20_U11 降序排名
        results.sort(key=lambda x: float(x["acc_W35_L20_U11"]) if x["acc_W35_L20_U11"] != "NA" else -1.0, reverse=True)
        for i, r in enumerate(results, 1):
            r["rank_W35_L20_U11"] = i if r["acc_W35_L20_U11"] != "NA" else "NA"
            
        # 3. 恢复成基于转录本 5' 到 3' 物理坐标排序
        results.sort(key=lambda x: x["start_pos"])

        # 写入 CSV
        fieldnames = [
            "transcript_version", "start_pos", "end_pos", "target_sequence", "target_region",
            "acc_W80_L40_U17", "rank_W80_L40_U17", 
            "acc_W35_L20_U11", "rank_W35_L20_U11"
        ]
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
            
        return ('SUCCESS', acc, len(results))
    
    return ('NO_TARGETS', acc, 0)

def main():
    log_msg("=== 03.extract_matrix_PspCas13b.py ===")
    os.makedirs(OUT_DIR, exist_ok=True)
    
    transcripts = list(csv.DictReader(open(CSV_FILE, "r", encoding="utf-8")))
    total = len(transcripts)
    log_msg(f"Loaded {total} transcripts. Starting {MAX_WORKERS} workers...")
    
    stats = {'SUCCESS': 0, 'SKIP': 0, 'MISSING_LUNP': 0, 'NO_TARGETS': 0, 'ERROR': 0}
    
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_transcript, row): row["transcript_accession"] for row in transcripts}
        
        for i, future in enumerate(as_completed(futures), 1):
            status, acc, count = future.result()
            stats[status] += 1
            if i % 1000 == 0 or i == total:
                log_msg(f"Progress: [{i}/{total}] | Success: {stats['SUCCESS']} | Skipped: {stats['SKIP']} | Missing LUNP: {stats['MISSING_LUNP']}")

    log_msg("=== Done ===")

if __name__ == "__main__":
    main()