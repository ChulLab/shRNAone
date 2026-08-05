#!/usr/bin/env python3
"""
partF.02.run_RNAplfold.py
基于 W80L40 和 W35L20 两种参数，为全转录组计算无配对概率。
"""

import os
import sys
import csv
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# ================= 配置区 =================
BASE_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partF"
CSV_FILE = os.path.join(BASE_DIR, "00.raw_data", "transcript_list.csv")
CACHE_DIR = os.path.join(BASE_DIR, "01.lunp_files")

# 设定的双参数配置 (统一设定 max_U 为 30 以涵盖 U11 和 U17)
PARAMS = [(80, 40), (35, 20)]
MAX_U = 30
RNAPLFOLD_BIN = "RNAplfold" 
MAX_WORKERS = 50

def log_msg(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def process_single_transcript(acc, seq):
    fasta_input = f">{acc}\n{seq}\n".encode('utf-8')
    status = 'SKIP'
    
    for W, L in PARAMS:
        target_lunp = os.path.join(CACHE_DIR, f"{acc}_W{W}_L{L}_lunp")
        if os.path.exists(target_lunp): continue
        
        cmd = [RNAPLFOLD_BIN, "-W", str(W), "-L", str(L), "-u", str(MAX_U)]
        try:
            subprocess.run(cmd, input=fasta_input, cwd=CACHE_DIR, 
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, 
                           timeout=300, check=True)
            os.rename(os.path.join(CACHE_DIR, f"{acc}_lunp"), target_lunp)
            status = 'SUCCESS'
        except Exception as e:
            return ('ERROR', acc, str(e))
        finally:
            ps_file = os.path.join(CACHE_DIR, f"{acc}_dp.ps")
            if os.path.exists(ps_file): os.remove(ps_file)
            
    return (status, acc, None)

def main():
    log_msg("=== 02.run_RNAplfold.py ===")
    os.makedirs(CACHE_DIR, exist_ok=True)
    
    transcripts = [(row["transcript_accession"], row["transcript_sequence"]) 
                   for row in csv.DictReader(open(CSV_FILE, "r", encoding="utf-8"))]
            
    total_tasks = len(transcripts)
    log_msg(f"Processing {total_tasks} transcripts with {MAX_WORKERS} workers...")
    
    stats = {'SUCCESS': 0, 'SKIP': 0, 'ERROR': 0}
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_acc = {executor.submit(process_single_transcript, acc, seq): acc for acc, seq in transcripts}
        
        for i, future in enumerate(as_completed(future_to_acc), 1):
            status, acc, err = future.result()
            stats[status] += 1
            if i % 1000 == 0 or i == total_tasks:
                log_msg(f"Progress: [{i}/{total_tasks}] | Success: {stats['SUCCESS']} | Skipped: {stats['SKIP']}")
            if status == 'ERROR':
                log_msg(f"  [FAIL] {acc}: {err}")

    log_msg("=== Done ===")

if __name__ == "__main__":
    main()