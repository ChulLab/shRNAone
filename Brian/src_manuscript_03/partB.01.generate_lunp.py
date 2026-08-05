#!/usr/bin/env python3
"""
01.generate_lunp.py — Generate RNAplfold _lunp files over (W, L) grid using multiprocessing.

Grid:
  W ∈ [20, 200]
  L ∈ [round(0.25*W), round(0.75*W)]
  u = 30
Output: /data/cai801/data/HKUcas/result/08.manuscript_03/partB/01.lunp_files/W{XXX}/lunp_W{XXX}_L{YYY}.txt
Cores: 40
"""

import os, sys, glob, shutil, subprocess, math
from datetime import datetime
from multiprocessing import Pool, cpu_count

# ================= Configuration =================
PROJECT = "/data/cai801/data/HKUcas"
BASE_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partB")

FASTA_DIR = os.path.join(BASE_DIR, "00.fasta")
LUNP_DIR  = os.path.join(BASE_DIR, "01.lunp_files")
os.makedirs(LUNP_DIR, exist_ok=True)

NUM_CORES = 40
U_VAL = 30
W_RANGE = range(20, 201)  # W from 20 to 200

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def get_fasta_info():
    """Find fasta file and extract transcript Accession ID."""
    fa_files = glob.glob(os.path.join(FASTA_DIR, "*.fasta")) + glob.glob(os.path.join(FASTA_DIR, "*.fa"))
    if not fa_files:
        sys.exit(f"[ERROR] No fasta file found in {FASTA_DIR}")
    fa_path = fa_files[0]
    
    with open(fa_path, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith(">"):
            acc = first_line[1:].split()[0]
        else:
            acc = os.path.basename(fa_path).split('.')[0]
            
    return fa_path, acc

def run_task(task):
    """Worker function for single (W, L) computation."""
    W, L, fasta_path, acc = task
    
    w_dir = os.path.join(LUNP_DIR, f"W{W:03d}")
    out_path = os.path.join(w_dir, f"lunp_W{W:03d}_L{L:03d}.txt")
    
    # Skip if already calculated and valid
    if os.path.exists(out_path) and os.path.getsize(out_path) > 50:
        return 0, 1  # done=0, skipped=1

    # Unique temporary directory per process to prevent File I/O collisions
    pid = os.getpid()
    tmpdir = os.path.join(LUNP_DIR, f"_tmp_p{pid}")
    os.makedirs(tmpdir, exist_ok=True)
    
    tmp_fa = os.path.join(tmpdir, f"{acc}.fa")
    if not os.path.exists(tmp_fa):
        shutil.copy(fasta_path, tmp_fa)

    try:
        cmd = ["RNAplfold", "-W", str(W), "-L", str(L), "-u", str(U_VAL)]
        with open(tmp_fa, "r") as stdin_f:
            subprocess.run(
                cmd,
                cwd=tmpdir,
                stdin=stdin_f,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )
            
        src = os.path.join(tmpdir, f"{acc}_lunp")
        if os.path.exists(src):
            os.makedirs(w_dir, exist_ok=True)
            shutil.move(src, out_path)
            
            # Clean up dotplot ps file
            dp = os.path.join(tmpdir, f"{acc}_dp.ps")
            if os.path.exists(dp):
                os.remove(dp)
            return 1, 0  # done=1, skipped=0
    except Exception as e:
        return 0, 0
    finally:
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)

    return 0, 0

def main():
    fasta_path, acc = get_fasta_info()
    
    # Generate list of tasks based on modified W, L bounds
    tasks = []
    for W in W_RANGE:
        # L in [0.25*W, 0.75*W] (正整数)
        l_min = max(1, math.ceil(0.25 * W))
        l_max = math.floor(0.75 * W)
        
        for L in range(l_min, l_max + 1):
            tasks.append((W, L, fasta_path, acc))

    total = len(tasks)
    log(f"=== 01.generate_lunp ({acc}) ===")
    log(f"FASTA: {fasta_path}")
    log(f"W range: [20, 200], L range: [0.25*W, 0.75*W], U fixed: {U_VAL}")
    log(f"Total tasks: {total} (W, L) pairs | Cores: {NUM_CORES}")

    t0 = datetime.now()
    
    # Parallel execution
    n_done, n_skip = 0, 0
    with Pool(processes=NUM_CORES) as pool:
        for done, skip in pool.imap_unordered(run_task, tasks, chunksize=10):
            n_done += done
            n_skip += skip
            completed = n_done + n_skip
            if completed % 1000 == 0 or completed == total:
                log(f" Progress: [{completed}/{total}] ({completed/total*100:.1f}%)")

    elapsed = (datetime.now() - t0).total_seconds()
    log(f"Done: {n_done} calculated, {n_skip} skipped. Elapsed: {elapsed:.1f}s")

if __name__ == "__main__":
    main()