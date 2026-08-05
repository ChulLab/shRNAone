#!/usr/bin/env python3
"""
04.run_deepcas13.py — DeepCas13 prediction on NM_000371.4 at multiple guide lengths.

DeepCas13 accepts a --length parameter for guide extraction from the target transcript.
All extracted guides are padded to 33 nt internally before CNN input, so any length works.
Guide lengths tested: 18, 23, 30 nt.

conda env: deepcas13
"""

import os
import sys
import subprocess
from datetime import datetime

PROJECT = "/data/cai801/data/HKUcas"
SOFTWARE = os.path.join(PROJECT, "software/deepcas13")
MODEL = os.path.join(SOFTWARE, "trained_model_3060")
RESULT_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/02.original_pred")
FASTA = os.path.join(PROJECT, "data", "NM_000371.4.fasta")
ACC = "NM_000371.4"
LENGTHS = [18, 23, 30]

os.makedirs(RESULT_DIR, exist_ok=True)


def log(msg):
    print(f"[{datetime.now():'%Y-%m-%d %H:%M:%S']} {msg}", flush=True)


# Flatten FASTA to single-line for DeepCas13
with open(FASTA) as f:
    lines = f.readlines()
flat_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
temp_fa = os.path.join(RESULT_DIR, f"temp_{ACC}.fasta")
with open(temp_fa, "w") as f:
    f.write(f">{ACC}\n{flat_seq}\n")

log("04. DeepCas13: NM_000371.4")
log(f"    Testing guide lengths: {LENGTHS}")
log(f"    Note: All guides padded to 33 nt internally before CNN input.")

for L in LENGTHS:
    out_csv = os.path.join(RESULT_DIR, f"DeepCas13_{L}_{ACC}.csv")
    cmd = [
        sys.executable, os.path.join(SOFTWARE, "deepcas13.py"),
        "--seq", temp_fa,
        "--model", MODEL,
        "--type", "target",
        "--length", str(L),
        "--output", out_csv
    ]
    log(f"    Running --length {L} ...")
    proc = subprocess.run(cmd, cwd=SOFTWARE, capture_output=True, text=True)
    if proc.returncode == 0 and os.path.exists(out_csv) and os.path.getsize(out_csv) > 100:
        log(f"    SUCCESS with guide length {L} -> {os.path.basename(out_csv)}")
    else:
        log(f"    FAILED with guide length {L}: {proc.stderr.strip()[:200]}")

if os.path.exists(temp_fa):
    os.remove(temp_fa)

log("04. Done.")
