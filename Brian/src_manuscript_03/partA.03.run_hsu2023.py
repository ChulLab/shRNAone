#!/usr/bin/env python3
"""
03.run_hsu2023.py — Hsu 2023 Web prediction on NM_000371.4 at 30 nt guide length.

The Hsu2023 web model (RNAtargeting) has a fixed guide length of 30 nt:
  - make_guide_library_features() in linearfold.py hardcodes 30nt extraction (trans_seq[i:i+30])
  - Keras model input shape is fixed at (30, 4) from training
  - The guidelength parameter in predict.py is unused during inference with saved_model

Therefore, this script simply runs the model at its native 30 nt length.
No hacking required — the model natively produces 30 nt guides.

conda env: casrx_wei2023_web
"""

import os
import sys
from datetime import datetime

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

PROJECT = "/data/cai801/data/HKUcas"
WEB_DIR = os.path.join(PROJECT, "software/RNAtargeting_web_custom")
DATA_DIR = os.path.join(PROJECT, "data")
RESULT_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/02.original_pred")
FASTA = os.path.join(DATA_DIR, "NM_000371.4.fasta")
ACC = "NM_000371.4"
LEN = 30
OUT_CSV = os.path.join(RESULT_DIR, f"Hsu2023_{LEN}_{ACC}.csv")

os.makedirs(RESULT_DIR, exist_ok=True)


def log(msg):
    print(f"[{datetime.now():'%Y-%m-%d %H:%M:%S']} {msg}", flush=True)


log("03. Hsu 2023 Web: NM_000371.4")
log(f"    Guide length: {LEN} nt (native, fixed by model architecture)")
log(f"    Note: linearfold.py hardcodes 30nt guide extraction.")
log(f"    Keras model input shape: (30, 4) — cannot be changed without retraining.")

os.chdir(WEB_DIR)
sys.path.insert(0, WEB_DIR)

from rnatargeting.predict import run_pred

with open(FASTA, "rb") as f:
    fasta_bytes = f.read()

log("    Running prediction ...")
pred_df = run_pred(fasta_bytes, outfile=OUT_CSV)

log(f"    {len(pred_df)} guides predicted -> {os.path.basename(OUT_CSV)}")
log("03. Done.")
