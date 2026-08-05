#!/usr/bin/env python3
"""
05.run_fareh2024.py — Fareh 2024 (Nature SMB) PspCas13b scoring on NM_000371.4.

This implements the exact heuristic scoring from the Nature SMB 2024 paper
(findScore.R) for 30 nt PspCas13b crRNAs. Not applicable to other guide lengths.

conda env: casrx_wei2023_web (has pandas; no special dependencies needed)
"""

import os
import sys
from datetime import datetime

PROJECT = "/data/cai801/data/HKUcas"
DATA_DIR = os.path.join(PROJECT, "data")
RESULT_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/02.original_pred")
FASTA = os.path.join(DATA_DIR, "NM_000371.4.fasta")
ACC = "NM_000371.4"
LEN = 30
OUT_CSV = os.path.join(RESULT_DIR, f"Fareh2024_{LEN}_{ACC}.csv")

os.makedirs(RESULT_DIR, exist_ok=True)


def log(msg):
    print(f"[{datetime.now():'%Y-%m-%d %H:%M:%S']} {msg}", flush=True)


def reverse_complement(seq):
    """Convert target RNA sequence to guide crRNA sequence."""
    mapping = str.maketrans("ACGTUacgtu", "TGCAATGCAA")
    return seq.translate(mapping)[::-1]


def calculate_fareh2024_score(seq):
    """
    Exact replica of Nature SMB 2024 findScore.R logic for 30 nt guides.
    """
    score = 0
    good_seq = "GGNNNNNNNNNNNNDDDNNNNNNNNNNNNN"
    bad_seq = "CCCCNNNNNNCCNNCCCHNNNNNNNNNNNN"

    iupac = {
        'G': {'G'}, 'C': {'C'}, 'A': {'A'}, 'T': {'T'},
        'N': {'A', 'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'}
    }

    # Positions 1-30
    for i in range(30):
        nt = seq[i]
        pos = i + 1

        # Favorable (Good)
        if nt in iupac.get(good_seq[i], set()):
            if pos in [1, 2]:
                score += 60
            if pos in [11, 12, 15, 16, 17]:
                score += 5

        # Penalised (Bad)
        if nt in iupac.get(bad_seq[i], set()):
            if pos in [1, 2]:
                score -= 60
            if pos in [3, 4, 11, 12, 15, 16, 17]:
                score -= 5

    # Secondary bonus/penalty at key positions 11, 12, 15, 16, 17
    for pos in [11, 12, 15, 16, 17]:
        if seq[pos - 1] == 'C':
            score -= 5
        else:
            score += 5

    return score


log("05. Fareh 2024 (Nature SMB): NM_000371.4")
log(f"    Guide length: {LEN} nt (native, designed for PspCas13b 30nt crRNA)")

# Read FASTA
with open(FASTA) as f:
    lines = f.readlines()
seq_flat = "".join(line.strip().upper().replace('U', 'T') for line in lines if not line.startswith(">"))

log(f"    Transcript length: {len(seq_flat)} nt")

records = []
# 30nt sliding window across full transcript
for i in range(len(seq_flat) - LEN + 1):
    target_sub = seq_flat[i:i + LEN]
    spacer = reverse_complement(target_sub)

    # Fatal structure filter: exclude sequences containing TTTT
    if "TTTT" in spacer:
        continue

    score = calculate_fareh2024_score(spacer)

    records.append({
        'transcript_id': ACC,
        'position': i + 1,
        'guide': spacer,
        'Nature2024_Score': score
    })

import pandas as pd
df = pd.DataFrame(records)
if not df.empty:
    df = df.sort_values(by=['Nature2024_Score'], ascending=False)
    df['rank'] = df['Nature2024_Score'].rank(method='min', ascending=False)
    df.to_csv(OUT_CSV, index=False)
    log(f"    {len(df)} compliant guides predicted -> {os.path.basename(OUT_CSV)}")
else:
    log("    WARNING: No compliant guides found (all contained TTTT?)")

log("05. Done.")
