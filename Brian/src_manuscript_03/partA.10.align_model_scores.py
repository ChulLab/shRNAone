#!/usr/bin/env python3
"""
12.align_model_scores.py — Align all model predictions to NM_000371.4 coordinates.

Takes score_result/*.csv files and normalizes them into a unified format:
  model_name, target_seq, start_pos, end_pos, prediction_score

Sorted by start_pos (ascending).

conda env: self_model
"""

import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime

PROJECT = "/data/cai801/data/HKUcas"
SCORE_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/02.original_pred")
RESULT_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/03.aligned_pred")
os.makedirs(RESULT_DIR, exist_ok=True)


def log(msg):
    print(f"[{datetime.now():'%Y-%m-%d %H:%M:%S']} {msg}", flush=True)


def reverse_complement(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def load_sanjana2020(filepath, length_label):
    """Sanjana2020: columns = GuideName, GuideSeq, MatchPos, GuideScores, ...
    
    MatchPos = end position (3' end) of the guide on the transcript (1-based).
    GuideSeq = crRNA sequence (already RC of target).
    Need to RC guide to get target_seq.
    """
    df = pd.read_csv(filepath)
    log(f"    Loaded Sanjana2020: {len(df)} rows")
    
    records = []
    for _, row in df.iterrows():
        guide = str(row["GuideSeq"]).strip()
        match_pos = int(row["MatchPos"])  # end_pos
        start_pos = match_pos - len(guide) + 1
        target_seq = reverse_complement(guide)
        score = float(row["GuideScores"])
        records.append({
            "model_name": f"Sanjana2020_{length_label}",
            "target_seq": target_seq.upper(),
            "start_pos": start_pos,
            "end_pos": match_pos,
            "prediction_score": score
        })
    return pd.DataFrame(records)


def load_tiger2023(filepath, length_label):
    """TIGER2023: columns = Transcript ID, Target Sequence, Guide Sequence, Guide Score
    
    Target Sequence is the target window on the transcript.
    Guide Sequence is the crRNA (RC of target).
    """
    df = pd.read_csv(filepath)
    log(f"    Loaded TIGER2023: {len(df)} rows")
    
    records = []
    for _, row in df.iterrows():
        target = str(row["Target Sequence"]).strip()
        guide = str(row["Guide Sequence"]).strip()
        score = float(row["Guide Score"])
        records.append({
            "model_name": f"TIGER2023_{length_label}",
            "target_seq": target.upper(),
            "start_pos": None,  # Will be filled after merging with FASTA
            "end_pos": None,
            "prediction_score": score
        })
    result = pd.DataFrame(records)
    result["guide_len"] = result["target_seq"].str.len()
    return result


def load_deepcas13(filepath, length_label):
    """DeepCas13: columns = sgrna, seq, deepscore
    
    sgrna format: sgRNA_X_Y where X=start(0-based), Y=end(0-based exclusive).
    seq is the sgRNA sequence (already RC'd from target).
    """
    df = pd.read_csv(filepath)
    log(f"    Loaded DeepCas13: {len(df)} rows")
    
    records = []
    for _, row in df.iterrows():
        sgrna_id = str(row["sgrna"])
        seq = str(row["seq"]).strip()
        score = float(row["deepscore"])
        
        # Parse sgRNA_0_23 -> start=0, end=22 (0-based)
        parts = sgrna_id.split("_")
        start_0 = int(parts[1])
        end_0 = int(parts[2])
        
        start_pos = start_0 + 1  # Convert to 1-based
        end_pos = end_0  # Already exclusive upper bound = end position
        
        target_seq = reverse_complement(seq)
        records.append({
            "model_name": f"DeepCas13_{length_label}",
            "target_seq": target_seq.upper(),
            "start_pos": start_pos,
            "end_pos": end_pos,
            "prediction_score": score
        })
    return pd.DataFrame(records)


def load_hsu2023(filepath, length_label):
    """Hsu2023 Web: columns = transcript id, guide, target_pos_list, predicted_value_sigmoid
    
    guide = crRNA sequence (RC of target).
    target_pos_list = list of positions (e.g., [181]).
    """
    df = pd.read_csv(filepath)
    log(f"    Loaded Hsu2023: {len(df)} rows")
    
    records = []
    for _, row in df.iterrows():
        guide = str(row["guide"]).strip()
        pos_str = str(row["target_pos_list"]).strip("[]").strip()
        if not pos_str or pos_str == "nan":
            continue
        target_pos = int(pos_str)  # This is the 0-based start position on the transcript
        score = float(row["predicted_value_sigmoid"])
        
        target_seq = reverse_complement(guide)
        guide_len = len(guide)
        
        # Correctly map coordinates based on 0-based target_pos list
        start_pos = target_pos + 1
        end_pos = target_pos + guide_len
        
        records.append({
            "model_name": f"Hsu2023_{length_label}",
            "target_seq": target_seq.upper(),
            "start_pos": start_pos,
            "end_pos": end_pos,
            "prediction_score": score
        })
    return pd.DataFrame(records)


def load_fareh2024(filepath, length_label):
    """Fareh2024: columns = transcript_id, position, guide, Nature2024_Score
    
    position = 1-based start position on the transcript.
    guide = crRNA sequence (RC of target).
    """
    df = pd.read_csv(filepath)
    log(f"    Loaded Fareh2024: {len(df)} rows")
    
    records = []
    for _, row in df.iterrows():
        position = int(row["position"])  # 1-based start
        guide = str(row["guide"]).strip()
        score = float(row["Nature2024_Score"])
        
        target_seq = reverse_complement(guide)
        guide_len = len(guide)
        start_pos = position
        end_pos = position + guide_len - 1
        
        records.append({
            "model_name": f"Fareh2024_{length_label}",
            "target_seq": target_seq.upper(),
            "start_pos": start_pos,
            "end_pos": end_pos,
            "prediction_score": score
        })
    return pd.DataFrame(records)


log("03. Aligning model scores")

# Load FASTA for coordinate mapping (needed for TIGER which has target seq but no coords)
fasta_path = os.path.join(PROJECT, "data/NM_000371.4.fasta")
with open(fasta_path) as f:
    lines = f.readlines()
transcript_seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
log(f"    Transcript length: {len(transcript_seq)} nt")

all_dfs = []
score_files = sorted(os.listdir(SCORE_DIR))

for fname in score_files:
    fpath = os.path.join(SCORE_DIR, fname)
    if not fname.endswith(".csv"):
        continue
    
    log(f"  Processing {fname}...")
    
    if fname.startswith("Sanjana2020_"):
        # Extract length label: Sanjana2020_23_NM_000371.4.csv -> 23
        length_label = fname.split("_")[1]
        df = load_sanjana2020(fpath, length_label)
    
    elif fname.startswith("TIGER2023_"):
        length_label = fname.split("_")[1]
        df = load_tiger2023(fpath, length_label)
        # Map target sequences to positions via exact match in transcript
        targets = []
        starts = []
        ends = []
        for _, row in df.iterrows():
            tseq = row["target_seq"]
            lg = row["guide_len"]
            pos = transcript_seq.find(tseq)
            if pos >= 0:
                targets.append(tseq)
                starts.append(pos + 1)  # 1-based
                ends.append(pos + lg)
            else:
                # Try reverse complement (in case it's already RC'd)
                rc_tseq = reverse_complement(tseq)
                pos = transcript_seq.find(rc_tseq)
                if pos >= 0:
                    targets.append(rc_tseq)
                    starts.append(pos + 1)
                    ends.append(pos + lg)
                else:
                    targets.append(tseq)
                    starts.append(None)
                    ends.append(None)
        df["target_seq"] = targets
        df["start_pos"] = starts
        df["end_pos"] = ends
    
    elif fname.startswith("DeepCas13_"):
        length_label = fname.split("_")[1]
        df = load_deepcas13(fpath, length_label)
    
    elif fname.startswith("Hsu2023_"):
        length_label = fname.split("_")[1]
        df = load_hsu2023(fpath, length_label)
    
    elif fname.startswith("Fareh2024_"):
        length_label = fname.split("_")[1]
        df = load_fareh2024(fpath, length_label)
    
    else:
        log(f"    Unknown format: {fname}, skipping")
        continue
    
    # Filter out rows with None positions
    df = df.dropna(subset=["start_pos", "end_pos"])
    df["start_pos"] = df["start_pos"].astype(int)
    df["end_pos"] = df["end_pos"].astype(int)
    
    # Keep only unified columns
    df = df[["model_name", "target_seq", "start_pos", "end_pos", "prediction_score"]]
    df = df.sort_values("start_pos").reset_index(drop=True)
    
    log(f"    -> {len(df)} aligned rows saved")
    all_dfs.append(df)

# Save each model separately
for df in all_dfs:
    model_name = df["model_name"].iloc[0]
    out_path = os.path.join(RESULT_DIR, f"{model_name}_aligned.csv")
    df.to_csv(out_path, index=False)
    log(f"    Saved {os.path.basename(out_path)} ({len(df)} rows)")

log("03. Done.")