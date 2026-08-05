#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
partA.07.run_SplashRNA.py — SplashRNA target site prediction on NM_000371.4.

This script automates the continuous scoring of SplashRNA for all possible 
22 nt shRNA windows across the transcript, aligning output formats with RNAxs.

conda env: siRNA_splash (has python 2.7, shogun, pandas)
"""

import os
import sys
import bz2
import cPickle
import pandas as pd
import shogun
from datetime import datetime

# ================= 1. Global Configurations =================
# 与 RNAxs 完全一致的目录和命名规则
PROJECT = "/data/cai801/data/HKUcas"
DATA_DIR = os.path.join(PROJECT, "data")
RESULT_DIR = os.path.join(PROJECT, "result/08.manuscript_03/partA/02.original_pred")
ACC = "NM_000371.4"
# SplashRNA 严格使用 22nt 滑动窗口
LEN = 22 
FASTA = os.path.join(DATA_DIR, "%s.fasta" % ACC)
OUT_CSV = os.path.join(RESULT_DIR, "SplashRNA_%d_%s.csv" % (LEN, ACC))

# SplashRNA 特有软件和模型目录
SPLASH_DIR = os.path.join(PROJECT, "software/SplashRNA")
MODEL_30_PATH = os.path.join(SPLASH_DIR, "mir30_libsvm.bz2")
MODEL_E_PATH = os.path.join(SPLASH_DIR, "mirE_libsvm.bz2")

# Python 2.7 兼容的目录创建方式
if not os.path.exists(RESULT_DIR):
    os.makedirs(RESULT_DIR)

def log(msg):
    # Python 2.7 兼容的时间戳打印
    print "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg)
    sys.stdout.flush()

# ================= 2. SplashRNA Feature Core =================
class ShRNAtor(object):
    def __init__(self):
        self.svm = None
    def predict(self, test_examples):
        test_feat = construct_features(test_examples)
        pred = self.svm.apply(test_feat)
        return pred.get_values()

def construct_features(features):
    feat_all = [inst for inst in features]
    feat_lhs = [inst[0:15] for inst in features]
    feat_rhs = [inst[15:] for inst in features]

    feat_wd = shogun.StringCharFeatures(shogun.DNA)
    feat_wd.set_features(feat_all)
    
    feat_spec_1 = get_spectrum_features(feat_lhs, 3)
    feat_spec_2 = get_spectrum_features(feat_rhs, 3)

    feat_comb = shogun.CombinedFeatures()
    feat_comb.append_feature_obj(feat_wd)
    feat_comb.append_feature_obj(feat_spec_1)
    feat_comb.append_feature_obj(feat_spec_2)
    return feat_comb

def get_spectrum_features(data, order):
    charfeat = shogun.StringCharFeatures(data, shogun.DNA)
    feat = shogun.StringWordFeatures(charfeat.get_alphabet())
    feat.obtain_from_char(charfeat, order-1, order, 0, True)
    preproc = shogun.SortWordString()
    preproc.init(feat)
    feat.add_preprocessor(preproc)
    feat.apply_preprocessor()
    return feat

# ================= 3. Main Execution =================
def main():
    log("07. SplashRNA: %s" % ACC)
    log("    Guide length: %d nt (native shRNA window length)" % LEN)

    # 1. Read FASTA (与 RNAxs 逻辑对齐)
    with open(FASTA, "r") as f:
        lines = f.readlines()
    seq_flat = "".join([line.strip().upper().replace('U', 'T') for line in lines if not line.startswith(">")])
    log("    Transcript length: %d nt" % len(seq_flat))

    # 2. Build sliding windows and records
    records = []
    candidates = []
    for i in range(len(seq_flat) - LEN + 1):
        target_sub = seq_flat[i:i + LEN]
        if "N" in target_sub:
            continue
        candidates.append(target_sub)
        records.append({
            'transcript_id': ACC,
            'position': i + 1, # 1-based start (5' end of target)
            'guide': target_sub
        })

    log("    Generated %d candidate %d-nt sequences." % (len(candidates), LEN))

    # 3. Execute Prediction
    log("    Loading LibSVM models from %s..." % SPLASH_DIR)
    try:
        with bz2.BZ2File(MODEL_30_PATH, 'rb') as f:
            svm_mir30 = cPickle.load(f)
        with bz2.BZ2File(MODEL_E_PATH, 'rb') as f:
            svm_mirE = cPickle.load(f)
    except Exception as e:
        log("    [ERROR] Model loading failed: %s" % str(e))
        sys.exit(1)

    log("    Predicting SplashRNA cascade scores...")
    scores_30 = svm_mir30.predict(candidates)
    scores_E = svm_mirE.predict(candidates)

    # Calculate final cascade scores
    for idx in range(len(records)):
        s30 = scores_30[idx]
        sE = scores_E[idx]
        
        # 1.1 threshold gatekeeper logic
        if s30 >= 1.1:
            final = 0.6 * s30 + 0.4 * sE
        else:
            final = s30
        
        records[idx]['miR30_Score'] = s30
        records[idx]['miR-E_Score'] = sE
        records[idx]['Splash_Score'] = final

    # ================= 4. Merge & Export =================
    df = pd.DataFrame(records)
    
    if not df.empty:
        # Sort primarily by final Splash_Score descending
        df = df.sort_values(by=['Splash_Score'], ascending=False)
        
        # Rank generation based on final score
        df['rank'] = df['Splash_Score'].rank(method='min', ascending=False)
        
        # Reorder columns for neatness
        cols = ['transcript_id', 'position', 'guide', 'miR30_Score', 'miR-E_Score', 'Splash_Score', 'rank']
        df = df[cols]
        
        df.to_csv(OUT_CSV, index=False)
        log("    %d predictions merged and evaluated -> %s" % (len(df), os.path.basename(OUT_CSV)))
    else:
        log("    WARNING: No predictions were generated.")

    log("07. Done.")

if __name__ == "__main__":
    main()