#!/usr/bin/env python3
# partC.10.siRNA_clean_acc.py
#

import os
import pandas as pd
import subprocess
from multiprocessing import Pool
from datetime import datetime

# ================= 0. 目录与初始化 =================
BASE_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partC"
GT_DIR = os.path.join(BASE_DIR, "02.GT_siRNA_info")
CAND_FILE = os.path.join(BASE_DIR, "01.Self_WL_list", "03.selected_WL_Candidates.csv")

OUT_DIR = os.path.join(BASE_DIR, "03.GT_siRNA_acc", "00.clean_siRNA")
LUNP_DIR = os.path.join(BASE_DIR, "03.GT_siRNA_acc", "01.lunp_files")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LUNP_DIR, exist_ok=True)

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

# ================= 1. 清洗 siRNA 数据 =================
log("Step 1: Cleaning siRNA info...")

# 读取原始 siRNA 注释 (仅保留 PMID 16025102 的文献数据)
df = pd.read_csv(os.path.join(GT_DIR, "siRNA_all.txt"), sep="\t")
df = df[df["PMID"] == 16025102].copy()

# 读取并格式化 FASTA 文件
fasta_file = os.path.join(GT_DIR, "sequence-siRNA.fasta")
gs = {}
cid = ""
cs = []
with open(fasta_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if cid: 
                gs[cid] = "".join(cs).upper().replace("U", "T")
            cid = line[1:].split()[0].split(".")[0]
            cs = []
        else:
            cs.append(line)
    if cid: 
        gs[cid] = "".join(cs).upper().replace("U", "T")

# 互补反向序列函数
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

rows = []
unmapped = 0
for _, r in df.iterrows():
    acc = r["Accession_number"]
    anti = str(r["Antisense_21mer"]).upper().replace("U", "T")
    gseq = gs.get(acc, "")
    
    if not gseq:
        unmapped += 1
        continue
        
    # 匹配 target 位置
    idx = gseq.find(anti)
    if idx == -1: 
        idx = gseq.find(rc(anti))
    if idx == -1:
        unmapped += 1
        continue
        
    rows.append({
        "ID": r["ID"],
        "Gene": r["Gene"],
        "Accession": acc,
        "i_score": r["i-score"],
        "target_start": idx + 1,
        "target_end": idx + 21
    })

all_siRNA = pd.DataFrame(rows)
log(f"Mapped {len(all_siRNA)} siRNAs. {unmapped} unmapped.")

# 输出 Continuous 版本
all_siRNA.to_csv(os.path.join(OUT_DIR, "siRNA_continuous.csv"), index=False)

# 输出 Binary 版本
func_thresh = 70
nonfunc_thresh = 30
all_siRNA["group"] = "intermediate"
all_siRNA.loc[all_siRNA["i_score"] > func_thresh, "group"] = "functional"
all_siRNA.loc[all_siRNA["i_score"] < nonfunc_thresh, "group"] = "nonfunctional"

binary_siRNA = all_siRNA[all_siRNA["group"] != "intermediate"]
binary_siRNA.to_csv(os.path.join(OUT_DIR, "siRNA_binary.csv"), index=False)
log(f"Saved continuous ({len(all_siRNA)}) and binary ({len(binary_siRNA)}) versions to {OUT_DIR}.")


# ================= 2. 多进程 RNAplfold 运算 =================
log("Step 2: Preparing RNAplfold tasks...")

# 提取唯一的 W 和 L 组合
cands = pd.read_csv(CAND_FILE)
wl_pairs = cands[["W", "L"]].drop_duplicates()
genes = list(all_siRNA["Accession"].unique())

tasks = []
for _, row in wl_pairs.iterrows():
    w = int(row["W"])
    l = int(row["L"])
    for acc in genes:
        tasks.append((acc, gs[acc], w, l))

log(f"Total RNAplfold tasks to compute: {len(tasks)}")

def run_plfold(task):
    acc, seq, w, l = task
    # 利用序列头命名避免同目录下的输出文件相互覆盖
    job_id = f"{acc}_W{w}_L{l}"
    fasta_path = os.path.join(LUNP_DIR, f"{job_id}.fa")
    out_lunp = os.path.join(LUNP_DIR, f"{job_id}_lunp")

    # 断点续传：若已存在对应的 lunp 文件，则跳过
    if os.path.exists(out_lunp):
        return True

    # 写入用于单次计算的 fasta
    with open(fasta_path, "w") as f:
        f.write(f">{job_id}\n{seq}\n")

    # 执行计算 (统一将 u 设定为 30，完全覆盖下游所需的范围)
    subprocess.run(
        ["RNAplfold", "-W", str(w), "-L", str(l), "-u", "30"],
        cwd=LUNP_DIR,
        stdin=open(fasta_path),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    # 清理运算残余临时文件
    ps_file = os.path.join(LUNP_DIR, f"{job_id}_dp.ps")
    if os.path.exists(fasta_path): 
        os.remove(fasta_path)
    if os.path.exists(ps_file): 
        os.remove(ps_file)
        
    return True

log("Starting multiprocessing pool with 20 cores...")
with Pool(20) as p:
    for i, _ in enumerate(p.imap_unordered(run_plfold, tasks)):
        if (i + 1) % 100 == 0:
            log(f"  Progress: {i+1}/{len(tasks)} completed")

log("Module partC.10 finished successfully. All _lunp files are ready for matrix extraction.")