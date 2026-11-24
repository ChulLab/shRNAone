import scipy.stats as stats
import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

shRNAoneDir = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/PCSK9_shRNAone_3nt'
Libdirectory = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/reference/PCSK9_shRNAone_3nt.fa'



Dfiles = [f for f in os.listdir(shRNAoneDir) if f.endswith('sorted.txt')]
Dfiles.sort()
print(Dfiles)

for file in Dfiles:
    file_path = os.path.join(shRNAoneDir, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'target_seq'])
    variable_name = file[:18].replace("-", "_")
    globals()[variable_name] = data
    print(file, variable_name)


sequence_name = []
sequence = []

for record in SeqIO.parse(Libdirectory, "fasta"):
    sequence_name.append(record.id)
    sequence.append(str(record.seq).upper())

df_ref = pd.DataFrame({"name": sequence_name, "target_seq": sequence})
df_ref


df_ref["start_nt_NM_174936.3"] = [int(x.split("|")[1].split("-")[0]) if "shRNAone" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_174936.3"] = [int(x.split("|")[1].split("-")[1]) if "shRNAone" in x else "NTC" for x in df_ref["name"]]
df_ref["region"] = ["NTC" if x == "NTC" else "5'UTR" if x <= 362 else "CDS" if (x > 362) & (x <= 2441) else "3'UTR" for x in df_ref["end_nt_NM_174936.3"]]


df_shRNAone_PCSK9_3nt = df_ref.merge(PCSK9_rep1_KD_bin1, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep1_bin1"}).dropna()

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep1_KD_bin2, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep1_bin2"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep1_KD_bin3, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep1_bin3"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep1_KD_bin4, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep1_bin4"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep2_KD_bin1, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep2_bin1"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep2_KD_bin2, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep2_bin2"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep2_KD_bin3, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep2_bin3"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep2_KD_bin4, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep2_bin4"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep3_KD_bin1, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep3_bin1"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep3_KD_bin2, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep3_bin2"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep3_KD_bin3, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep3_bin3"})

df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.merge(PCSK9_rep3_KD_bin4, how='left')
df_shRNAone_PCSK9_3nt = df_shRNAone_PCSK9_3nt.rename(columns={"Count": "rep3_bin4"})

df_shRNAone_PCSK9_3nt['target_seq'] = [x[:-1] for x in df_shRNAone_PCSK9_3nt['target_seq']]

columns_range = df_shRNAone_PCSK9_3nt.columns[5:]
column_sums = df_shRNAone_PCSK9_3nt[columns_range].sum()

for column in range(len(columns_range)):
    df_shRNAone_PCSK9_3nt[columns_range[column] + "_CPM"] = df_shRNAone_PCSK9_3nt[columns_range[column]] / column_sums.iloc[column] * 1000000


df_shRNAone_PCSK9_3nt ["rep1_Log2FC"] = [np.log2((df_shRNAone_PCSK9_3nt["rep1_bin1_CPM"][x]) / (df_shRNAone_PCSK9_3nt["rep1_bin4_CPM"][x])) if (df_shRNAone_PCSK9_3nt["rep1_bin1_CPM"][x] > 40 and df_shRNAone_PCSK9_3nt["rep1_bin4_CPM"][x] > 40) else np.nan for x in range(len(df_shRNAone_PCSK9_3nt))]
df_shRNAone_PCSK9_3nt ["rep2_Log2FC"] = [np.log2((df_shRNAone_PCSK9_3nt["rep2_bin1_CPM"][x]) / (df_shRNAone_PCSK9_3nt["rep2_bin4_CPM"][x])) if (df_shRNAone_PCSK9_3nt["rep2_bin1_CPM"][x] > 40 and df_shRNAone_PCSK9_3nt["rep2_bin4_CPM"][x] > 40) else np.nan for x in range(len(df_shRNAone_PCSK9_3nt))]
df_shRNAone_PCSK9_3nt ["rep3_Log2FC"] = [np.log2((df_shRNAone_PCSK9_3nt["rep3_bin1_CPM"][x]) / (df_shRNAone_PCSK9_3nt["rep3_bin4_CPM"][x])) if (df_shRNAone_PCSK9_3nt["rep3_bin1_CPM"][x] > 40 and df_shRNAone_PCSK9_3nt["rep3_bin4_CPM"][x] > 40) else np.nan for x in range(len(df_shRNAone_PCSK9_3nt))]

numeric_cols = ['rep2_Log2FC', 'rep3_Log2FC']
df_shRNAone_PCSK9_3nt["Log2FC_mean"] = df_shRNAone_PCSK9_3nt[numeric_cols].mean(axis=1)
df_shRNAone_PCSK9_3nt["Log2FC_std"] = df_shRNAone_PCSK9_3nt[numeric_cols].std(axis=1)

df_shRNAone_PCSK9_3nt.to_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_PCSK9_shRNAone_3nt.csv", index=False)