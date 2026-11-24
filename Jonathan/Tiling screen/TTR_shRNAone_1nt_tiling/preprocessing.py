import scipy.stats as stats
import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


shRNAoneDir = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/TTR_shRNAone_1nt'
Libdirectory = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/reference/hkuoligopool9.fa'



Dfiles = [f for f in os.listdir(shRNAoneDir) if f.endswith('sorted.txt')]
Dfiles.sort()
print(Dfiles)

for file in Dfiles:
    file_path = os.path.join(shRNAoneDir, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'target_seq'])
    variable_name = 'shRNAone_rep' + file[7:13]
    if variable_name[-1] == 'u':
        variable_name += file[13:14]
    globals()[variable_name] = data
    print(file, variable_name)


sequence_name = []
sequence = []

for record in SeqIO.parse(Libdirectory, "fasta"):
    if "shRNAone" not in record.id:
        continue
    else:
        sequence_name.append(record.id)
        sequence.append(str(record.seq).upper())

sequence = [str(Seq(i[51:70]).reverse_complement()) for i in sequence]

df_ntc = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/reference/hkuoligopool2.txt", sep="\t")

df_ntc = df_ntc[df_ntc["Library"] == "control.1U"]
list_ntc_name = [l.split("|")[1] for l in df_ntc["OligoName"]]
list_ntc_seq = [str(Seq(l[51:70].upper()).reverse_complement()) for l in df_ntc["OligoSeq"]]
sequence_name = sequence_name + list_ntc_name
sequence = sequence + list_ntc_seq

df_ref = pd.DataFrame({"name": sequence_name, "target_seq": sequence})
df_ref

df_ref["start_nt_NM_000371.3"] = [int(x.split("|")[1].split("-")[0]) if "shRNAone-TTR" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.3"] = [int(x.split("|")[1].split("-")[1]) if "shRNAone-TTR" in x else "NTC" for x in df_ref["name"]]
df_ref["start_nt_NM_000371.4"] = [int(x.split("|")[1].split("-")[0]) - 110 if "shRNAone-TTR" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.4"] = [int(x.split("|")[1].split("-")[1]) - 110 if "shRNAone-TTR" in x else "NTC" for x in df_ref["name"]]
df_ref["region"] = ["NTC" if x == "NTC" else "5'UTR" if (x >= 18) & (x <= 26) else "CDS" if (x >= 27) & (x <= 470) else "3'UTR" if (x >= 471) & (x <= 616) else "Misannotated" for x in df_ref["end_nt_NM_000371.4"]]

df_shRNAone = df_ref.merge(shRNAone_rep1_Bin1, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_bin1"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep1_Bin2, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_bin2"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep1_Bin3, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_bin3"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep1_Bin4, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_bin4"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep1_Input, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_input"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep2_Bin1, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_bin1"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep2_Bin2, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_bin2"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep2_Bin3, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_bin3"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep2_Bin4, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_bin4"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep2_Input, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_input"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep3_Bin1, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_bin1"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep3_Bin2, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_bin2"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep3_Bin3, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_bin3"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep3_Bin4, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_bin4"})

df_shRNAone = df_shRNAone.merge(shRNAone_rep3_Input, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_input"})

df_shRNAone['target_seq'] = [x[:-1] for x in df_shRNAone['target_seq']]

columns_range = df_shRNAone.columns[7:]
column_sums = df_shRNAone[columns_range].sum()

for column in range(len(columns_range)):
    df_shRNAone[columns_range[column] + "_CPM"] = df_shRNAone[columns_range[column]] / column_sums.iloc[column] * 1000000

df_shRNAone["rep1_Log2FC"] = [np.log2((df_shRNAone["rep1_bin1_CPM"][x]) / (df_shRNAone["rep1_bin4_CPM"][x])) if (df_shRNAone["rep1_bin1_CPM"][x] > 100 and df_shRNAone["rep1_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_shRNAone))]
df_shRNAone["rep2_Log2FC"] = [np.log2((df_shRNAone["rep2_bin1_CPM"][x]) / (df_shRNAone["rep2_bin4_CPM"][x])) if (df_shRNAone["rep2_bin1_CPM"][x] > 100 and df_shRNAone["rep2_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_shRNAone))]
df_shRNAone["rep3_Log2FC"] = [np.log2((df_shRNAone["rep3_bin1_CPM"][x]) / (df_shRNAone["rep3_bin4_CPM"][x])) if (df_shRNAone["rep3_bin1_CPM"][x] > 100 and df_shRNAone["rep3_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_shRNAone))]

numeric_cols = ['rep1_Log2FC', 'rep2_Log2FC']
df_shRNAone["Log2FC_mean"] = df_shRNAone[numeric_cols].mean(axis=1)
df_shRNAone["Log2FC_std"] = df_shRNAone[numeric_cols].std(axis=1)

df_shRNAone.to_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_1nt.csv", index=False)