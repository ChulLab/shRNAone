import scipy.stats as stats
import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


shRNAoneDir = '/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/miR_shRNAone'
Libdirectory = '/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/miR_screening_ref.txt'



Dfiles = [f for f in os.listdir(shRNAoneDir) if f.endswith('sorted.txt')]
Dfiles.sort()
print(Dfiles)

for file in Dfiles:
    file_path = os.path.join(shRNAoneDir, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'complementary_seq'])
    if "SF015" in file:
        variable_name = 'rep1_day' + file[7:].split("_L1_2")[0].replace("-", "_")
    else:
        variable_name = 'rep' + file[8:9] + "_day" + file[11:].split("_L")[0]
    globals()[variable_name] = data
    print(file, variable_name)


sequence_name = []
sequence = []

df_ref = pd.read_table("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/miR_screening_ref.txt", sep="\t")
df_ref['guide_seq'] = [str(Seq(x).reverse_complement()) for x in df_ref['complementary_seq']]
df_ref["type"] = [x.split("|")[0] for x in df_ref["name"]]
df_ref["miRNA"] = [x.split("|")[1] for x in df_ref["name"]]


df_shRNAone = df_ref.merge(rep1_day00, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_day00"})

df_shRNAone = df_shRNAone.merge(rep1_day05, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_day05"})

df_shRNAone = df_shRNAone.merge(rep1_day12_DMSO, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_day12_DMSO"})

df_shRNAone = df_shRNAone.merge(rep1_day12_Sorafenib, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep1_day12_Sorafenib"})

df_shRNAone = df_shRNAone.merge(rep2_day00, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_day00"})

df_shRNAone = df_shRNAone.merge(rep2_day05, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_day05"})

df_shRNAone = df_shRNAone.merge(rep2_day12_DMSO, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_day12_DMSO"})

df_shRNAone = df_shRNAone.merge(rep2_day12_Sorafenib, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep2_day12_Sorafenib"})

df_shRNAone = df_shRNAone.merge(rep3_day00, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_day00"})

df_shRNAone = df_shRNAone.merge(rep3_day05, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_day05"})

df_shRNAone = df_shRNAone.merge(rep3_day12_DMSO, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_day12_DMSO"})

df_shRNAone = df_shRNAone.merge(rep3_day12_Sorafenib, how='left')
df_shRNAone = df_shRNAone.rename(columns={"Count": "rep3_day12_Sorafenib"})

df_shRNAone = df_shRNAone.drop(columns=['complementary_seq'])

columns_range = df_shRNAone.columns[4:]
column_sums = df_shRNAone[columns_range].sum()

for column in range(len(columns_range)):
    df_shRNAone[columns_range[column] + "_CPM"] = df_shRNAone[columns_range[column]] / column_sums.iloc[column] * 1000000

df_shRNAone["rep1_Log2FC"] = [np.log2((df_shRNAone["rep1_day12_DMSO"][x]) / (df_shRNAone["rep1_day12_Sorafenib"][x])) if (df_shRNAone["rep1_day12_DMSO"][x] > 50 and df_shRNAone["rep1_day12_Sorafenib"][x] > 50) else np.nan for x in range(len(df_shRNAone))]
df_shRNAone["rep2_Log2FC"] = [np.log2((df_shRNAone["rep2_day12_DMSO"][x]) / (df_shRNAone["rep2_day12_Sorafenib"][x])) if (df_shRNAone["rep2_day12_DMSO"][x] > 50 and df_shRNAone["rep2_day12_Sorafenib"][x] > 50) else np.nan for x in range(len(df_shRNAone))]
df_shRNAone["rep3_Log2FC"] = [np.log2((df_shRNAone["rep3_day12_DMSO"][x]) / (df_shRNAone["rep3_day12_Sorafenib"][x])) if (df_shRNAone["rep3_day12_DMSO"][x] > 50 and df_shRNAone["rep3_day12_Sorafenib"][x] > 50) else np.nan for x in range(len(df_shRNAone))]

numeric_cols = ['rep1_Log2FC', 'rep2_Log2FC', 'rep3_Log2FC']
df_shRNAone["Log2FC_mean"] = df_shRNAone[numeric_cols].mean(axis=1)
df_shRNAone["Log2FC_std"] = df_shRNAone[numeric_cols].std(axis=1)

df_shRNAone.to_csv("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/Log2FC_miR_shRNAone_biological_replicates.csv", index=False)
matrix_DESeq2 = df_shRNAone[["name", "rep1_day00", "rep1_day05", "rep1_day12_DMSO", "rep1_day12_Sorafenib", "rep2_day00", "rep2_day05", "rep2_day12_DMSO", "rep2_day12_Sorafenib", "rep3_day00", "rep3_day05", "rep3_day12_DMSO", "rep3_day12_Sorafenib"]]
matrix_DESeq2.to_csv("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/DESeq2_input/shRNAone_miR_matrix.csv", index = False)