import os
import glob
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl


Ddirectory = '/home/wkkongaa/data/notebooks/260416_novaseq99/data_read2'
Dfiles = sorted([f for f in os.listdir(Ddirectory) if f.endswith('sorted.txt')])


for file in Dfiles:
    
    file_path = os.path.join(Ddirectory, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'shRNA_seq'])
    variable_name = 'lib_' + file[3:].split('-')[0] + "_" + file[3:].split('-')[2].split('_')[0] + "_" + file[3:].split('-')[1]
    globals()[variable_name] = data
    print(file, variable_name)

shRNA_pool = pd.read_csv('/home/wkkongaa/data/notebooks/260416_novaseq99/shRNA_lib/miR494_ref_lib.txt', sep='\t').drop_duplicates()

df_shRNA_pool = shRNA_pool.merge(lib_WT_R1_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R1_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R1_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R1_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R1_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R1_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R1_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R1_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R2_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R2_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R2_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R2_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R2_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R2_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R2_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R2_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R3_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R3_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R3_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R3_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R3_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R3_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_WT_R3_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "WT_R3_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R1_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R1_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R1_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R1_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R1_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R1_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R1_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R1_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R2_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R2_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R2_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R2_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R2_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R2_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R2_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R2_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R3_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R3_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R3_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R3_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R3_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R3_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_599_R3_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "599_R3_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R1_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R1_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R1_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R1_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R1_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R1_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R1_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R1_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R2_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R2_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R2_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R2_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R2_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R2_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R2_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R2_Veh"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R3_day0, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R3_day0"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R3_day5, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R3_day5"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R3_Sora, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R3_Sora"})

df_shRNA_pool = df_shRNA_pool.merge(lib_600_R3_Veh, how='left')
df_shRNA_pool = df_shRNA_pool.rename(columns={"Count": "600_R3_Veh"})

df_shRNA_pool.to_excel("novaseq99_read_count.xlsx", index=False)

df_shRNA_pool["WT_R1_day0"] = [i * 1000000 / sum(df_shRNA_pool["WT_R1_day0"]) for i in df_shRNA_pool["WT_R1_day0"]]
df_shRNA_pool["WT_R1_day5"] = [i * 1000000 / sum(df_shRNA_pool["WT_R1_day5"]) for i in df_shRNA_pool["WT_R1_day5"]]
df_shRNA_pool["WT_R1_Sora"] = [i * 1000000 / sum(df_shRNA_pool["WT_R1_Sora"]) for i in df_shRNA_pool["WT_R1_Sora"]]
df_shRNA_pool["WT_R1_Veh"] = [i * 1000000 / sum(df_shRNA_pool["WT_R1_Veh"]) for i in df_shRNA_pool["WT_R1_Veh"]]
df_shRNA_pool["WT_R2_day0"] = [i * 1000000 / sum(df_shRNA_pool["WT_R2_day0"]) for i in df_shRNA_pool["WT_R2_day0"]]
df_shRNA_pool["WT_R2_day5"] = [i * 1000000 / sum(df_shRNA_pool["WT_R2_day5"]) for i in df_shRNA_pool["WT_R2_day5"]]
df_shRNA_pool["WT_R2_Sora"] = [i * 1000000 / sum(df_shRNA_pool["WT_R2_Sora"]) for i in df_shRNA_pool["WT_R2_Sora"]]
df_shRNA_pool["WT_R2_Veh"] = [i * 1000000 / sum(df_shRNA_pool["WT_R2_Veh"]) for i in df_shRNA_pool["WT_R2_Veh"]]
df_shRNA_pool["WT_R3_day0"] = [i * 1000000 / sum(df_shRNA_pool["WT_R3_day0"]) for i in df_shRNA_pool["WT_R3_day0"]]
df_shRNA_pool["WT_R3_day5"] = [i * 1000000 / sum(df_shRNA_pool["WT_R3_day5"]) for i in df_shRNA_pool["WT_R3_day5"]]
df_shRNA_pool["WT_R3_Sora"] = [i * 1000000 / sum(df_shRNA_pool["WT_R3_Sora"]) for i in df_shRNA_pool["WT_R3_Sora"]]
df_shRNA_pool["WT_R3_Veh"] = [i * 1000000 / sum(df_shRNA_pool["WT_R3_Veh"]) for i in df_shRNA_pool["WT_R3_Veh"]]
df_shRNA_pool["599_R1_day0"] = [i * 1000000 / sum(df_shRNA_pool["599_R1_day0"]) for i in df_shRNA_pool["599_R1_day0"]]
df_shRNA_pool["599_R1_day5"] = [i * 1000000 / sum(df_shRNA_pool["599_R1_day5"]) for i in df_shRNA_pool["599_R1_day5"]]
df_shRNA_pool["599_R1_Sora"] = [i * 1000000 / sum(df_shRNA_pool["599_R1_Sora"]) for i in df_shRNA_pool["599_R1_Sora"]]
df_shRNA_pool["599_R1_Veh"] = [i * 1000000 / sum(df_shRNA_pool["599_R1_Veh"]) for i in df_shRNA_pool["599_R1_Veh"]]
df_shRNA_pool["599_R2_day0"] = [i * 1000000 / sum(df_shRNA_pool["599_R2_day0"]) for i in df_shRNA_pool["599_R2_day0"]]
df_shRNA_pool["599_R2_day5"] = [i * 1000000 / sum(df_shRNA_pool["599_R2_day5"]) for i in df_shRNA_pool["599_R2_day5"]]
df_shRNA_pool["599_R2_Sora"] = [i * 1000000 / sum(df_shRNA_pool["599_R2_Sora"]) for i in df_shRNA_pool["599_R2_Sora"]]
df_shRNA_pool["599_R2_Veh"] = [i * 1000000 / sum(df_shRNA_pool["599_R2_Veh"]) for i in df_shRNA_pool["599_R2_Veh"]]
df_shRNA_pool["599_R3_day0"] = [i * 1000000 / sum(df_shRNA_pool["599_R3_day0"]) for i in df_shRNA_pool["599_R3_day0"]]
df_shRNA_pool["599_R3_day5"] = [i * 1000000 / sum(df_shRNA_pool["599_R3_day5"]) for i in df_shRNA_pool["599_R3_day5"]]
df_shRNA_pool["599_R3_Sora"] = [i * 1000000 / sum(df_shRNA_pool["599_R3_Sora"]) for i in df_shRNA_pool["599_R3_Sora"]]
df_shRNA_pool["599_R3_Veh"] = [i * 1000000 / sum(df_shRNA_pool["599_R3_Veh"]) for i in df_shRNA_pool["599_R3_Veh"]]
df_shRNA_pool["600_R1_day0"] = [i * 1000000 / sum(df_shRNA_pool["600_R1_day0"]) for i in df_shRNA_pool["600_R1_day0"]]
df_shRNA_pool["600_R1_day5"] = [i * 1000000 / sum(df_shRNA_pool["600_R1_day5"]) for i in df_shRNA_pool["600_R1_day5"]]
df_shRNA_pool["600_R1_Sora"] = [i * 1000000 / sum(df_shRNA_pool["600_R1_Sora"]) for i in df_shRNA_pool["600_R1_Sora"]]
df_shRNA_pool["600_R1_Veh"] = [i * 1000000 / sum(df_shRNA_pool["600_R1_Veh"]) for i in df_shRNA_pool["600_R1_Veh"]]
df_shRNA_pool["600_R2_day0"] = [i * 1000000 / sum(df_shRNA_pool["600_R2_day0"]) for i in df_shRNA_pool["600_R2_day0"]]
df_shRNA_pool["600_R2_day5"] = [i * 1000000 / sum(df_shRNA_pool["600_R2_day5"]) for i in df_shRNA_pool["600_R2_day5"]]
df_shRNA_pool["600_R2_Sora"] = [i * 1000000 / sum(df_shRNA_pool["600_R2_Sora"]) for i in df_shRNA_pool["600_R2_Sora"]]
df_shRNA_pool["600_R2_Veh"] = [i * 1000000 / sum(df_shRNA_pool["600_R2_Veh"]) for i in df_shRNA_pool["600_R2_Veh"]]
df_shRNA_pool["600_R3_day0"] = [i * 1000000 / sum(df_shRNA_pool["600_R3_day0"]) for i in df_shRNA_pool["600_R3_day0"]]
df_shRNA_pool["600_R3_day5"] = [i * 1000000 / sum(df_shRNA_pool["600_R3_day5"]) for i in df_shRNA_pool["600_R3_day5"]]
df_shRNA_pool["600_R3_Sora"] = [i * 1000000 / sum(df_shRNA_pool["600_R3_Sora"]) for i in df_shRNA_pool["600_R3_Sora"]]
df_shRNA_pool["600_R3_Veh"] = [i * 1000000 / sum(df_shRNA_pool["600_R3_Veh"]) for i in df_shRNA_pool["600_R3_Veh"]]

df_shRNA_pool.to_excel("novaseq99_CPM.xlsx", index=False)
