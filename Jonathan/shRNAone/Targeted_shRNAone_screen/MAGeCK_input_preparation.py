df_shRNA_pool = pd.read_csv('/home/wkkongaa/data/notebooks/260416_novaseq99/shRNA_lib/miR494_ref_lib.txt', sep='\t').drop_duplicates()
df_shRNA_pool["Gene"] = shRNA_pool["shRNA_name"].str.split("_").str[0]
df_MAGeCK_input_lib = df_shRNA_pool.rename(columns={"shRNA_name": "sgRNA", "shRNA_seq": "Sequence", "Gene": "Gene"}).to_csv("MAGeCK_input_lib.csv", sep="\t", index=False)

df_WT = df_shRNA_pool.merge(lib_WT_R1_Sora, how='left')
df_WT = df_WT.rename(columns={"shRNA_name": "sgRNA", "Count": "WT_R1_Sora"})

df_WT = df_WT.merge(lib_WT_R1_Veh, how='left')
df_WT = df_WT.rename(columns={"Count": "WT_R1_Veh"})

df_WT = df_WT.merge(lib_WT_R2_Sora, how='left')
df_WT = df_WT.rename(columns={"Count": "WT_R2_Sora"})

df_WT = df_WT.merge(lib_WT_R2_Veh, how='left')
df_WT = df_WT.rename(columns={"Count": "WT_R2_Veh"})

df_WT = df_WT.merge(lib_WT_R3_Sora, how='left')
df_WT = df_WT.rename(columns={"Count": "WT_R3_Sora"})

df_WT = df_WT.merge(lib_WT_R3_Veh, how='left')
df_WT = df_WT.rename(columns={"Count": "WT_R3_Veh"})

df_WT.drop(columns=['shRNA_seq']).to_csv("MAGeCK_input_WT_count.txt",sep = '\t', index=False)

df_599 = df_shRNA_pool.merge(lib_599_R1_Sora, how='left')
df_599 = df_599.rename(columns={"shRNA_name": "sgRNA", "Count": "599_R1_Sora"})

df_599 = df_599.merge(lib_599_R1_Veh, how='left')
df_599 = df_599.rename(columns={"Count": "599_R1_Veh"})

df_599 = df_599.merge(lib_599_R2_Sora, how='left')
df_599 = df_599.rename(columns={"Count": "599_R2_Sora"})

df_599 = df_599.merge(lib_599_R2_Veh, how='left')
df_599 = df_599.rename(columns={"Count": "599_R2_Veh"})

df_599 = df_599.merge(lib_599_R3_Sora, how='left')
df_599 = df_599.rename(columns={"Count": "599_R3_Sora"})

df_599 = df_599.merge(lib_599_R3_Veh, how='left')
df_599 = df_599.rename(columns={"Count": "599_R3_Veh"})

df_599.drop(columns=['shRNA_seq']).to_csv("MAGeCK_input_599_count.txt",sep = '\t', index=False)

df_600 = df_shRNA_pool.merge(lib_600_R1_Sora, how='left')
df_600 = df_600.rename(columns={"shRNA_name": "sgRNA", "Count": "600_R1_Sora"})

df_600 = df_600.merge(lib_600_R1_Veh, how='left')
df_600 = df_600.rename(columns={"Count": "600_R1_Veh"})

df_600 = df_600.merge(lib_600_R2_Sora, how='left')
df_600 = df_600.rename(columns={"Count": "600_R2_Sora"})

df_600 = df_600.merge(lib_600_R2_Veh, how='left')
df_600 = df_600.rename(columns={"Count": "600_R2_Veh"})

df_600 = df_600.merge(lib_600_R3_Sora, how='left')
df_600 = df_600.rename(columns={"Count": "600_R3_Sora"})

df_600 = df_600.merge(lib_600_R3_Veh, how='left')
df_600 = df_600.rename(columns={"Count": "600_R3_Veh"})

df_600.drop(columns=['shRNA_seq']).to_csv("MAGeCK_input_600_count.txt",sep = '\t', index=False)

df_Veh = df_shRNA_pool.merge(lib_599_R1_Veh, how='left')
df_Veh = df_Veh.rename(columns={"shRNA_name": "sgRNA", "Count": "599_R1_Veh"})

df_Veh = df_Veh.merge(lib_599_R2_Veh, how='left')
df_Veh = df_Veh.rename(columns={"Count": "599_R2_Veh"})

df_Veh = df_Veh.merge(lib_599_R3_Veh, how='left')
df_Veh = df_Veh.rename(columns={"Count": "599_R3_Veh"})

df_Veh = df_Veh.merge(lib_600_R1_Veh, how='left')
df_Veh = df_Veh.rename(columns={"Count": "600_R1_Veh"})

df_Veh = df_Veh.merge(lib_600_R2_Veh, how='left')
df_Veh = df_Veh.rename(columns={"Count": "600_R2_Veh"})

df_Veh = df_Veh.merge(lib_600_R3_Veh, how='left')
df_Veh = df_Veh.rename(columns={"Count": "600_R3_Veh"})

df_Veh.drop(columns=['shRNA_seq']).to_csv("MAGeCK_input_count_Veh.txt",sep = '\t', index=False)