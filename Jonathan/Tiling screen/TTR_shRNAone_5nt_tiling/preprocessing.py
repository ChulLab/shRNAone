shRNAoneDir = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/TTR_shRNAone_5nt'
Libdirectory = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/reference/TTR_shRNAone_5nt.fa'



Dfiles = [f for f in os.listdir(shRNAoneDir) if f.endswith('sorted.txt')]
Dfiles.sort()
print(Dfiles)

for file in Dfiles:
    file_path = os.path.join(shRNAoneDir, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'target_seq'])
    variable_name = file[:15].replace("-", "_")
    globals()[variable_name] = data
    print(file, variable_name)


sequence_name = []
sequence = []

for record in SeqIO.parse(Libdirectory, "fasta"):
    sequence_name.append(record.id)
    sequence.append(str(record.seq).upper())

df_ref = pd.DataFrame({"name": sequence_name, "target_seq": sequence})
df_ref

df_ref["tech_rep"] = ["U_starting" if "1U" in x else "A_starting" if "1A" in x else "NTC" for x in df_ref["name"]]
df_ref["start_nt_NM_000371.3"] = [int(x.split("_")[1].split("-")[0]) + 136 if "shTTR" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.3"] = [int(x.split("_")[1].split("-")[1]) + 135 if "shTTR" in x else "NTC" for x in df_ref["name"]]
df_ref["start_nt_NM_000371.4"] = [int(x.split("_")[1].split("-")[0]) + 136 - 110 if "shTTR" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.4"] = [int(x.split("_")[1].split("-")[1]) + 135 - 110 if "shTTR" in x else "NTC" for x in df_ref["name"]]
df_ref["name"] = ["shRNAone-TTR" + x[-3:] + "|" + str(df_ref.loc[df_ref["name"] == x, "start_nt_NM_000371.3"].iloc[0]) + "-" + str(df_ref.loc[df_ref["name"] == x, "end_nt_NM_000371.3"].iloc[0]) if "shTTR" in x else x for x in df_ref["name"]]
df_ref["region"] = ["NTC" if x == "NTC" else "5'UTR" if (x >= 23) & (x <= 26) else "CDS" if (x >= 27) & (x <= 470) else "3'UTR" if (x >= 471) & (x <= 616) else "Misannotated" for x in df_ref["end_nt_NM_000371.4"]]

df_shRNAone_TTR_5nt = df_ref.merge(TTR_br1_KD_bin1, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep1_bin1"}).dropna()

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br1_KD_bin2, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep1_bin2"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br1_KD_bin3, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep1_bin3"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br1_KD_bin4, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep1_bin4"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br2_KD_bin1, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep2_bin1"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br2_KD_bin2, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep2_bin2"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br2_KD_bin3, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep2_bin3"})

df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.merge(TTR_br2_KD_bin4, how='left')
df_shRNAone_TTR_5nt = df_shRNAone_TTR_5nt.rename(columns={"Count": "rep2_bin4"})

df_shRNAone_TTR_5nt['target_seq'] = [x[:-1] for x in df_shRNAone_TTR_5nt['target_seq']]

columns_range = df_shRNAone_TTR_5nt.columns[8:]
column_sums = df_shRNAone_TTR_5nt[columns_range].sum()

for column in range(len(columns_range)):
    df_shRNAone_TTR_5nt[columns_range[column] + "_CPM"] = df_shRNAone_TTR_5nt[columns_range[column]] / column_sums.iloc[column] * 1000000

df_shRNAone_TTR_5nt ["rep1_Log2FC"] = [np.log2((df_shRNAone_TTR_5nt["rep1_bin1_CPM"][x]) / (df_shRNAone_TTR_5nt["rep1_bin4_CPM"][x])) if (df_shRNAone_TTR_5nt["rep1_bin1_CPM"][x] > 100 and df_shRNAone_TTR_5nt["rep1_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_shRNAone_TTR_5nt))]
df_shRNAone_TTR_5nt ["rep2_Log2FC"] = [np.log2((df_shRNAone_TTR_5nt["rep2_bin1_CPM"][x]) / (df_shRNAone_TTR_5nt["rep2_bin4_CPM"][x])) if (df_shRNAone_TTR_5nt["rep2_bin1_CPM"][x] > 100 and df_shRNAone_TTR_5nt["rep2_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_shRNAone_TTR_5nt))]

numeric_cols = ['rep1_Log2FC', 'rep2_Log2FC']
df_shRNAone_TTR_5nt["Log2FC_mean"] = df_shRNAone_TTR_5nt[numeric_cols].mean(axis=1)
df_shRNAone_TTR_5nt["Log2FC_std"] = df_shRNAone_TTR_5nt[numeric_cols].std(axis=1)

df_shRNAone_TTR_5nt.to_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_5nt.csv", index=False)