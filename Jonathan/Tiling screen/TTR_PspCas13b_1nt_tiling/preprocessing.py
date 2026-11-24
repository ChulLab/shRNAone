PspCas13bDir = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/TTR_PspCas13b_1nt'

Dfiles = [f for f in os.listdir(PspCas13bDir) if f.endswith('sorted.txt')]
Dfiles.sort()


for file in Dfiles:
    file_path = os.path.join(PspCas13bDir, file)
    data = pd.read_csv(file_path, delimiter='\t', names=['Count', 'target_seq'])
    index = file.find("U1Psp-R")
    length = len("U1Psp-R")
    variable_name = ('PspCas13b_rep' + file[index + length:index + length+ 6]).replace("-", "_")
    globals()[variable_name] = data
    print(file, variable_name)


sequence_name = []
sequence = []

for record in SeqIO.parse(Libdirectory, "fasta"):
    if "U1_Psp" not in record.id:
        continue
    else:
        if "rc_" in record.id:
            sequence_name.append(record.id[-9:])
        else:
            sequence_name.append(record.id)
        sequence.append(str(Seq(str(record.seq).upper()[38:-38]).reverse_complement()))
print(sequence_name)

df_ref = pd.DataFrame({"name": sequence_name, "target_seq": sequence})

df_ref["start_nt_NM_000371.3"] = [int(x.split("|")[1].split("-")[0]) if "U1_Psp-crRNA" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.3"] = [int(x.split("|")[1].split("-")[1]) if "U1_Psp-crRNA" in x else "NTC" for x in df_ref["name"]]
df_ref["start_nt_NM_000371.4"] = [int(x.split("|")[1].split("-")[0]) - 110 if "U1_Psp-crRNA" in x else "NTC" for x in df_ref["name"]]
df_ref["end_nt_NM_000371.4"] = [int(x.split("|")[1].split("-")[1]) - 110 if "U1_Psp-crRNA" in x else "NTC" for x in df_ref["name"]]
df_ref["region"] = ["NTC" if x == "NTC" else "CDS" if (x >= 30) & (x <= 470) else "3'UTR" if (x >= 471) & (x <= 616) else "Misannotated" for x in df_ref["end_nt_NM_000371.4"]]

df_PspCas13b = df_ref.merge(PspCas13b_rep1_Bin1, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep1_bin1"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep1_Bin2, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep1_bin2"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep1_Bin3, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep1_bin3"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep1_Bin4, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep1_bin4"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep2_Bin1, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep2_bin1"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep2_Bin2, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep2_bin2"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep2_Bin3, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep2_bin3"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep2_Bin4, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep2_bin4"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep3_Bin1, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep3_bin1"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep3_Bin2, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep3_bin2"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep3_Bin3, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep3_bin3"})

df_PspCas13b = df_PspCas13b.merge(PspCas13b_rep3_Bin4, how='left')
df_PspCas13b = df_PspCas13b.rename(columns={"Count": "rep3_bin4"})

columns_range = df_PspCas13b.columns[7:]
column_sums = df_PspCas13b[columns_range].sum()

for column in range(len(columns_range)):
    df_PspCas13b[columns_range[column] + "_CPM"] = df_PspCas13b[columns_range[column]] / column_sums.iloc[column] * 1000000

df_PspCas13b ["rep1_Log2FC"] = [np.log2((df_PspCas13b["rep1_bin1_CPM"][x]) / (df_PspCas13b["rep1_bin4_CPM"][x])) if (df_PspCas13b["rep1_bin1_CPM"][x] > 100 and df_PspCas13b["rep1_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_PspCas13b))]
df_PspCas13b ["rep2_Log2FC"] = [np.log2((df_PspCas13b["rep2_bin1_CPM"][x]) / (df_PspCas13b["rep2_bin4_CPM"][x])) if (df_PspCas13b["rep2_bin1_CPM"][x] > 100 and df_PspCas13b["rep2_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_PspCas13b))]
df_PspCas13b ["rep3_Log2FC"] = [np.log2((df_PspCas13b["rep3_bin1_CPM"][x]) / (df_PspCas13b["rep3_bin4_CPM"][x])) if (df_PspCas13b["rep3_bin1_CPM"][x] > 100 and df_PspCas13b["rep3_bin4_CPM"][x] > 100) else np.nan for x in range(len(df_PspCas13b))]

numeric_cols = ['rep1_Log2FC', 'rep2_Log2FC', 'rep3_Log2FC']
df_PspCas13b["Log2FC_mean"] = df_PspCas13b[numeric_cols].mean(axis=1)
df_PspCas13b["Log2FC_std"] = df_PspCas13b[numeric_cols].std(axis=1)

df_PspCas13b.to_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_PspCas13b_1nt.csv", index=False)
