df_shRNAone = pd.read_excel("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/TRC-like shRNAone.xlsx")
df_shRNAone['NCBI_geneId'] = df_shRNAone['NCBI_geneId'].astype('str')

df_shRNAone_upreg = df_shRNAone[df_shRNAone["NCBI_geneId"].isin([ID + ".0" for ID in pd.Series(NCBI_geneIds_up, dtype='str')])]
df_shRNAone_downreg = df_shRNAone[df_shRNAone["NCBI_geneId"].isin([ID + ".0" for ID in pd.Series(NCBI_geneIds_down, dtype='str')])]

df_shRNAone_upreg.sort_values(by='NCBI_geneId', ascending=True, inplace=True)
df_shRNAone_upreg = df_shRNAone_upreg.sort_values(by=['NCBI_geneId', 'cloneId'], ascending=[True, True]).groupby('NCBI_geneId').head(5)
df_shRNAone_upreg.to_excel("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/shRNAone_lib_upreg_gene.xlsx", index=False)
shRNAone_upreg_seq = df_shRNAone_upreg['shRNAone'].tolist()
shRNAone_upreg_ID = [ID.split("|")[0] for ID in df_shRNAone_upreg['cloneId'].tolist()]
shRNAone_upreg_name = [df_shRNAone_upreg['symbol'].tolist()[ID] + "_" + shRNAone_upreg_ID[ID] for ID in range(len(shRNAone_upreg_ID))]
forward_PBS = ["ATTTACTATCATCCCAAG", "GTTTACATTCACCAGGAT", "TGATCCTCCAAAGAAAGC"]
reverse_PBS = ["CGTTGTTGATAAGGTTTA", "CATCGTCAATTTAAGGGA", "TATAGAAAAACCAAAACC"]
with open("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/order_shRNAone_upreg.fa", "w") as f:
    for idx in range(len(shRNAone_upreg_seq)):
        if idx < 1000:
            PBS = 0
        elif idx < 2000 and idx >= 1000:
            PBS = 1
        else:
            PBS = 2
        f.write(f'>{shRNAone_upreg_name[idx]}\n{forward_PBS[PBS]}ctcgag{shRNAone_upreg_seq[idx]}gaattc{reverse_PBS[PBS]}\n')

df_shRNAone_downreg.sort_values(by='NCBI_geneId', ascending=True, inplace=True)
df_shRNAone_downreg = df_shRNAone_downreg.sort_values(by=['NCBI_geneId', 'cloneId'], ascending=[True, True]).groupby('NCBI_geneId').head(5)
df_shRNAone_downreg.to_excel("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/shRNAone_lib_downreg_gene.xlsx", index=False)
shRNAone_downreg_seq = df_shRNAone_downreg['shRNAone'].tolist()
shRNAone_downreg_ID = [ID.split("|")[0] for ID in df_shRNAone_downreg['cloneId'].tolist()]
shRNAone_downreg_name = [df_shRNAone_downreg['symbol'].astype(str).tolist()[ID] + "_" + shRNAone_downreg_ID[ID] for ID in range(len(shRNAone_downreg_ID))]
forward_PBS = "TAGAAAACGATAGAAGTA"
reverse_PBS = "CAATAACAATCGCTCCAG"
with open("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/order_shRNAone_downreg.fa", "w") as f:
    for idx in range(len(shRNAone_downreg_seq)):
        f.write(f'>{shRNAone_downreg_name[idx]}\n{forward_PBS}ctcgag{shRNAone_downreg_seq[idx]}gaattc{reverse_PBS}\n')

with open("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/filtered_controls.txt", "r") as cntrl:
    control_genes = [line.strip() for line in cntrl]
df_shRNAone_controls = df_shRNAone[df_shRNAone["symbol"].isin(control_genes)]
df_shRNAone_controls = df_shRNAone_controls.sort_values(by=['NCBI_geneId', 'cloneId'], ascending=[True, True]).groupby('NCBI_geneId').tail(5)
df_shRNAone_controls.to_excel("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/shRNAone_lib_controls.xlsx", index=False)
shRNAone_controls_seq = df_shRNAone_controls['shRNAone'].tolist()
shRNAone_controls_ID = [ID.split("|")[0] for ID in df_shRNAone_controls['cloneId'].tolist()]
shRNAone_controls_name = [df_shRNAone_controls['symbol'].astype(str).tolist()[ID] + "_" + shRNAone_controls_ID[ID] for ID in range(len(shRNAone_controls_ID))]
forward_PBS = "ACCTACAAGACACCCATG"
reverse_PBS = "CTCAGTCCTTTCGCATAT"
with open("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/order_shRNAone_controls.fa", "w") as f:
    for idx in range(len(shRNAone_controls_seq)):
        f.write(f'>{shRNAone_controls_name[idx]}\n{forward_PBS}ctcgag{shRNAone_controls_seq[idx]}gaattc{reverse_PBS}\n')