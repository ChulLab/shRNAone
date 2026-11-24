import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import itertools

BASE_DIR = "/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2"

df = pd.read_csv(f"{BASE_DIR}/miR_494_DESeq2_output.csv").dropna()
df.columns = ['gene', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj']
df = df[(df["baseMean"] > 50) & (df["pvalue"] < 0.05)].sort_values(
    by='log2FC', ascending=True
).reset_index(drop=True)

with open(f"{BASE_DIR}/shRNA_screen/filtered_downreg_jonathan.txt", "w") as f:
    for item in df[df["log2FC"] < -0.5]["gene"]:
        f.write(f'{item.split("|")[0]}\n')

with open(f"{BASE_DIR}/shRNA_screen/filtered_upreg_jonathan.txt", "w") as f:
    for item in df[df["log2FC"] > 0.5]["gene"]:
        f.write(f'{item.split("|")[0]}\n')

UTR_DIR = "/home/wkkongaa/data/transcriptome_datasets/MANE_Nature_2022/3utr_identification"
input_file = f"{UTR_DIR}/MANE_GRCh38_v1.4.3utr.fa"
output_file = f"{BASE_DIR}/shRNA_screen/gene_with_seed_target_MANE.fa"

# Define the 7nt seed (miRNA positions 2-8, RNA with U) — example: miR-494
seed = "GAAACAU".upper().replace("U", "T")
site_8mer    = str(Seq(seed).reverse_complement()) + "A"
site_7mer_m8 = str(Seq(seed).reverse_complement())
site_7mer_A1 = str(Seq(seed[0:6]).reverse_complement()) + "A"
site_6mer    = str(Seq(seed[0:6]).reverse_complement())
all_complementary = [site_8mer, site_7mer_m8, site_7mer_A1, site_6mer]
print(all_complementary)

target_genes_with_seed_target = []

with open(output_file, "w") as out_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq).upper()[15:]
        include = False
        for motif in all_complementary:
            if include:
                break  # Early exit if already found for one motif
            if seq.find(motif) != -1:
                target_genes_with_seed_target.append(record.id.split("|")[0])
                include = True
                break
        if include:
            SeqIO.write(record, out_handle, "fasta")

upreg_genes, downreg_genes = [], []

with open(f"{BASE_DIR}/shRNA_screen/filtered_downreg_jonathan.txt", 'r') as file:
    downreg_genes = [line.strip() for line in file]

with open(f"{BASE_DIR}/shRNA_screen/filtered_upreg_jonathan.txt", 'r') as file:
    upreg_genes = [line.strip() for line in file]

upreg_genes_target = set(upreg_genes) & set(target_genes_with_seed_target)
downreg_genes_target = set(downreg_genes) & set(target_genes_with_seed_target)

print(len(upreg_genes), len(upreg_genes_target), len(downreg_genes), len(downreg_genes_target))

with open(f"{BASE_DIR}/shRNA_screen/filtered_upreg_target_jonathan.txt", "w") as f:
    for item in upreg_genes_target:
        f.write(f'{item}\n')

with open(f"{BASE_DIR}/shRNA_screen/filtered_downreg_target_jonathan.txt", "w") as f:
    for item in downreg_genes_target:
        f.write(f'{item}\n')

file_path = [
    f"{BASE_DIR}/shRNA_screen/filtered_upreg_target_jonathan.txt",
    f"{BASE_DIR}/shRNA_screen/filtered_downreg_target_jonathan.txt",
    f"{BASE_DIR}/shRNA_screen/filtered_upreg_chul.txt",
    f"{BASE_DIR}/shRNA_screen/filtered_downreg_chul.txt",
]

list_deg = []
for path in file_path:
    with open(path, 'r') as f:
        list_deg.append(f.read().splitlines())

up_reg_genes = set(list_deg[0]) | set(list_deg[2])
down_reg_genes = set(list_deg[1]) | set(list_deg[3])

df_TRC = pd.read_csv(f"{UTR_DIR}/MANE.GRCh38.v1.5.summary.txt", sep="\t")

mane_genes_upreg = set(up_reg_genes) & set(df_TRC["symbol"].tolist())
with open(f"{BASE_DIR}/shRNA_screen/filtered_upreg_with_NCBIgeneID.txt", 'w') as f:
    for item in mane_genes_upreg:
        f.write(f'{item}\n')

mane_genes_downreg = set(down_reg_genes) & set(df_TRC["symbol"].tolist())
with open(f"{BASE_DIR}/shRNA_screen/filtered_downreg_with_NCBIgeneID.txt", 'w') as f:
    for item in mane_genes_downreg:
        f.write(f'{item}\n')

df_MANE_down = df_TRC[df_TRC["symbol"].isin(mane_genes_downreg)].reset_index(drop=True)
NCBI_geneIds_down = [ID.split("GeneID:")[1] for ID in set(df_MANE_down["NCBI_GeneID"].tolist())]

df_MANE_up = df_TRC[df_TRC["symbol"].isin(mane_genes_upreg)].reset_index(drop=True)
NCBI_geneIds_up = [ID.split("GeneID:")[1] for ID in set(df_MANE_up["NCBI_GeneID"].tolist())]