import pandas as pd
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

cntrl_gene = []
types = []
with open("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/targeted_shRNA_screen/control.txt", 'r') as file:
    for line in file:
        line = line.strip()
        cntrl_gene.append(line.split("_")[0])
        types.append(line.split("_")[1])

down_gene = []
for record in SeqIO.parse("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/order_shRNAone_downreg.fa", "fasta"):
    down_gene.append(record.id.split("_")[0])
down_gene = list(set(down_gene))
types.extend(["down_reg"] * len(down_gene))

pos_gene = []
for record in SeqIO.parse("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/shRNA_screen/order_shRNAone_upreg.fa", "fasta"):
    pos_gene.append(record.id.split("_")[0])
pos_gene = list(set(pos_gene))
types.extend(["up_reg"] * len(pos_gene))

genes = cntrl_gene + down_gene + pos_gene
df_scatter = pd.DataFrame({"gene": genes, "type": types})

df_MAGeCK_599 = pd.read_table("/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/599/599.gene_summary.txt")[["id", "neg|lfc"]].rename(columns={"id": "gene"})
df_MAGeCK_600 = pd.read_table("/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/600/600.gene_summary.txt")[["id", "gene_summary.txt" if False else "neg|lfc"]].rename(columns={"id": "gene"})
df_MAGeCK_WT = pd.read_table("/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/WT/WT.gene_summary.txt")[["id", "neg|lfc"]].rename(columns={"id": "gene"})

df_DESeq2 = pd.read_csv("/home/wkkongaa/data/notebooks/250715_novaseq50_reanalysis/DESeq2/miR_494_DESeq2_output.csv")[["Unnamed: 0", "log2FoldChange"]].rename(columns={"Unnamed: 0": "gene"})
df_DESeq2["gene"] = [x.split("|")[0] for x in df_DESeq2["gene"]]
df_scatter = pd.merge(df_scatter, df_DESeq2, how= "left", on="gene")

y_low_lim = min(df_MAGeCK_599["neg|lfc"].to_list() + df_MAGeCK_600["neg|lfc"].to_list() + df_MAGeCK_WT["neg|lfc"].to_list()) - 0.1
y_up_lim = max(df_MAGeCK_599["neg|lfc"].to_list() + df_MAGeCK_600["neg|lfc"].to_list() + df_MAGeCK_WT["neg|lfc"].to_list()) + 0.1

custom_colors = {
    "up_reg": "#e09fa1",
    "down_reg": "#a4adea",
    "PC": "#d3d3d3",
    "anti-ferroptotic": "#b8e987",
    "pro-ferroptotic": "#5b8e7d",
    "anti-apoptotic": "#f3e49c",
    "pro-apoptotic": "#b6a44c",
}

custom_sizes = {
    "up_reg": 15,
    "down_reg": 15,
    "PC": 30,
    "anti-ferroptotic": 30,
    "pro-ferroptotic": 30,
    "anti-apoptotic": 30,
    "pro-apoptotic": 30,
}

desired_order_layered = [
    "up_reg",
    "down_reg",
    "anti-ferroptotic",
    "pro-ferroptotic",
    "anti-apoptotic",
    "pro-apoptotic",
    "PC",
]

dataset_names = {0: "WT", 1: "599", 2: "600"}

for idx, df_raw in enumerate([df_MAGeCK_WT, df_MAGeCK_599, df_MAGeCK_600]):
    df_plot = pd.merge(df_scatter, df_raw, how="left", on="gene")
    df_plot["type"] = pd.Categorical(
        df_plot["type"], categories=desired_order_layered, ordered=True
    )
    df_plot = df_plot.sort_values("type").reset_index(drop=True)

    print(f"--- Dataset: {dataset_names[idx]} ---")
    print(df_plot[df_plot["neg|lfc"] > 0.5])

    x_low_lim = df_plot["log2FoldChange"].min() - 0.3
    x_up_lim = df_plot["log2FoldChange"].max() + 0.3

    plt.figure(figsize=(5, 5.06))
    plt.xlim(x_low_lim, x_up_lim)
    plt.ylim(y_low_lim, y_up_lim)

    sns.scatterplot(
        data=df_plot,
        x="log2FoldChange",
        y="neg|lfc",
        hue="type",
        palette=custom_colors,
        size="type",  
        sizes=custom_sizes,  
        edgecolor="#000000",
        linewidth=0.5,
    )

    # --- ADD LABELS FOR SIGNIFICANT GENES ---
    df_labels = df_plot[(df_plot["neg|lfc"] > 0.5) | (df_plot["neg|lfc"] < -0.5)]
    for _, row in df_labels.iterrows():
        # Ensure values exist to avoid plotting errors
        if pd.notna(row["log2FoldChange"]) and pd.notna(row["neg|lfc"]):
            plt.text(
                x=row["log2FoldChange"] + 0.02, # Small horizontal gap offset
                y=row["neg|lfc"] + 0.01,        # Small vertical gap offset
                s=str(row["gene"]),
                fontsize=8,                     # Clean readable size for PDF
                color="#000000",
                ha="left",                      # Align text to the left side
                va="bottom"                     # Align text above the dot
            )

    plt.title(f"MAGeCK Screening - {dataset_names[idx]}")
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(f"/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/targeted_shRNA_screen/RNAseq_vs_shRNA_screen_{dataset_names[idx]}.pdf")
    plt.show()
