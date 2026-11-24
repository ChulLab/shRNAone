import math
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

df_DESeq2 = pd.read_excel("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/Log2FC_miR_shRNAone_DESeq2.xlsx")

df_DESeq2_miR = df_DESeq2[~(df_DESeq2["type"] == "NTC") & ~(df_DESeq2["type"] == "PC")]
miRs = list(df_DESeq2_miR["miRNA"].unique())
filtered_miRs = []
for miR in miRs:
    df_miR = df_DESeq2_miR[df_DESeq2_miR["miRNA"] == miR]
    if any(n < 0.05 for n in df_miR["padj"].to_list()):
        filtered_miRs.append(miR)


df = pd.read_excel("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/Log2FC_miR_shRNAone_byType.xlsx")
df_heatmap = df[df["miRNA"].isin(filtered_miRs)][["Log2FC_miR1U", "Log2FC_miR1A", "Log2FC_miR_7merSupp", "Log2FC_miR_7mer", "Log2FC_miR_6mer"]]
df_ref = df[df["miRNA"].isin(filtered_miRs)][["miRNA"]]

df_DESeq2_miR_filtered = df_DESeq2_miR[df_DESeq2_miR["miRNA"].isin(filtered_miRs)]

def get_type(t):
    return df_DESeq2_miR_filtered[df_DESeq2_miR_filtered["type"] == t]

highlighted_miRs = [["Hsa-Mir-935-v1_3p", "Hsa-Mir-430-P5_5p", "Hsa-Mir-154-P18_3p", "Hsa-Mir-8-P2a_3p", 
                    "Hsa-Mir-1914_5p", "Hsa-Mir-3140_3p", "Hsa-Mir-541_3p", "Hsa-Mir-5579_5p"],
                    ["Hsa-Mir-935-v1_3p", "Hsa-Mir-430-P5_5p", "Hsa-Mir-154-P18_3p", "Hsa-Mir-8-P2a_3p", 
                    "Hsa-Mir-1914_5p", "Hsa-Mir-3140_3p", "Hsa-Mir-541_3p", "Hsa-Mir-5579_5p"],
                    ["Hsa-Mir-935-v1_3p", "Hsa-Mir-430-P5_5p", "Hsa-Mir-154-P18_3p", "Hsa-Mir-8-P2a_3p", 
                    "Hsa-Mir-1914_5p", "Hsa-Mir-3140_3p", "Hsa-Mir-541_3p", "Hsa-Mir-5579_5p"],
                    ["Hsa-Mir-935-v1_3p", "Hsa-Mir-430-P5_5p", "Hsa-Mir-154-P18_3p", "Hsa-Mir-8-P2a_3p", 
                    "Hsa-Mir-1914_5p", "Hsa-Mir-3140_3p", "Hsa-Mir-541_3p", "Hsa-Mir-5579_5p"],
                    ["Hsa-Mir-935-v1_3p", "Hsa-Mir-134_5p", "Hsa-Mir-197_3p", "Hsa-Mir-430-P2_3p", "Hsa-Mir-744_5p", "Hsa-Mir-154-P18_3p"],]

pairs = [
    ("miR1U", "miR1A", highlighted_miRs[0]),
    ("miR1U", "miR_7merSupp", highlighted_miRs[1]),
    ("miR1U", "miR_7mer", highlighted_miRs[2]),
    ("miR1U", "miR_6mer", highlighted_miRs[3]),
    ("miR_7merSupp", "miR_7mer", highlighted_miRs[4])
]

p_values = []
sig_miR = []

for type_x, type_y, highlight in pairs:
    dfx = get_type(type_x)[["miRNA", "log2FoldChange", "lfcSE"]]
    dfx["miRNA"] = pd.Categorical(dfx["miRNA"], categories=filtered_miRs, ordered=True)
    dfx = dfx.sort_values("miRNA").reset_index(drop=True)
    dfy = get_type(type_y)[["miRNA", "log2FoldChange", "lfcSE"]]
    dfy["miRNA"] = pd.Categorical(dfy["miRNA"], categories=filtered_miRs, ordered=True)
    dfy = dfy.sort_values("miRNA").reset_index(drop=True)

    # merge on miRNA identity, not row position -- this is the actual fix
    merged = dfx.merge(dfy, on="miRNA", suffixes=("_x", "_y"), how="inner")

    # drop rows where DESeq2 couldn't fit an estimate (NA log2FC/SE)
    merged = merged.dropna(subset=["log2FoldChange_x", "lfcSE_x",
                                   "log2FoldChange_y", "lfcSE_y"])

    merged["delta_FC"] = merged["log2FoldChange_x"] - merged["log2FoldChange_y"]
    merged["sqrt_SE"] = np.sqrt(merged["lfcSE_x"]**2 + merged["lfcSE_y"]**2)
    merged["z"] = merged["delta_FC"] / merged["sqrt_SE"]
    merged["p"] = 2 * stats.norm.sf(np.abs(merged["z"]))
    rejected, p_corrected, _, _ = multipletests(merged["p"].fillna(1.0), alpha=0.05, method='fdr_bh')
    merged["fdr"] = p_corrected
    merged.to_excel(f"/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/pairwise_comparison_{type_x}_vs_{type_y}.xlsx", index=False)
    merged["neglog10p"] = -np.log10(merged["fdr"])

    # clip values so points beyond the display range pile up at the edge
    # instead of being cut off / excluded from the plot
    X_LIM = 1.1
    Y_LIM = 26
    merged["delta_FC_clipped"] = np.clip(merged["delta_FC"], -X_LIM+0.1, X_LIM-0.1)
    merged["neglog10p_clipped"] = np.clip(merged["neglog10p"], -1, Y_LIM - 1)  # clip to 25 max


    merged_sig = merged[(merged["neglog10p"] > -np.log10(0.05)) & (np.abs(merged["delta_FC"]) > 0.5)]
    print(f'there are {len(merged_sig)} significant miRNAs')
    sig_miR.append(merged_sig["miRNA"].to_list())
    p_values.append(merged_sig["neglog10p_clipped"].to_list())
    merged_unsig = merged[(merged["neglog10p"] <= -np.log10(0.05)) | (np.abs(merged["delta_FC"]) <= 0.5)]
    merged_sig_validated = merged_sig[merged_sig["miRNA"].isin(highlight)]
    merged_unsig_validated = merged_unsig[merged_unsig["miRNA"].isin(highlight)]
    merged_validated = merged[merged["miRNA"].isin(highlight)]

    plt.figure(figsize=(6, 6))

    plt.xlim(-X_LIM, X_LIM)
    plt.ylim(-1, Y_LIM)
    plt.axhline(-np.log10(0.05), color="red", linestyle="--", linewidth=1)
    plt.axvline(-0.5, color="grey", linestyle="--", linewidth=1)
    plt.axvline(0.5, color="grey", linestyle="--", linewidth=1)

    sns.scatterplot(x=merged_unsig["delta_FC_clipped"], y=merged_unsig["neglog10p_clipped"], color="lightgrey", edgecolor="black", linewidth=0.5, s=50)
    sns.scatterplot(x=merged_sig["delta_FC_clipped"], y=merged_sig["neglog10p_clipped"], color="coral", edgecolor="black", linewidth=0.5, s=50)
    sns.scatterplot(x=merged_sig_validated["delta_FC_clipped"], y=merged_sig_validated["neglog10p_clipped"], color="coral", edgecolor="black", linewidth=0.5, s=150)
    sns.scatterplot(x=merged_unsig_validated["delta_FC_clipped"], y=merged_unsig_validated["neglog10p_clipped"], color="lightgrey", edgecolor="black", linewidth=0.5, s=150)
    plt.xlabel("Delta Log2 Fold Change")
    plt.ylabel("-Log10(p-value)")
    plt.title(f"Pairwise Comparison: {type_x} vs {type_y}")

    texts = []
    for index, row in merged_validated.iterrows():
        texts.append(
            plt.text(
                row["delta_FC_clipped"],
                row["neglog10p_clipped"],
                row["miRNA"],  # <-- Your label column
                size=8,
            )
        )
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
    plt.savefig(f"/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_figure/figures/main_text/pairwise_comparison_{type_x}_vs_{type_y}.pdf", bbox_inches="tight")
    plt.show()
    plt.close()

    print(f"{type_x} vs {type_y}: {len(merged)} miRNAs compared "
          f"(dropped {len(dfx)+len(dfy)-2*len(merged)} unmatched/NA rows)")

p_values
sig_miR
sig_miR_index = []
for comparison in sig_miR:
    indexes = []
    for miR in comparison:
        indexes.append(filtered_miRs.index(miR))
    sig_miR_index.append(indexes)
sig_miR_index

LoL_heatmap = []
for comparison in range(len(sig_miR_index)):
    p_sig = [0]*434
    for p in range(len(sig_miR_index[comparison])):
        p_sig[sig_miR_index[comparison][p]] = p_values[comparison][p]
    LoL_heatmap.append(p_sig)
    
    

df_heatmap = pd.DataFrame(LoL_heatmap)
plt.figure(figsize=(500,7))
colors = [(0.0, "white"), (1, "#bd262a")]
custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", colors)
sns.heatmap(df_heatmap, square=True, cmap=custom_cmap)
plt.title("Heatmap with Missing Values Replaced by 0")
plt.xlabel("Columns")
plt.ylabel("Rows")
plt.savefig("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_figure/figures/main_text/heatmap_p_values.pdf")
plt.show()