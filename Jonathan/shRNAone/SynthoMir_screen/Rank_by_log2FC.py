import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd_log2FC = pd.read_excel("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/Log2FC_miR_shRNAone_byType.xlsx")
pd_log2FC = pd_log2FC[["miRNA", "Log2FC_miR1U", "Log2FC_miR1A"]]
pd_log2FC["Log2FC_mean"] = pd_log2FC[["Log2FC_miR1U", "Log2FC_miR1A"]].mean(axis=1)
pd_log2FC = pd_log2FC.sort_values(by="Log2FC_mean", ascending=False).copy()
pd_log2FC["Rank"] = range(1, len(pd_log2FC) + 1)
pd_log2FC

interested = [["Hsa-Mir-181-P1a_5p", "Hsa-Mir-181-P1c_5p", "Hsa-Mir-181-P2a_5p", "Hsa-Mir-181-P2c_5p"]]
interested.append(["Hsa-Let-7-P1b_5p", "Hsa-Let-7-P1c_5p", "Hsa-Let-7-P1d_5p", 
                   "Hsa-Let-7-P2a3_5p", "Hsa-Let-7-P2b2_5p", "Hsa-Let-7-P2b3_5p", 
                   "Hsa-Let-7-P2c1_5p", "Hsa-Let-7-P2c2_5p", "Hsa-Let-7-P2c3_5p"])
interested.append(["Hsa-Mir-30-P1a_5p", "Hsa-Mir-30-P1b_5p", "Hsa-Mir-30-P1d_5p", "Hsa-Mir-30-P2a_5p", "Hsa-Mir-30-P2b_5p"])
interested.append(["Hsa-Mir-430-P1_3p", "Hsa-Mir-430-P2_3p", "Hsa-Mir-430-P3_3p", "Hsa-Mir-430-P4_3p"])
interested.append(["Hsa-Mir-34-P1_5p", "Hsa-Mir-34-P2a_5p", "Hsa-Mir-34-P2b_5p", "Hsa-Mir-34-P3b_5p", "Hsa-Mir-34-P3a-v1_5p"])
interested_label = ["miR-181 family", "let-7 family", "miR-30 family", "miR-302 family", "miR-34 family"]
highlight_colors = ["#bc2429", "#8cb368", "#f4e182", "#0fb4d6", "#5b8e7d"]

plt.figure(figsize=(6, 6))
plt.scatter(
    pd_log2FC["Rank"], 
    pd_log2FC["Log2FC_mean"], 
    color="coral", 
    edgecolor="w",
    linewidth=0.2,
    s=20, 
    alpha=0.4
)

for idx in range(len(interested)):
    subset = pd_log2FC[pd_log2FC["miRNA"].isin(interested[idx])]
    plt.scatter(
        subset["Rank"], 
        subset["Log2FC_mean"], 
        color=highlight_colors[idx], 
        s=60, 
        edgecolor="black",
        linewidth=0.5,
        alpha=1,
        label=interested_label[idx]
    )
    # for _, row in subset.iterrows():
    #     plt.text(
    #         row["Rank"], 
    #         row["Log2FC_mean"], 
    #         row["miRNA"],
    #         fontsize=8, 
    #         weight="bold", 
    #         color=highlight_colors[idx],
    #         ha='right',
    #         va='bottom'
    #     )
plt.title("Ranked Log2 Fold Change of miRNA Variants\n(shRNAone Screen)", fontsize=14, weight="bold", pad=20)
plt.xlabel("Rank", fontsize=12)
plt.ylabel("Mean Log2 Fold Change", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.4)
plt.savefig("/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_figure/figures/preprocessing/ranked_Log2FC_miR_variants.pdf", bbox_inches="tight")
plt.show()

