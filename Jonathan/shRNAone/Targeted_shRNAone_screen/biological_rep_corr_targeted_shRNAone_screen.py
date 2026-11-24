import pandas as pd
import numpy as np

dirs = ["/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/WT/WT.sgrna_summary.txt",
        "/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/599/599.sgrna_summary.txt",
        "/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/600/600.sgrna_summary.txt"]
output_file_name = ["wildtype", "miR494_7merSupp", "miR494_7mer"]

plot_pairs = [
    ("log2FC_rep1", "log2FC_rep2"),
    ("log2FC_rep1", "log2FC_rep3"),
    ("log2FC_rep2", "log2FC_rep3"),
]

for dir, out in zip(dirs, output_file_name):
    df = pd.read_table(dir)
    df = df[["sgrna", "Gene", "control_count", "treatment_count"]]
    df["log2FC_rep1"] = [np.log2(float(y.split("/")[0]) / float(x.split("/")[0])) for x, y in zip(df["control_count"], df["treatment_count"])]
    df["log2FC_rep2"] = [np.log2(float(y.split("/")[1]) / float(x.split("/")[1])) for x, y in zip(df["control_count"], df["treatment_count"])]
    df["log2FC_rep3"] = [np.log2(float(y.split("/")[2]) / float(x.split("/")[2])) for x, y in zip(df["control_count"], df["treatment_count"])]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.3))
    min_val = df[["log2FC_rep1", "log2FC_rep2", "log2FC_rep3"]].min().min() - 0.1
    max_val = df[["log2FC_rep1", "log2FC_rep2", "log2FC_rep3"]].max().max() + 0.1

    for i, (x_col, y_col) in enumerate(plot_pairs):
        sns.scatterplot(
            data=df,
            x=x_col,
            y=y_col,
            color = "coral",
            ax=axes[i],
            #alpha=0.7,
            edgecolor = "black",
            s = 40
        )
        corr, p_val = pearsonr(df[x_col], df[y_col])
        p_str = f"{p_val:.2e}" if p_val < 0.001 else f"{p_val:.4f}"
        axes[i].set_xlim(min_val, max_val)
        axes[i].set_ylim(min_val, max_val)
        axes[i].set_title(f"{x_col} vs {y_col}\nr = {corr:.3f}, p = {p_str}")
        axes[i].grid(True, linestyle="--", alpha=0.5)

    # plt.tight_layout()
    plt.savefig(f"/home/wkkongaa/data/notebooks/260416_novaseq99/_figure/correlation_bio_reps_{out}.pdf")
    plt.show()