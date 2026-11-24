import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text



df = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_1nt_260702_3reps.csv")
df = df[~(df["region"] == "NTC")].astype({'start_nt_NM_000371.3': int})

plt.figure(figsize=(30, 10))

custom_palette = {
    "Misannotated": "#bd262a",
    "5'UTR": "#537dbf",
    "CDS": "#adb6be",
    "3'UTR": "#f4e181",
}

sns.scatterplot(
    data=df, 
    x="start_nt_NM_000371.4", 
    y="Log2FC_mean", 
    hue="region", 
    palette=custom_palette, 
    s = 100, 
    edgecolor="black", 
    linewidth = 0.5
)


validated = [613, 618, 608, 418, 468, 463, 483, 668, 508, 628, 658, 603, 148, 598, 663, 258, 568]
df_validated = df[df["start_nt_NM_000371.3"].isin(validated)].sort_values(by=['Log2FC_mean'], ascending=True)
df_validated

sns.scatterplot(
    data=df_validated,
    x="start_nt_NM_000371.4",
    y="Log2FC_mean",
    hue="region",
    palette=custom_palette,
    s=400,
    edgecolor="black",
    linewidth=0.5,
)
texts = []
for index, row in df_validated.iterrows():
    texts.append(
        plt.text(
            row["start_nt_NM_000371.4"],
            row["Log2FC_mean"],
            row["start_nt_NM_000371.4"],  # <-- Your label column
            size=20,
        )
    )
adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

plt.savefig("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/main_text/TTR_shRNAone_1nt_scatter_region.pdf", dpi=300, bbox_inches='tight')