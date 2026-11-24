import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_PCSK9_shRNAone_3nt.csv")
df = df[~(df["region"] == "NTC")].astype({"start_nt_NM_174936.3": int})

plt.figure(figsize=(25, 6))

custom_palette = {
    "5'UTR": "#537dbf",
    "CDS": "#adb6be",
    "3'UTR": "#f4e181",
}

sns.scatterplot(
    data=df, 
    x="start_nt_NM_174936.3", 
    y="Log2FC_mean", 
    hue="region", 
    palette=custom_palette, 
    s = 100, 
    edgecolor="black", 
    linewidth = 0.5
)

validated = [3335, 3395, 3638, 704, 2684, 3077, 425, 3590, 3035, 2783, 2774, 911, 3176, 260, 1688, 1067, 1580, 224, 3605]
df_validated = df[df["start_nt_NM_174936.3"].isin(validated)].sort_values(by=['Log2FC_mean'], ascending=True)
df_validated

sns.scatterplot(
    data=df_validated,
    x="start_nt_NM_174936.3",
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
            row["start_nt_NM_174936.3"],
            row["Log2FC_mean"],
            row["start_nt_NM_174936.3"],  # <-- Your label column
            size=20,
        )
    )
adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

plt.savefig("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/main_text/PCSK9_shRNAone_3nt_scatter_region.pdf", dpi=300, bbox_inches='tight')