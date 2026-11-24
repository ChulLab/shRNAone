import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from adjustText import adjust_text

df_TTR_shRNAone_1nt = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_PCSK9_shRNAone_3nt.csv")
gene = ["PCSK9_shRNAone_3nt"]

validated = [[3335, 3395, 3638, 704, 2684, 3077, 425, 3590, 3035, 2783, 2774, 911, 3176, 260, 1688, 1067, 1580, 224, 3605]]
for transcript, df_gene, valid in zip(gene, [df_PCSK9_shRNAone_3nt], validated):
    df_gene = df_gene.sort_values(by=['Log2FC_mean'], ascending=False)
    df_gene["rank"] = [i + 1 for i in range(len(df_gene))]
    plt.figure(figsize=(4, 4))
    sns.set_theme(style="ticks")
    sns.scatterplot(
        data=df_gene, 
        x="rank", 
        y="Log2FC_mean", 
        s = 12, 
        color="coral", 
        alpha=0.5, 
        edgecolor="#ffffff", 
        linewidth = 0.4
        )
    
    custom_palette = {
        "Misannotated": "#bd262a",
        "5'UTR": "#537dbf",
        "CDS": "#adb6be",
        "3'UTR": "#f4e181",
    }

    if transcript == "TTR_shRNAone_5nt" or transcript == "TTR_shRNAone_1nt":
        df_gene = df_gene[~(df_gene["region"] == "NTC")].astype({"start_nt_NM_000371.3": int})
        df_valid = df_gene[df_gene["start_nt_NM_000371.3"].isin(valid)].sort_values(by=['Log2FC_mean'], ascending=True)
        print(df_valid)
        sns.scatterplot(
            data=df_valid,
            x="rank",
            y="Log2FC_mean",
            hue="region",
            palette=custom_palette,
            s=36,
            edgecolor="black",
            linewidth=0.5,
        )
        texts = []
        for index, row in df_valid.iterrows():
            texts.append(
                plt.text(
                    row["rank"],
                    row["Log2FC_mean"],
                    row["start_nt_NM_000371.4"],  # <-- Your label column
                    size=9,
                )
            )
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
    elif transcript == "PCSK9_shRNAone_3nt":
        df_gene = df_gene[~(df_gene["region"] == "NTC")].astype({"start_nt_NM_174936.3": int})
        df_valid = df_gene[df_gene["start_nt_NM_174936.3"].isin(valid)].sort_values(by=['Log2FC_mean'], ascending=True)
        print(df_valid)
        sns.scatterplot(
            data=df_valid,
            x="rank",
            y="Log2FC_mean",
            hue="region",
            palette=custom_palette,
            s=36,
            edgecolor="black",
            linewidth=0.5,
        )
        texts = []
        for index, row in df_valid.iterrows():
            texts.append(
                plt.text(
                    row["rank"],
                    row["Log2FC_mean"],
                    row["start_nt_NM_174936.3"],  # <-- Your label column
                    size=9,
                )
            )
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
    

    # plt.savefig(f"/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/main_text/{transcript}_ranked_scatter.pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    region_matrix = pd.crosstab(df_gene["region"], df_gene["rank"])
    region_order = ["5'UTR", "CDS", "3'UTR", "Misannotated", "NTC"]

    data_list = []
    for type in region_order:
        data_list.append([1 if region == type else 0 for region in df_gene["region"]])    

    for l in data_list:
        print(l.count(1))
    plt.figure(figsize=(16, 2))

    custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_gradient", [(0.0, "#ffffff"),(1.0, "#bd262a")])

    sns.heatmap(
    data=data_list,
    cmap = custom_cmap,
    yticklabels = region_order,
    rasterized=False
    )
    plt.title("Heatmap from List of Lists", fontsize=14, pad=10)
    plt.tight_layout()
    # plt.savefig(f"/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/main_text/{transcript}_heatmap_region.pdf", dpi=300, bbox_inches='tight')
    plt.show()

    