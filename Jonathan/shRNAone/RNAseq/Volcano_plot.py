
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


df = pd.read_csv("/home/wkkongaa/data/notebooks/250916_RNAseq_analysis_miR-494/DESeq2/results_miR494.csv", index_col=0)


LFC_THRESH = 0.5
SIG_THRESH = 0.05
COLOR_BY   = "padj"     # "padj" (recommended) or "pvalue"
N_LABEL    = 8          # top-N significant genes to label (0 = off)

HIGHLIGHT_GENES = ["RPLP0|ENSG00000089157|ENST00000392514"]   # exact index match(es)
HIGHLIGHT_COLOR = "#FF9400"   # orange


d = df.dropna(subset=["log2FoldChange", "pvalue"]).copy()
min_pos = d["pvalue"][d["pvalue"] > 0].min()
d["pvalue"] = d["pvalue"].clip(lower=min_pos)
d["neglog10p"] = -np.log10(d["pvalue"])

sig = d[COLOR_BY].fillna(1)
d["category"] = "NS"
d.loc[(sig < SIG_THRESH) & (d["log2FoldChange"] >=  LFC_THRESH), "category"] = "Up"
d.loc[(sig < SIG_THRESH) & (d["log2FoldChange"] <= -LFC_THRESH), "category"] = "Down"
palette = {"Up": "#b0202a", "Down": "#0279EE", "NS": "#b7b7b7"}


present = [g for g in HIGHLIGHT_GENES if g in d.index]
missing = [g for g in HIGHLIGHT_GENES if g not in d.index]
if missing:
    print("WARNING: highlight ID(s) not found in results index:", missing)


fig, ax = plt.subplots(figsize=(6, 5))

ns = d[d["category"] == "NS"]
ax.scatter(ns["log2FoldChange"], ns["neglog10p"],
           s=14, color=palette["NS"], edgecolor="black", alpha=0.4,
           rasterized=False, zorder=1, label="NS")

for cat in ["Down", "Up"]:
    sub = d[d["category"] == cat]
    ax.scatter(sub["log2FoldChange"], sub["neglog10p"],
               s=20, color=palette[cat], edgecolor="black", linewidth=0.4,
               alpha=1.0, rasterized=False, zorder=2, label=cat)

ax.axhline(-np.log10(SIG_THRESH), ls="--", lw=0.9, color="#555")
ax.axvline( LFC_THRESH, ls="--", lw=0.9, color="#555")
ax.axvline(-LFC_THRESH, ls="--", lw=0.9, color="#555")


ymax = d["neglog10p"].max() * 1.15
ax.set_xlim(d["log2FoldChange"].min() - 0.2, d["log2FoldChange"].max() + 0.2); ax.set_ylim(-0.5, ymax)

hl = d.loc[present]
ax.scatter(hl["log2FoldChange"], hl["neglog10p"],
           s=50, color=HIGHLIGHT_COLOR, edgecolor="black", linewidth=0.8,
           zorder=5, label="Highlight")

label_ids = list(dict.fromkeys(
    d[d["category"] != "NS"].nsmallest(N_LABEL, "pvalue").index.tolist() + present))
texts = []
for gid in label_ids:
    r = d.loc[gid]
    is_hl = gid in present
    # show only the gene symbol (before first '|') to keep labels short
    lbl = str(gid).split("|")[0]
    texts.append(ax.text(r["log2FoldChange"], r["neglog10p"], lbl,
                         fontsize=9 if is_hl else 8,
                         fontweight="bold" if is_hl else "normal",
                         color=HIGHLIGHT_COLOR if is_hl else "black",
                         clip_on=False))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#888", lw=0.6, shrinkA=5),
            expand_points=(1.5, 1.6), force_text=(0.5, 0.8))

ax.set_xlabel(r"$\log_2$ fold change")
ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$")
ax.set_title("DESeq2 volcano plot")
ax.legend(title=f"{COLOR_BY} < {SIG_THRESH}\n|LFC| ≥ {LFC_THRESH}",
          frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=9)
sns.despine(ax=ax)
plt.tight_layout()

plt.savefig("/home/wkkongaa/data/notebooks/250916_RNAseq_analysis_miR-494/DESeq2/volcano_plot.pdf")
plt.show()
plt.close()