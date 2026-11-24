"""
Pairwise biological-replicate correlation plots for miRNA screen Log2FC data.

For each screen, this loads an Excel file of Log2FC columns, builds a
seaborn PairGrid (scatter off-diagonal, histogram on the diagonal), overlays
Pearson/Spearman stats on each scatter panel, and optionally highlights a
set of named miRNAs.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import numpy as np

# ====================== CONFIGURATION ======================
ID_COL = "miRNA"  # column in the source Excel file that identifies each row

# Example usage - customize these per screen or globally
# highlight_names = ["Hsa-Mir-483_5p", "Hsa-Mir-493_5p", "Hsa-Mir-9-P1_3p", "Hsa-Mir-28-P3_5p", 
#                    "Hsa-Mir-770_5p", "Hsa-Mir-430-P47_5p", "Hsa-Mir-541_3p", 
#                    "Hsa-Mir-542-v1_3p", "Hsa-Mir-452-v1_5p", "Hsa-Mir-335_3p", "Hsa-Mir-154_P18_3p"]

highlight_names = ["Hsa-Mir-23-P2_3p", "Hsa-Mir-28-P1_5p", "Hsa-Mir-155_5p", "Hsa-Mir-210_3p", "Hsa-Mir-154-P18_3p"]

highlight_colors = ["#bc2429", "#8cb368", "#f4e182", "#f4a259", "#5b8e7d"]#, "#3a86ff", "#ff6fb0", "#00b4d8", "#adb5bd", "#bde0fe", "#000000"]  # one color per name, or a single color for all

# label_names = ["miR-483-5p", "miR-493-5p", "miR-9-3p", "miR-708-5p", 
#                "miR-770-5p", "miR-526b-5p", "miR-541-3p", 
#                "miR-542-3p", "miR-452-5p", "miR-335-3p", "miR-494-3p"]
label_names = ["miR-23b-3p", "miR-28-5p", "miR-155-5p", "miR-210-3p", "miR-494-3p"]

# Optional: size and style of highlighted points
highlight_size = 80
highlight_alpha = 0.95
highlight_edgecolor = "black"
highlight_linewidth = 1.2

dataDirs = ["/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_preprocessing/Log2FC_miR_shRNAone_byType.xlsx"]
screens = ["shRNAone"]
# ===========================================================

if len(highlight_colors) > 1 and len(highlight_colors) != len(highlight_names):
    raise ValueError("highlight_colors must have length 1 or match highlight_names")

if label_names is not None and len(label_names) != len(highlight_names):
    raise ValueError("label_names must be None or match the length of highlight_names")


def annotate_stats(x, y, **kwargs):
    """Add Pearson/Spearman correlation stats as text in the top-left of the axis."""
    mask = x.notna() & y.notna()
    x_clean, y_clean = x[mask], y[mask]
    if len(x_clean) < 2:
        return

    p_coef, p_val = pearsonr(x_clean, y_clean)
    s_coef, s_val = spearmanr(x_clean, y_clean)

    stats_text = (
        f"Pearson r: {p_coef:.3f} (p: {p_val:.2e})\n"
        f"Spearman ρ: {s_coef:.3f} (p: {s_val:.2e})"
    )
    ax = plt.gca()
    ax.text(
        0.05, 0.95, stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )


def scatter_with_highlight(x, y, highlight_df, highlight_names, highlight_colors,
                            label_names=None, id_col=ID_COL, **kwargs):
    """Scatter plot for a PairGrid off-diagonal panel, with selected points highlighted."""
    ax = plt.gca()

    # PairGrid.map_offdiag injects its own 'color' kwarg into every call; drop it
    # so it doesn't collide with the explicit color="coral" below.
    kwargs.pop("color", None)

    # Background points
    sns.scatterplot(
        x=x, y=y,
        color="coral",
        # alpha=0.5,
        edgecolor="black",
        s=40,
        ax=ax,
        **kwargs
    )

    if not highlight_names:
        return

    mask = highlight_df[id_col].isin(highlight_names)
    highlight_subset = highlight_df[mask]
    if highlight_subset.empty:
        return

    # Resolve a color per highlighted row, preserving highlight_names order
    if len(highlight_colors) == 1:
        colors = [highlight_colors[0]] * len(highlight_subset)
    else:
        color_map = dict(zip(highlight_names, highlight_colors))
        colors = highlight_subset[id_col].map(color_map).tolist()

    # Use matplotlib's scatter directly since it accepts a per-point color list
    ax.scatter(
        highlight_subset[x.name],
        highlight_subset[y.name],
        c=colors,
        s=highlight_size,
        alpha=highlight_alpha,
        edgecolor=highlight_edgecolor,
        linewidth=highlight_linewidth,
        zorder=5,
    )

    # Offset labels relative to the axis range so they scale with the data
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    x_off = 0.01 * x_range if x_range else 0
    y_off = 0.01 * y_range if y_range else 0

    # Map each highlighted name to its custom label, falling back to the raw ID
    if label_names:
        label_map = dict(zip(highlight_names, label_names))
    else:
        label_map = {}

    for _, row in highlight_subset.iterrows():
        display_label = label_map.get(row[id_col]) or row[id_col]
        ax.text(
            row[x.name] + x_off, row[y.name] + y_off,
            display_label,
            fontsize=8,
            color="black",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.5),
        )


def plot_screen_correlations(screen_name, xlsx_path, id_col=ID_COL, save=False):
    """Load one screen's Log2FC data and render its pairwise correlation grid."""
    df_screen = pd.read_excel(xlsx_path)

    log2fc_cols = [c for c in df_screen.columns if c.startswith("Log2FC")]
    if len(log2fc_cols) < 2:
        print(f"Skipping {screen_name}: not enough Log2FC columns")
        return None

    df_screen = df_screen[[id_col] + log2fc_cols].copy()

    g = sns.PairGrid(df_screen, vars=log2fc_cols, height=3.8, diag_sharey=False)
    
    g.set(xlim=(df_screen[log2fc_cols].min().min() - 0.2, df_screen[log2fc_cols].max().max() + 0.2),
          ylim=(df_screen[log2fc_cols].min().min() - 0.2, df_screen[log2fc_cols].max().max() + 0.2))
    
    for ax in g.axes.flat:
        if ax is not None:
            ax.set_box_aspect(1)

    g.map_diag(sns.histplot, kde=True, color="skyblue", edgecolor="black")

    g.map_offdiag(
        scatter_with_highlight,
        highlight_df=df_screen,
        highlight_names=highlight_names,
        highlight_colors=highlight_colors,
        label_names=label_names,
        id_col=id_col,
    )
    for i, j in zip(*np.triu_indices_from(g.axes, k=1)):
        # Upper triangle
        g.axes[i, j].axline((0, 0), slope=1, color='gray', linestyle='--', zorder=0)
        # Lower triangle
        g.axes[j, i].axline((0, 0), slope=1, color='gray', linestyle='--', zorder=0)
    g.map_offdiag(annotate_stats)

    g.figure.suptitle(
        f"Pairwise Biological Replicate Correlations - {screen_name} Screen",
        y=1.02,
        fontsize=15,
        weight="bold",
    )

    for ax in g.axes.flat:
        ax.tick_params(labelsize=9)

    g.figure.tight_layout()

    if save:
        g.figure.savefig(f"correlations_{screen_name}.pdf", bbox_inches="tight", dpi=300)

    return g


if __name__ == "__main__":
    for screen, data_path in zip(screens, dataDirs):
        grid = plot_screen_correlations(screen, data_path, save=True)
        if grid is not None:
            plt.show()