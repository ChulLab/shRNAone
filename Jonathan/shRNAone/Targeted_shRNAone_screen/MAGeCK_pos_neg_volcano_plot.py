import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')


# Phylo color palette
COLOR_NS   = '#BDBDBD'   # non-significant
COLOR_POS  = '#0279EE'   # positive selection (resistance)
COLOR_NEG  = '#FF9400'   # negative selection (sensitivity)
COLOR_HL   = '#FD9BED'   # highlighted hit


# ---- Adjustable parameters ----
FDR_THRESHOLD   = 0.05
LFC_THRESHOLD   = 0.5
HIGHLIGHT_GENES = ['RPLP0']
TOP_N_LABEL     = 10

print("Styling and parameters loaded.")
print(f"  FDR threshold: {FDR_THRESHOLD}")
print(f"  LFC threshold: {LFC_THRESHOLD}")
print(f"  Highlight genes: {HIGHLIGHT_GENES}")

def load_mageck_rra(filepath):
    df = pd.read_csv(filepath, sep='\t')
    # Normalize column names (MAGeCK uses 'pos|lfc', 'pos|fdr', etc.)
    col_map = {}
    for c in df.columns:
        cl = c.lower().replace(' ', '')
        if cl == 'pos|lfc':  col_map[c] = 'pos_lfc'
        elif cl == 'pos|fdr': col_map[c] = 'pos_fdr'
        elif cl == 'neg|lfc': col_map[c] = 'neg_lfc'
        elif cl == 'neg|fdr': col_map[c] = 'neg_fdr'
    df = df.rename(columns=col_map)
    return df
Directories = ["/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/WT/WT.gene_summary.txt", \
               "/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/599/599.gene_summary.txt", \
               "/home/wkkongaa/data/notebooks/260416_novaseq99/analysis/MAGeCK_output/600/600.gene_summary.txt"]
names = ["WT", "miR494_7merSupp", "miR494_7mer"]

for dir, name in zip(Directories, names):
    df = load_mageck_rra(dir)
    df['combined_lfc'] = np.where(df['pos_fdr'] <= df['neg_fdr'],
                                df['pos_lfc'], df['neg_lfc'])
    df['combined_fdr'] = np.minimum(df['pos_fdr'], df['neg_fdr'])
    df['neg_log10_fdr'] = -np.log10(df['combined_fdr'] + 1e-10)

    # Significance flags
    df['sig_pos'] = (df['pos_fdr'] < FDR_THRESHOLD) & (df['pos_lfc'] > LFC_THRESHOLD)
    df['sig_neg'] = (df['neg_fdr'] < FDR_THRESHOLD) & (df['neg_lfc'] < -LFC_THRESHOLD)
    df['sig_combined_pos'] = (df['combined_fdr'] < FDR_THRESHOLD) & (df['combined_lfc'] > LFC_THRESHOLD)
    df['sig_combined_neg'] = (df['combined_fdr'] < FDR_THRESHOLD) & (df['combined_lfc'] < -LFC_THRESHOLD)


    fig, (ax_pos, ax_neg) = plt.subplots(1, 2, figsize=(13, 6.9))

    # --- Positive selection panel ---
    neg_log10_pos = -np.log10(df['pos_fdr'] + 1e-10)
    sig_p = df['sig_pos']

    ax_pos.scatter(df.loc[~sig_p, 'pos_lfc'], neg_log10_pos[~sig_p],
                c=COLOR_NS, s=10, alpha=0.3)
    ax_pos.scatter(df.loc[sig_p, 'pos_lfc'], neg_log10_pos[sig_p],
                c=COLOR_POS, s=18, alpha=0.6)

    # Highlight RPLP0
    row = df[df['id'] == 'RPLP0']
    if len(row) > 0:
        ax_pos.scatter(row['pos_lfc'], -np.log10(row['pos_fdr'] + 1e-10),
                    c=COLOR_HL, s=120, edgecolors='black', linewidth=1.5, zorder=5)
        ax_pos.annotate('RPLP0', (row['pos_lfc'].values[0], -np.log10(row['pos_fdr'].values[0] + 1e-10)),
                        textcoords="offset points", xytext=(12, 5),
                        fontsize=11, fontweight='bold', color=COLOR_HL)

    # Auto-label top pos hits
    for _, r in df[df['sig_pos']].nlargest(5, 'pos_lfc').iterrows():
        if r['id'] != 'RPLP0':
            ax_pos.annotate(r['id'], (r['pos_lfc'], -np.log10(r['pos_fdr'] + 1e-10)),
                            textcoords="offset points", xytext=(6, 3), fontsize=6, alpha=0.7)

    ax_pos.axhline(-np.log10(FDR_THRESHOLD), color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax_pos.axvline(LFC_THRESHOLD, color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax_pos.set_xlabel('pos|lfc (log2 fold change)', fontsize=11)
    ax_pos.set_ylabel('-log10(pos|fdr)', fontsize=11)
    ax_pos.set_xlim(-2, 2)
    ax_pos.set_ylim(-0.1, 3)
    ax_pos.set_title('Positive Selection\n(knockout → resistance)', fontsize=11, fontweight='bold', color=COLOR_POS)

    # --- Negative selection panel ---
    neg_log10_neg = -np.log10(df['neg_fdr'] + 1e-10)
    sig_n = df['sig_neg']

    ax_neg.scatter(df.loc[~sig_n, 'neg_lfc'], neg_log10_neg[~sig_n],
                c=COLOR_NS, s=10, alpha=0.3)
    ax_neg.scatter(df.loc[sig_n, 'neg_lfc'], neg_log10_neg[sig_n],
                c=COLOR_NEG, s=18, alpha=0.6)

    # Auto-label top neg hits
    for _, r in df[df['sig_neg']].nsmallest(5, 'neg_lfc').iterrows():
        if r['id'] != 'PHGDH':
            ax_neg.annotate(r['id'], (r['neg_lfc'], -np.log10(r['neg_fdr'] + 1e-10)),
                            textcoords="offset points", xytext=(6, 3), fontsize=6, alpha=0.7)

    ax_neg.axhline(-np.log10(FDR_THRESHOLD), color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax_neg.axvline(-LFC_THRESHOLD, color='red', linestyle='--', alpha=0.4, linewidth=1)
    ax_neg.set_xlabel('neg|lfc (log2 fold change)', fontsize=11)
    ax_neg.set_ylabel('-log10(neg|fdr)', fontsize=11)
    ax_neg.set_xlim(-2, 2)
    ax_neg.set_ylim(-0.1, 3)
    ax_neg.set_title('Negative Selection\n(knockout → sensitivity)', fontsize=11, fontweight='bold', color=COLOR_NEG)

    plt.suptitle('MAGeCK RRA Two-Panel Volcano Plots', fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'/home/wkkongaa/data/notebooks/260612_miR_screen_reanalysis/_figure/figures/main_text/MAGeCK_two_panel_volcano_{name}.pdf')
    plt.show()
