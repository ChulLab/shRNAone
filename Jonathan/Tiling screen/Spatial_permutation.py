"""
Spatial Permutation Analysis for Tiling shRNA Screen
=====================================================
Corrected version of the original script.

FIXES APPLIED (see review for full rationale):
  [CRITICAL] 1. KeyError on 'region' — column was not included in the
                .select subset before the filter. Now included.
  [CRITICAL] 2. PCSK9 data was labelled "TTR" (copy-paste bug). Fixed.
  [METHOD]   3. Added optional Benjamini-Hochberg FDR correction on
                per-position p-values (controls family-wise error across
                hundreds of positions).
  [METHOD]   4. Added a diagnostic print of the inferred window_rows and
                effective nt span, so you can sanity-check the row→nt
                conversion for gappy tiling designs.
  [MINOR]    5. Relative_Pos now optionally uses a user-supplied transcript
                length instead of max(start).
  [MINOR]    6. Guarded against empty gene group (n=0).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from scipy.ndimage import label as scipy_label
from scipy.stats import false_discovery_control  # BH correction (SciPy >=1.11)

mpl.rcParams['pdf.fonttype'] = 42

# ==========================================
# FILE PATHS  (update to your environment)
# ==========================================
sTtrTable   = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_1nt_260702_3reps.csv'
sPcsk9Table = '/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_PCSK9_shRNAone_3nt.csv'
sOutputFig  = '/home/wkkongaa/data/notebooks/260522_RNAplfold/spatial_permutation/TTR_PCSK9_permutation.pdf'

# ==========================================
# CONFIGURATION
# ==========================================
WINDOW_NT    = 21       # desired window width in nucleotides
N_ITERATIONS = 1000     # permutation iterations
APPLY_FDR    = True     # Benjamini-Hochberg on per-position p-values
# Optional: supply known transcript lengths for accurate Relative_Pos.
#   Set to None to fall back to max(start) per gene (original behaviour).
TRANSCRIPT_LEN_NT = {'TTR': None, 'PCSK9': None}   # e.g. {'TTR': 735, 'PCSK9': 3747}

np.random.seed(42)

C_TTR_PT    = '#f3a359'
C_TTR_LN    = '#000000'
C_PC_PT     = '#f3a359'
C_PC_LN     = '#000000'
C_HIGHLIGHT = '#2ECC71'
C_NOISE     = '#95A5A6'

# ==========================================
# DATA LOADING
# ==========================================
# FIX 1: include 'region' in the column subset so the filter below works.
df_TTR = pd.read_csv(sTtrTable)[
    ['name', 'target_seq', 'region',
     'start_nt_NM_000371.4', 'end_nt_NM_000371.4', 'Log2FC_mean']
]
df_TTR = df_TTR[(df_TTR["region"] != "NTC") & (df_TTR["region"] != "Misannotated")].copy()
df_TTR["Gene"] = "TTR"                                   # FIX: correct label
df_TTR = df_TTR.rename(columns={
    "target_seq": "seq",
    "start_nt_NM_000371.4": "start",
    "end_nt_NM_000371.4": "end",
    "Log2FC_mean": "Log2FC",
})

df_PCSK9 = pd.read_csv(sPcsk9Table)[
    ['name', 'target_seq', 'region',
     'start_nt_NM_174936.3', 'end_nt_NM_174936.3', 'Log2FC_mean']
]
df_PCSK9 = df_PCSK9[(df_PCSK9["region"] != "NTC") & (df_PCSK9["region"] != "Misannotated")].copy()
df_PCSK9["Gene"] = "PCSK9"                               # FIX 2: was "TTR"
df_PCSK9 = df_PCSK9.rename(columns={
    "target_seq": "seq",
    "start_nt_NM_174936.3": "start",
    "end_nt_NM_174936.3": "end",
    "Log2FC_mean": "Log2FC",
})

df = pd.concat([df_TTR, df_PCSK9], ignore_index=True)

# FIX 7: Coerce coordinate / value columns to numeric.
# A single non-numeric token (e.g. "NA", "-", "1,234", a duplicated header
# row) makes pandas read the whole column as object/str, which breaks
# division, .diff(), and .sort_values() (lexicographic sort: "100" < "20").
# errors='coerce' converts unparseable tokens to NaN; we then drop those rows.
for col in ['start', 'end', 'Log2FC']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

n_before = len(df)
df = df.dropna(subset=['start', 'end', 'Log2FC']).reset_index(drop=True)
n_dropped = n_before - len(df)
if n_dropped > 0:
    print(f"  Dropped {n_dropped} row(s) with non-numeric start/end/Log2FC.")

# ==========================================
# 1. PERMUTATION TEST
# ==========================================
df = df.sort_values(by=['Gene', 'start']).reset_index(drop=True)

# FIX 5: Relative_Pos — use transcript length if provided, else max(start).
def _rel_pos(group, gene):
    tx_len = TRANSCRIPT_LEN_NT.get(gene)
    denom = tx_len if tx_len else group['start'].max()
    if denom == 0 or pd.isna(denom):
        return pd.Series(np.zeros(len(group)), index=group.index)
    return group['start'] / denom * 100

df['Relative_Pos'] = np.nan
for gene in ['TTR', 'PCSK9']:
    mask = df['Gene'] == gene
    df.loc[mask, 'Relative_Pos'] = _rel_pos(df[mask], gene).values

processed_data = []
stats_results  = {}

for gene in ['TTR', 'PCSK9']:
    gene_df   = df[df['Gene'] == gene].copy().reset_index(drop=True)
    raw_vals  = gene_df['Log2FC'].values
    n         = len(raw_vals)

    if n == 0:
        print(f"  WARNING [{gene}]: no guides after filtering. Skipping.")
        continue

    # ---- Infer tiling resolution (nt per guide) from consecutive starts ----
    diffs = gene_df['start'].diff().dropna()
    diffs = diffs[diffs > 0]

    if len(diffs) == 0 or diffs.median() == 0:
        print(f"  WARNING [{gene}]: could not infer tiling resolution from 'start' "
              f"(n_guides={n}, unique_starts={gene_df['start'].nunique()}). "
              f"Falling back to window_rows=WINDOW_NT ({WINDOW_NT}).")
        window_rows = WINDOW_NT
    else:
        res         = diffs.median()
        window_rows = max(1, int(round(WINDOW_NT / res)))

    # DIAGNOSTIC: report the conversion so you can sanity-check it
    actual_span = gene_df['start'].max() - gene_df['start'].min()
    avg_nt_per_row = actual_span / max(1, n - 1)
    print(f"  [{gene}] n_guides={n}  res={res if 'res' in dir() else 'N/A'} nt/guide  "
          f"window_rows={window_rows}  "
          f"~{window_rows * avg_nt_per_row:.0f} nt effective span  "
          f"(avg {avg_nt_per_row:.1f} nt/row)")

    # ---- Observed rolling mean ----
    # NOTE: rolling(window=window_rows) operates on CONSECUTIVE ROWS of the
    # sorted DataFrame, not on nucleotide coordinates.  center=True places
    # the current row in the centre; min_periods=1 shrinks the window at the
    # first/last (window_rows//2) rows.
    gene_df['Rolling_Log2FC'] = (
        pd.Series(raw_vals)
        .rolling(window=window_rows, center=True, min_periods=1)
        .mean()
        .values
    )

    # ---- Null model: N_ITERATIONS permutations ----
    simulated_means = np.empty((N_ITERATIONS, n), dtype=float)
    for i in range(N_ITERATIONS):
        scrambled = np.random.permutation(raw_vals)
        sim_roll  = (
            pd.Series(scrambled)
            .rolling(window=window_rows, center=True, min_periods=1)
            .mean()
            .values
        )
        simulated_means[i] = sim_roll
        if i == 0:
            gene_df['Rep_Scrambled_Pts']  = scrambled
            gene_df['Rep_Scrambled_Line'] = sim_roll

    # ---- Position-wise 99% CI ----
    upper_bound = np.nanpercentile(simulated_means, 99.5, axis=0)
    lower_bound = np.nanpercentile(simulated_means,  0.5, axis=0)

    obs = gene_df['Rolling_Log2FC'].values

    # ---- Empirical two-tailed p-value per position ----
    p_vals = np.mean(np.abs(simulated_means) >= np.abs(obs)[np.newaxis, :], axis=0)
    p_vals = np.clip(p_vals, 1 / N_ITERATIONS, 1.0)

    # FIX 3: optional Benjamini-Hochberg FDR correction
    if APPLY_FDR and n > 1:
        p_adj = false_discovery_control(p_vals, method='bh')
    else:
        p_adj = p_vals.copy()

    gene_df['pval']         = p_vals
    gene_df['pval_adj']     = p_adj
    gene_df['neg_log10p']   = -np.log10(p_vals)
    gene_df['neg_log10p_adj'] = -np.log10(p_adj)

    # ---- Significance mask (position-wise CI, two-tailed α=0.01) ----
    sig_above = obs > upper_bound          # NEW: keep the two tails apart
    sig_below = obs < lower_bound
    sig_mask  = sig_above | sig_below      # unchanged, still used downstream
    pct_sig = sig_mask.mean() * 100

    # ---- Cluster summary ----
    # Label each tail separately. An above-band run that happens to touch a
    # below-band run would otherwise fuse into one cluster with no usable
    # direction and a mean_abs_FC that averages opposite-signed excursions.
    starts  = gene_df['start'].to_numpy()          # positional, index-safe
    rel_pos = gene_df['Relative_Pos'].to_numpy()

    cluster_info = []
    for direction, mask in (('above', sig_above), ('below', sig_below)):
        labeled, n_lab = scipy_label(mask)
        for lab in range(1, n_lab + 1):
            idx = np.where(labeled == lab)[0]
            cluster_info.append({
                'cluster': 0,                                 # renumbered below
                'direction': direction,                       # NEW
                'start_nt': int(starts[idx].min()),           # NEW
                'end_nt':   int(starts[idx].max()),           # NEW
                'span_nt':   int(starts[idx].max() - starts[idx].min() + 1),           # NEW
                'n_guides': len(idx),
                'span_pct': round(rel_pos[idx].max() - rel_pos[idx].min(), 2),
                'mean_FC': round(obs[idx].mean(), 3),
                'min_pval': round(p_vals[idx].min(), 4),
                'min_pval_adj': round(p_adj[idx].min(), 4),
                'mean_pval': round(p_vals[idx].mean(), 4),
            })

    df_clusters = pd.DataFrame(cluster_info)
    if not df_clusters.empty:
        df_clusters = df_clusters.sort_values('start_nt').reset_index(drop=True)
        df_clusters['cluster'] = np.arange(1, len(df_clusters) + 1)
        df_clusters = df_clusters[['cluster', 'direction', 'start_nt', 'end_nt', 'span_nt',
                                   'n_guides', 'span_pct', 'mean_FC',
                                   'min_pval', 'min_pval_adj', 'mean_pval']]
    n_clusters = len(df_clusters)

    _, n_fused = scipy_label(sig_mask)
    if n_clusters != n_fused:
        print(f"  note: {n_fused} cluster(s) with the tails merged -> "
              f"{n_clusters} after separating above/below")

    stats_results[gene] = {
        'upper':       upper_bound,
        'lower':       lower_bound,
        'pct_sig':     pct_sig,
        'all_sims':    simulated_means,
        'df_clusters': df_clusters,
        'n_clusters':  n_clusters,
    }

    processed_data.append(gene_df)

    print(f"\n{'='*60}")
    print(f"  {gene}  |  {n_clusters} significant cluster(s)  |  {pct_sig:.1f}% of transcript")
    print(f"{'='*60}")
    if not df_clusters.empty:
        print(df_clusters.to_string(index=False))
    else:
        print("  No significant clusters detected.")

df_clean = pd.concat(processed_data).reset_index(drop=True)

# ==========================================
# 2. DYNAMIC Y-LIMITS
# ==========================================
def gene_ylim(gene_name):
    m = df_clean[df_clean['Gene'] == gene_name]['Log2FC'].abs().max()
    return (-m * 1.15, m * 1.15)

ylim_TTR   = gene_ylim('TTR')
ylim_PCSK9 = gene_ylim('PCSK9')

# ==========================================
# 3. PLOTTING HELPERS  (unchanged logic)
# ==========================================
sns.set_theme(style="ticks")

def plot_spatial_single(ax, gene_df, y_pts_col, y_line_col,
                        pt_color, ln_color, stats, title, y_limits, is_original=True):
    x   = gene_df['Relative_Pos']
    obs = gene_df[y_line_col]

    sns.scatterplot(data=gene_df, x='Relative_Pos', y=y_pts_col, color=pt_color,
                    alpha=0.4, edgecolor="black", s=25, label='Raw Guides', ax=ax, zorder=1)

    ax.fill_between(x, stats['lower'], stats['upper'],
                    color=C_NOISE, alpha=0.15,
                    label=f'99% Null Band (N={N_ITERATIONS})', zorder=2)

    sns.lineplot(data=gene_df, x='Relative_Pos', y=y_line_col, color=ln_color,
                 linewidth=1, label=f'{WINDOW_NT}nt Rolling Mean', ax=ax, zorder=4)

    if is_original:
        ax.fill_between(x, obs, stats['upper'],
                        where=(obs > stats['upper']), color=C_HIGHLIGHT, alpha=0.6,
                        label='Significant Cluster ($p<0.01$)', zorder=5)
        ax.fill_between(x, obs, stats['lower'],
                        where=(obs < stats['lower']), color=C_HIGHLIGHT, alpha=0.6, zorder=5)

        txt = (f"Clusters: {stats['n_clusters']}\n"
               f"Transcript coverage: {stats['pct_sig']:.1f}%")
        ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=10, fontweight='bold',
                va='top', bbox=dict(boxstyle='round,pad=0.4', fc='white',
                                    ec=C_HIGHLIGHT, alpha=0.9), zorder=6)

    ax.set_title(title, loc='left', fontsize=13, fontweight='bold', pad=8)
    ax.set_xlabel('Relative Transcript Position (%)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Log$_2$FC', fontsize=11, fontweight='bold')
    ax.axhline(0, color='black', linewidth=1, zorder=0)
    ax.set_ylim(y_limits);  ax.set_xlim(-2, 102)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(), dict(zip(labels, handles)).keys(),
              loc='lower right', framealpha=0.9, fontsize=9)


def plot_permutation_aggregate(ax, gene_df, stats, ln_color, title, y_limits):
    x     = gene_df['Relative_Pos'].values
    y_obs = gene_df['Rolling_Log2FC'].values

    for i, sim in enumerate(stats['all_sims']):
        ax.plot(x, sim, color='grey', alpha=0.04, linewidth=0.8,
                label=f'{N_ITERATIONS} Permutations' if i == 0 else None, zorder=1)

    ax.plot(x, y_obs, color=ln_color, linewidth=1.5,
            label='Observed Biological Mean', zorder=4)

    ax.plot(x, stats['upper'], color='black', linestyle='--', linewidth=1,
            label='99% CI Boundary', zorder=3)
    ax.plot(x, stats['lower'], color='black', linestyle='--', linewidth=1, zorder=3)

    ax.fill_between(x, y_obs, stats['upper'], where=(y_obs > stats['upper']),
                    color=C_HIGHLIGHT, alpha=0.6, label='Significant Divergence', zorder=5)
    ax.fill_between(x, y_obs, stats['lower'], where=(y_obs < stats['lower']),
                    color=C_HIGHLIGHT, alpha=0.6, zorder=5)

    ax.set_title(title, loc='left', fontsize=13, fontweight='bold', pad=8)
    ax.set_xlabel('Relative Transcript Position (%)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Rolling Mean Log$_2$FC', fontsize=11, fontweight='bold')
    ax.axhline(0, color='black', linewidth=1, zorder=0)
    ax.set_ylim(y_limits);  ax.set_xlim(-2, 102)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(dict(zip(labels, handles)).values(), dict(zip(labels, handles)).keys(),
              loc='lower right', framealpha=0.9, fontsize=9)


def plot_empirical_pval(ax, gene_df, stats, ln_color, title, use_adj=False):
    """
    Panel G/H — empirical p-value track.
    use_adj=False (default): raw per-position p-values.
    use_adj=True: BH-adjusted p-values (controls FDR across positions).
    """
    x    = gene_df['Relative_Pos'].values
    col  = 'neg_log10p_adj' if use_adj else 'neg_log10p'
    nlp  = gene_df[col].values
    pval = gene_df['pval_adj' if use_adj else 'pval'].values

    threshold = -np.log10(0.01)

    ax.fill_between(x, 0, nlp, where=(nlp >= threshold),
                    color=C_HIGHLIGHT, alpha=0.5, label='$p < 0.01$', zorder=3)
    ax.fill_between(x, 0, nlp, where=(nlp < threshold),
                    color=ln_color, alpha=0.25, zorder=2)

    ax.plot(x, nlp, color=ln_color, linewidth=1, zorder=4)

    ax.axhline(threshold, color='black', linestyle='--', linewidth=0.5,
               label=f'$p = 0.01$ ($-\\log_{{10}}$ = {threshold:.1f})', zorder=5)

    n_sig = (nlp >= threshold).sum()
    label_txt = 'FDR' if use_adj else 'raw'
    ax.text(0.02, 0.97,
            f"Guides $p<0.01$ ({label_txt}): {n_sig}/{len(x)} ({n_sig/len(x)*100:.1f}%)",
            transform=ax.transAxes, fontsize=10, fontweight='bold', va='top',
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=ln_color, alpha=0.9), zorder=6)

    ax.set_title(title, loc='left', fontsize=13, fontweight='bold', pad=8)
    ax.set_xlabel('Relative Transcript Position (%)', fontsize=11, fontweight='bold')
    ax.set_ylabel('$-\\log_{10}$ (empirical $p$)', fontsize=11, fontweight='bold')
    ax.set_xlim(-2, 102);  ax.set_ylim(0, nlp.max() * 1.25)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9)

# ==========================================
# 4. RENDER (4 rows x 2 cols)
# ==========================================
fig, axes = plt.subplots(4, 2, figsize=(15, 18))

ttr_df   = df_clean[df_clean['Gene'] == 'TTR'].reset_index(drop=True)
pcsk9_df = df_clean[df_clean['Gene'] == 'PCSK9'].reset_index(drop=True)

# Row 1 — Biological data
plot_spatial_single(axes[0, 0], ttr_df,   'Log2FC', 'Rolling_Log2FC',
                    C_TTR_PT, C_TTR_LN, stats_results['TTR'],
                    'A. TTR: Biological Spatial Distribution', ylim_TTR)
plot_spatial_single(axes[0, 1], pcsk9_df, 'Log2FC', 'Rolling_Log2FC',
                    C_PC_PT,  C_PC_LN,  stats_results['PCSK9'],
                    'B. PCSK9: Biological Spatial Distribution', ylim_PCSK9)

# Row 2 — Single representative null permutation
plot_spatial_single(axes[1, 0], ttr_df,   'Rep_Scrambled_Pts', 'Rep_Scrambled_Line',
                    C_NOISE, 'black', stats_results['TTR'],
                    'C. TTR Null Model: Single Permutation', ylim_TTR, is_original=False)
plot_spatial_single(axes[1, 1], pcsk9_df, 'Rep_Scrambled_Pts', 'Rep_Scrambled_Line',
                    C_NOISE, 'black', stats_results['PCSK9'],
                    'D. PCSK9 Null Model: Single Permutation', ylim_PCSK9, is_original=False)

# Row 3 — All permutations aggregate
plot_permutation_aggregate(axes[2, 0], ttr_df,   stats_results['TTR'],
                           C_TTR_LN, f'E. TTR: Observed vs {N_ITERATIONS} Permutations', ylim_TTR)
plot_permutation_aggregate(axes[2, 1], pcsk9_df, stats_results['PCSK9'],
                           C_PC_LN,  f'F. PCSK9: Observed vs {N_ITERATIONS} Permutations', ylim_PCSK9)

# Row 4 — Empirical p-value track (uses BH-adjusted p-values if APPLY_FDR=True)
plot_empirical_pval(axes[3, 0], ttr_df,   stats_results['TTR'],
                    C_TTR_LN, 'G. TTR: Empirical $p$-value Track', use_adj=APPLY_FDR)
plot_empirical_pval(axes[3, 1], pcsk9_df, stats_results['PCSK9'],
                    C_PC_LN,  'H. PCSK9: Empirical $p$-value Track', use_adj=APPLY_FDR)

sns.despine()
fig.tight_layout(pad=2.5)

# plt.savefig(sOutputFig, dpi=300, bbox_inches='tight')
plt.close()
# print(f"\nFigure saved -> {sOutputFig}")