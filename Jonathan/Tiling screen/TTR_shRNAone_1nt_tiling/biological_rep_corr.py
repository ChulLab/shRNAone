df_screen = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_1nt.csv")
columns = [x for x in df_screen.columns.tolist() if x.endswith("Log2FC")]
df_screen = df_screen[columns]

def annotate_stats(x, y, **kwargs):
    # Drop NaNs dynamically for calculating correlations if any rows are missing data
    mask = x.notna() & y.notna()
    x_clean, y_clean = x[mask], y[mask]

    if len(x_clean) < 2:
        return  # Not enough data points to compute correlation

    # Calculate Pearson and Spearman stats
    p_coef, p_val = pearsonr(x_clean, y_clean)
    s_coef, s_val = spearmanr(x_clean, y_clean)

    # Format text labels (scientific notation for p-values is safest for small numbers)
    stats_text = (
        f"Pearson r: {p_coef:.3f} (p: {p_val:.2e})\n"
        f"Spearman ρ: {s_coef:.3f} (p: {s_val:.2e})"
    )

    # Grab current axes object
    ax = plt.gca()

    # Place text box in the upper-left or lower-right corner depending on preference
    ax.text(
        0.05,
        0.95,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="left",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
    )


# 3. Initialize the PairGrid object
# We use 'vars' to strictly limit the plotting to your 3 target replicate columns
g = sns.PairGrid(df_screen, vars=columns, height=3.5)

# 4. Map plots to the Grid
# Diagonal: Distribution histograms
g.map_diag(sns.histplot, kde=True, color="skyblue")

# Off-diagonal (upper and lower triangles): Scatter plots
g.map_offdiag(sns.scatterplot, alpha=0.6, color="coral", edgecolor="w")

# 5. Map the annotation function to the off-diagonal plots
g.map_offdiag(annotate_stats)

# 6. Aesthetic adjustment & show
plt.suptitle(
    "Pairwise Biological Replicate Correlations (TTR_shRNAone_1nt Screen)",
    y=1.02,
    fontsize=14,
    weight="bold",
)
plt.tight_layout()
plt.savefig("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/preprocessing/correlations_TTR_shRNAone_1nt.pdf", bbox_inches="tight")
plt.show()