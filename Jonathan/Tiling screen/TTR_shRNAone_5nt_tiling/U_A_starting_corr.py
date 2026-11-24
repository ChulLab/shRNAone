import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import scipy.stats as stats

df = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/Log2FC_TTR_shRNAone_5nt.csv")
U_log2FC = df[df["tech_rep"] == "U_starting"]['Log2FC_mean'].tolist()
A_log2FC = df[df["tech_rep"] == "A_starting"]['Log2FC_mean'].tolist()
PearsonR, PearsonP = pearsonr(U_log2FC, A_log2FC)
df = pd.DataFrame({"U_log2FC": U_log2FC, "A_log2FC": A_log2FC})
plt.figure(figsize=(5, 5))
plt.grid(True, linestyle='--', alpha=0.5)

sns.scatterplot(data=df, x='U_log2FC', y='A_log2FC', 
            color= 'coral', edgecolors =  'black')


text_str = f'Pearson r = {PearsonR:.3f}\nPearson p = {PearsonP:.3e}\nn={len(df)}'
plt.gca().text(0.05, 0.95, text_str, transform=plt.gca().transAxes, 
               fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
plt.xlim(min(min(U_log2FC), min(A_log2FC)) - 0.1, max(max(U_log2FC), max(A_log2FC)) + 0.1)
plt.ylim(min(min(U_log2FC), min(A_log2FC)) - 0.1, max(max(U_log2FC), max(A_log2FC)) + 0.1)

plt.title('Correlation of log2FC between A- and U-starting 5nt-tiling shRNA screens')

plt.savefig("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/preprocessing/correlations_TTR_shRNAone_5nt_tech_reps.pdf", bbox_inches='tight')
plt.show()
plt.close()
