import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import scipy.stats as stats

df = pd.read_csv("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_analysis/output/log2FC_validated.csv")
qPCR_shRNAone = df['shRNAone'].tolist()
log2FC = df['log2FC'].tolist()
PearsonR, PearsonP = pearsonr(qPCR_shRNAone, log2FC)
df = pd.DataFrame({"qPCR_shRNAone": qPCR_shRNAone, "log2FC": log2FC})
plt.figure(figsize=(5, 5))
plt.grid(True, linestyle='--', alpha=0.5)

sns.scatterplot(data=df, x='qPCR_shRNAone', y='log2FC', 
            color= 'coral', edgecolors =  'black')


text_str = f'Pearson r = {PearsonR:.3f}\nPearson p = {PearsonP:.3e}\nn={len(df)}'
plt.gca().text(0.05, 0.95, text_str, transform=plt.gca().transAxes, 
               fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
plt.xlabel("Normalised TTR mRNA transcript level")
plt.ylabel("-log2FC (bin1/bin4)")

plt.title('Correlation of log2FC and knockdown in qPCR')

plt.savefig("/home/wkkongaa/data/notebooks/260608_tiling_screen_reanalysis/_figure/figures/preprocessing/correlations_screen_vs_shRNAone.pdf", bbox_inches='tight')
plt.show()
plt.close()