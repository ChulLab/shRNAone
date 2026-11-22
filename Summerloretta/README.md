# Figure 5. Comparison of shRNAone, CasRx and PspCas13b tiling profiles reveals guide-architecture-dependent positional offsets.

Notebooks for the Figure 5 panels and the supplementary figures. They compare shRNA, PspCas13b, and CasRx (RfxCas13d with U1-crRNA) activity at 1-nt resolution along the TTR transcript.

## Notebooks

| Notebook                                                                                                                   | Purpose                                                             |
| -------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| [`260728_Figure5_Main_a_dummy.ipynb`](260728_Figure5_Main_a_dummy.ipynb)                                                   | Synthetic dummy panel for Figure 5a.                                |
| [`260724_Figure5_Main_b_supa_shRNA_Psp_CasRx_GLine_Dot.ipynb`](260724_Figure5_Main_b_supa_shRNA_Psp_CasRx_GLine_Dot.ipynb) | Gaussian-smoothed shRNA, PspCas13b, and CasRx day 1/day 5 profiles. |
| [`260724_Figure5_Main_c_d_Offset_shRNA_CasRx.ipynb`](260724_Figure5_Main_c_d_Offset_shRNA_CasRx.ipynb)                     | Offset analysis of shRNA versus CasRx.                              |
| [`260805_Figure5_Main_c_shRNA_Psp_Offset.ipynb`](260805_Figure5_Main_c_shRNA_Psp_Offset.ipynb)                             | Offset analysis of shRNA versus PspCas13b.                          |
| [`260805_Figure5_Main_e.ipynb`](260805_Figure5_Main_e.ipynb)                                                               | shRNA versus PspCas13b profile panels at two offsets.               |
| [`260807_VERSUS_SumL_supplementary_figure.ipynb`](260807_VERSUS_SumL_supplementary_figure.ipynb)                           | VERSUS/SumL supplementary figure panels.                            |

The two `Main_c` notebooks are separate analyses: `260724_..._c_d...` compares shRNA with CasRx, whereas `260805_..._c...` compares shRNA with PspCas13b.

## Input data

The experimental notebooks read three CSV tables:

```text
Log2FC_TTR_shRNAone_1nt_260702_3reps.csv
Log2FC_TTR_PspCas13b_1nt.csv
Log2FC_TTR_RfxCas13d_U1-crRNA_1nt.csv
```

Expected columns:

- `start_nt_NM_000371.3` in all tables
- `Log2FC_mean` in the shRNA and PspCas13b tables
- `rep3_day1_Log2FC` and `rep3_day5_Log2FC` in the CasRx table

Edit the input paths (`shrna_in`, `psp_in`, `u1_in`) and `output_dir` in each experimental notebook before running. The data files are not included.

## Requirements

Python 3.11 with `pandas`, `numpy`, `matplotlib`, `scipy`, and Jupyter.

```bash
conda create -n versus -c conda-forge python=3.11 pandas numpy matplotlib scipy jupyter
conda activate versus
```

## Analysis conventions

- Coordinates are shifted by 110 nt for plotting.
- Profiles are smoothed with a 2-nt Gaussian filter.
- Offset analyses scan −10 to +30 nt and report the best Pearson correlation.

## How to run

1. Open the notebook in Jupyter.
2. Set the input paths and `output_dir`.
3. Run the cells in order.
4. Uncomment the relevant `plt.savefig(...)` line to export a PDF.

## Notes

- Most `savefig` calls are commented out.
- `260805_Figure5_Main_c_shRNA_Psp_Offset.ipynb` currently raises `NameError` because `pdf_path` is undefined; uncomment its definition or remove the save call.
- Several notebooks assign the U1 position column from `psp_df` instead of `u1_df`; correct this if the U1 range must be preserved.
- The input column uses `NM_000371.3`, while the plot label uses `NM_000371.4`.

## Documentation and verification

This README was written by DeepSeek and checked by S. Chul Kwon.

# 
