# Brian (Guojia Wu) — Computational Analysis of Target RNA Accessibility

## Role

Computational dissection of RNA "super-interference regions" — identification of target-site accessibility as the universal determinant of guide efficacy across shRNAone, CasRx, and PspCas13b systems (Figure 6 and Supplementary Figures SF1–SF2).

## Contents

- `src_manuscript_03/` — Complete computational pipeline (Parts A–F)
  - **Part A**: Data preprocessing, external model benchmarking, offset alignment, spatial permutation controls
  - **Part B**: Combinatorial RNAplfold grid search (W × L × U × direction × region), direction selection, parameter impact quantification
  - **Part C**: U plateau definition, WL performance matrix, siRNA ground-truth validation
  - **Part D**: Optimal parameter selection (W35L20 vs W80L40), cross-system validation, grand performance comparison
  - **Part E**: Permutation-based spatial cluster profiling and accessibility comparison
  - **Part F**: Genome-wide PspCas13b guide accessibility computation (MANE RefSeq v1.5)
  - `Manuscript_FigData/` — Figure-generation scripts for all main and supplementary panels
- `01.Brian_Fig_Legend.md` — Figure legends for Figure 6, SF1, and SF2
- `02.Brian_Method.md` — Computational methods
- `03.1.Brian_result_outline.md` — Results narrative outline
- `03.2.Brian_result_main.md` — Results section draft

## Key findings

- W = 35, L = 20 is the optimal RNAplfold parameter set for predicting knockdown efficacy from 3′-end accessibility
- Accessibility correlates with experimental log2FC across shRNAone (r ≈ 0.33), CasRx (r ≈ 0.50), and PspCas13b (r ≈ 0.38) in the TTR 3′UTR
- Our accessibility model outperforms all published PspCas13b prediction tools, including the sole PspCas13b-specific tool (Fareh2024, r ≈ 0.09)
- Released the first genome-wide, accessibility-ranked PspCas13b guide set for the human transcriptome

## Contact

Guojia Wu — School of Biomedical Sciences, LKS Faculty of Medicine, The University of Hong Kong
