# Small RNA-seq preprocessing, figure generation, non-redundant mRNA reference construction, and miRNA target-site analysis

This repository contains the code used to:

1. preprocess small RNA-seq FASTQ files and quantify shRNA-derived reads,
2. generate 5′ end heatmaps and small RNA metric figures from the quantified reads,
3. construct a non-redundant human mRNA reference set from UCSC/GENCODE annotations, and
4. search mRNAs for miRNA seed and 3′-supplementary target sites.

## Repository contents

| File | Purpose |
|---|---|
| [`Snakefile_miseq4_v5-no-collapsing`](Snakefile_miseq4_v5-no-collapsing) | Snakemake workflow for small RNA-seq preprocessing, Bowtie2 mapping, SAM filtering, BAM/BED generation, and strand-specific counting. |
| [`remove_multiple_mapped_reads.py`](remove_multiple_mapped_reads.py) | SAM filter called by the Snakemake workflow. Removes multi-mapped reads according to the shRNA-specific FLAG rules used here. |
| [`261123_FIGURE_novaseq49_share_v1.ipynb`](261123_FIGURE_novaseq49_share_v1.ipynb) | Draws 5′ end heatmaps for the seven miseq4 shRNA backbones (144 variants each) and the corresponding small RNA metric box plots and statistics. |
| [`261123_FIGURE_miseq1and2_share_v1.ipynb`](261123_FIGURE_miseq1and2_share_v1.ipynb) | Draws 5′ end heatmaps and per-sample metric bar plots for the Round 1 (21 samples) and Round 2 (12 samples) shRNA libraries. |
| [`260630_FIGURE_250114_salmon_index_knownCanonical.ipynb`](260630_FIGURE_250114_salmon_index_knownCanonical.ipynb) | Builds the non-redundant human mRNA and 3′ UTR tables from UCSC knownCanonical and foldUtr3. The mRNA table is converted to FASTA for Salmon indexing. |
| [`260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb`](260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb) | Scans the non-redundant mRNAs for miR-494 seed and 3′-supplementary target sites. |

## Requirements

- Python 3.8 or later
- Jupyter Notebook or JupyterLab
- `pandas`, `numpy`, `biopython`
- `matplotlib`, `seaborn`, `scipy` (figure notebooks)
- `snakemake`, `cutadapt`, `bowtie2`, `samtools`, `bedtools`
- `seqkit`
- Optional: `salmon`, if you want to build a Salmon index from the generated FASTA

A conda environment can be created with:

```bash
conda create -n mirna -c conda-forge -c bioconda \
  python=3.8 pandas numpy biopython matplotlib seaborn scipy jupyter snakemake \
  cutadapt bowtie2 samtools bedtools seqkit
conda activate mirna
```

The SAM-filtering script uses only the Python standard library.

## Workflow overview

```text
Small RNA-seq FASTQ
  └─ Snakefile_miseq4_v5-no-collapsing
       └─ Bowtie2 BAM/BED files and per-sample count/position tables
            ├─ 261123_FIGURE_novaseq49_share_v1.ipynb
            │    └─ 5′ end heatmaps + metric figures (7 backbones)
            └─ 261123_FIGURE_miseq1and2_share_v1.ipynb
                 └─ 5′ end heatmaps + metric figures (Round 1/2)

UCSC knownCanonical + foldUtr3
  └─ 260630_FIGURE_250114_salmon_index_knownCanonical.ipynb
       ├─ non-redundant mRNA table
       │    └─ 260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb
       │         └─ miRNA seed / 3′-supplementary site table
       ├─ 3′ UTR table
       └─ mRNA FASTA (for Salmon indexing)
```

---

## 1. Small RNA-seq preprocessing

The Snakemake workflow performs the following steps:

1. Trims the 3′ TruSeq small RNA adapter with `cutadapt` and retains only adapter-trimmed reads (`--trimmed-only`).
2. Applies a minimum quality score of 20 and a length filter of 18–26 nt.
3. Removes four random nucleotides from both ends of each read.
4. Maps reads to sample-specific Bowtie2 indexes with no mismatches (`--score-min C,0,-1`) and reports up to two alignments (`-k 2`).
5. Removes multi-mapped reads with `remove_multiple_mapped_reads.py`.
6. Converts alignments to BAM (dropping unmapped reads with `samtools view -F 4`), sorts and indexes them, and converts them to BED.
7. Counts only alignments on the `+` strand and reports per-reference counts and read positions.

### Configure the workflow

Edit the configuration block at the top of `Snakefile_miseq4_v5-no-collapsing`:

```python
RAW_DATA_DIR = '/path/to/fastq'          # default: '..'
BOWTIE2INDEX = '/path/to/bowtie2_indexes/'
SAMPLES = ['sample1', 'sample2']

ADAPTER_3P = 'TGGAATTCTCGGGTGCCAAGG'     # TruSeq Small RNA
MIN_QUAL = 20
MIN_INSERT_LEN = 18
MAX_INSERT_LEN = 26
RANDOM_SEQ_5p = 4
RANDOM_SEQ_3p = -4
```

The workflow expects:

- one gzipped FASTQ file per sample, named `{sample}.fq.gz`, in `RAW_DATA_DIR`;
- one Bowtie2 index per sample, named `{sample}_ref.*`, in `BOWTIE2INDEX` (include the trailing slash in the configured path);
- sample names in `SAMPLES` must match both the FASTQ prefix and the Bowtie2 index prefix;
- an existing `logs/` directory.

Create the output directories before running:

```bash
mkdir -p logs data annotated
```

### Run

Run the commands from the repository root; the default `RAW_DATA_DIR = '..'` is relative to that directory.

```bash
# Check the workflow without executing it
snakemake -s Snakefile_miseq4_v5-no-collapsing -n

# Run with 8 cores
snakemake -s Snakefile_miseq4_v5-no-collapsing --cores 8
```

### Main outputs

| Path | Description |
|---|---|
| `data/{sample}.trima.fasta.gz` | Adapter-trimmed reads. |
| `data/{sample}.insertonly.fasta.gz` | Reads after random-end trimming and length filtering. |
| `data/{sample}.shrna.bt2output.sam` | Bowtie2 alignments. |
| `data/{sample}.shrna.bt2output.uniq.sam` | Alignments after the multi-mapping filter. |
| `data/{sample}.shrna.bt2output.unsorted.bam` | Unsorted BAM file. |
| `data/{sample}.shrna.bt2output.sorted.bam` | Sorted BAM file. |
| `data/{sample}.shrna.bt2output.sorted.bam.bai` | BAM index. |
| `annotated/{sample}.shrna.bed` | BED file of aligned reads. |
| `annotated/{sample}.shrna.count.txt` | Per-reference read counts (`+` strand only). |
| `annotated/{sample}.shrna.position.txt` | Per-read positions (`+` strand only). |

### The `remove_multiple_mapped_reads.py` script

The script is called automatically by the `remove_multiple_mapped_reads` rule. It can also be run manually:

```bash
python3 remove_multiple_mapped_reads.py input.sam output.sam
```

The script preserves SAM headers and applies the following shRNA-specific rules; it is tailored to the Bowtie2 `-k 2` output produced by this workflow.

- keeps alignments with FLAG 0;
- keeps alignments with FLAG 256 only when the same read name also has a FLAG 16 alignment;
- removes other records, including FLAG 4, FLAG 16, and FLAG 272.

Review this logic before applying the script to other assays.

---

## 2. 5′ end heatmaps and small RNA metrics

Two notebooks turn the per-read position tables (from section 1 or equivalent runs) into heatmaps and metric figures. Both write PDF and PNG figures plus the underlying tables as CSV, and both have input/output paths hard-coded at the top of the notebook.

### Metric definitions

For each sample (miseq1/miseq2) or shRNA variant (novaseq49), reads starting at the intended guide and passenger 5′ ends are counted, and a ±6 nt window is centred on each intended 5′ end:

| Metric | Definition |
|---|---|
| Small RNA homogeneity | `100 x N_guide / N_total` |
| Guide heterogeneity | `100 x (N_guide_window - N_guide) / N_guide_window` |
| Passenger retention | `100 x N_passenger_window / (N_guide_window + N_passenger_window)` |
| Outliers | `100 x (N_total - N_guide_window - N_passenger_window) / N_total` |

`N_guide` is the number of reads starting exactly at the intended guide 5′ end; `N_guide_window` and `N_passenger_window` are the numbers of reads starting within ±6 nt of the intended guide and passenger 5′ ends, respectively; `N_total` is the total number of reads for the sample or variant.

### 2.1 miseq1/miseq2: Round 1 and Round 2 samples

[`261123_FIGURE_miseq1and2_share_v1.ipynb`](261123_FIGURE_miseq1and2_share_v1.ipynb) analyses 33 shRNA samples from two position files: 21 Round 1 samples (18 from the first MiSeq run plus three pLKO controls from the second run, corresponding to seven designs with three guide sequences each) and 12 individual Round 2 samples.

- **Input:** `papername.txt`, `miseq1.shrna.position.txt`, and `miseq2.shrna.position.txt`.
- **Heatmaps:** rows are samples and columns are 5′ end positions; colour is the percentage of a sample's reads starting at that position (log scale 1–100%, with 0% and values below 1% in grey). Three heatmaps are drawn: all 33 samples (with a divider between Round 1 and Round 2), Round 1 only, and Round 2 only. Heatmap height is 50% of the previous version (`HEATMAP_HEIGHT_FRACTION = 0.50`).
- **Metric bar plots:** four panels per figure. For Round 1, each design shows three bars, one per guide sequence (light to dark); guides are not averaged and Kong-A is highlighted in red. For Round 2, each bar is one sample and 26Kong-J3 is highlighted in red. No statistical tests are performed.
- **Default output directory:** `miseq1_261124fixed/`.
- **Output files:** `heatmap_reads_all_samples.{pdf,png}`, `heatmap_reads_round1.*`, `heatmap_reads_round2.*`, `metrics_round1_barplots.*`, `metrics_round2_barplots.*`, `miseq12_metrics_per_sample.csv`, `heatmap_percent_all_samples.csv`, `heatmap_percent_round1.csv`, `heatmap_percent_round2.csv`.

### 2.2 novaseq49 (miseq4): shRNA backbone panel

[`261123_FIGURE_novaseq49_share_v1.ipynb`](261123_FIGURE_novaseq49_share_v1.ipynb) analyses the miseq4 panel: the same 144 guide sequences expressed from seven shRNA backbones (1mir33, 2mirGE, 3pGIPZ, 4J3, 6mir26, 8mir16, 9mir16ghg; 7mir26ghg is not included).

- **Input:** `papername2.txt` and the per-sample `{sample}.shrna.position.txt` tables from section 1.
- **Heatmaps:** one per backbone, with the 144 variants as rows (serial numbers 001–144, every 20th labelled) and 5′ end positions as columns; colour is the percentage of a variant's reads at that position (log scale 1–100%; 0% and values below 1% in grey). Dashed red and dotted blue lines mark the intended guide and passenger 5′ ends. Heatmaps are displayed inline only (`SAVE_FIGURES = False`).
- **Metric box plots:** a combined 2 × 2 panel plus four individual plots (homogeneity, heterogeneity, passenger retention, outliers), with every variant overlaid. 6mir26 is the pink reference and is compared with each of the other six backbones by two-sided Mann–Whitney U tests, with Benjamini–Hochberg adjustment across the six comparisons per metric.
- **Pairwise comparison:** the same four metrics for 3pGIPZ vs 6mir26, with raw Mann–Whitney U p values on the plots and BH-adjusted q values in the table.
- **Default output directory:** `heatmap/metrics_5end_v2/`.
- **Output files:** `metrics_5end_combined.{pdf,png}`, `metrics_5end_{homogeneity,heterogeneity,passenger_retention,outliers}.*`, `metrics_5end_schematic.*`, `3pGIPZ_vs_6mir26/metrics_5end_3pGIPZ_vs_6mir26_*`, and the CSV tables (`shRNA_5end_position_percent.csv`, `shRNA_metrics_5end_all_variants.csv`, `shRNA_metrics_5end_plot_data.csv`, `shRNA_metrics_5end_summary.csv`, `shRNA_metrics_5end_stats_vs_6mir26.csv`, and `3pGIPZ_vs_6mir26/shRNA_metrics_5end_stats_3pGIPZ_vs_6mir26.csv`).

---

## 3. Build the non-redundant human mRNA reference

Open `260630_FIGURE_250114_salmon_index_knownCanonical.ipynb` and run the cells in order after editing the input and output paths.

### Input

The notebook expects a CSV table exported from the UCSC Table Browser on 14 January 2025, named `Canonical_hg38_250114_v6.csv` in the original analysis. The table combines the hg38 knownCanonical transcript collection (GENCODE v47) with the foldUtr3 track; it can be assembled from the UCSC `knownCanonical`, `knownGene`, `kgXref`, `knownGeneMrna`, and `foldUtr3` tables.

The notebook uses these columns:

- `hg38.knownCanonical.protein`
- `hg38.foldUtr3.seq`
- `hg38.knownGene.chrom`
- `hg38.kgXref.geneSymbol`
- `hg38.knownCanonical.transcript`
- `hg38.knownGeneMrna.seq`
- `hg38.kgXref.refseq`

Adjust the column names if your UCSC export uses a different prefix.

### Processing steps

1. Discard records lacking an annotated protein-coding gene or an annotated 3′ UTR.
2. Remove transcripts on alternative or variant contigs.
3. Construct identifiers in the form `gene symbol|Ensembl gene ID|Ensembl transcript ID`.
4. Keep one transcript per gene symbol, selecting the isoform with the longest 3′ UTR.
5. Remove exact mRNA sequence duplicates after sorting by gene symbol.
6. Export the mRNA and 3′ UTR tables.

The final reference contains 19,491 non-redundant mRNAs, one per gene symbol. The 3′ UTR sequences are converted to the DNA alphabet (U→T) and upper-cased.

### Output

| File | Columns |
|---|---|
| `Canonical_hg38_260630_mrna.txt` | `id`, `mrna` |
| `Canonical_hg38_260630_utr3.txt` | `id`, `utr3` |
| `Canonical_hg38_260630_mrna.fa` | mRNA FASTA |

The mRNA table can be converted to FASTA with SeqKit:

```bash
seqkit tab2fx -p id Canonical_hg38_260630_mrna.txt > Canonical_hg38_260630_mrna.fa
```

The notebook name includes `salmon_index` because the FASTA is intended for Salmon indexing. The notebook itself does not build the index; this can be done separately:

```bash
salmon index -t Canonical_hg38_260630_mrna.fa -i salmon_hg38_canonical
```

---

## 4. Search miRNA seed and 3′-supplementary sites

Open `260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb` and run the cells in order after editing the input and output paths.

### Input

- `Canonical_hg38_260630_mrna.txt` from the previous notebook.
- `Canonical_hg38_260630_utr3.txt` is loaded by the notebook but is not used in the current version.

Two 19-nt variants of miR-494-3p are analysed:

| Name | Sequence (5′→3′) |
|---|---|
| `miR-494-7merSupp` | `TGAAACATTGTGGGGATTG` |
| `miR-494-7mer` | `TGAAACATTGTGCCCTTTG` |

The sequences are written in the DNA alphabet (T instead of U).

### Site definitions

miRNA nucleotides are numbered from the 5′ end. Target motifs are written as reverse complements in the mRNA 5′→3′ orientation.

Seed motifs:

- **8mer:** `rc(m2–m8) + A`
- **7mer-m8:** `rc(m2–m8)`
- **7mer-A1:** `rc(m2–m7) + A`
- **6mer:** `rc(m2–m7)`

where `A` is the target nucleotide opposite miRNA position 1.

3′-supplementary motifs:

- `rc(m12–m15)`
- `rc(m13–m16)`
- `rc(m14–m17)`

A composite 3′-supplementary site requires a supplementary motif followed by a 1–20 nt spacer and a 6mer seed match in the mRNA 5′→3′ orientation. This is detected with the regular expression:

```text
(SUPPL.{1,20}SEED6)
```

where `SUPPL` is one of the three supplementary motifs and `SEED6` is the 6mer motif (`TGTTTC` for both variants). For example, the 13–16 window of `miR-494-7merSupp` is searched with:

```text
(TCCC.{1,20}TGTTTC)
```

Matching is exact; no mismatches, G:U wobble, or degenerate positions are allowed. The 8mer, 7mer-m8, and 7mer-A1 seed sites are also recorded separately.

### Output

The notebook writes:

```text
nonredun_mir494_seedsuppl_in_mRNA_260630.csv
```

The table contains one row per mRNA (19,491 rows) and 28 columns: the identifier, the mRNA sequence, and motif/presence columns for each miRNA variant. A period (`.`) indicates no match; multiple matches are comma-separated.

### 3′ UTR-only search

The current code searches the full-length mRNA sequence (`mrna` column). To restrict the search to 3′ UTRs, merge the 3′ UTR table by `id` first (for example, `dfChosen = dfMrna.merge(dfUtr3, on='id')`), then use `x['utr3']` instead of `x['mrna']` in the `contains_seed` and `contains_suppl_region` calls.

---

## Notes

- All file paths in the Snakefile and notebooks are hard-coded to the original analysis environment. Edit them before running.
- The notebooks were developed with Python 3.8.16.
- The two figure notebooks depend on the position tables produced by the Snakemake workflow (or equivalent tables); the mRNA target-site notebooks are independent analyses.
- Figure export is controlled by `SAVE_FIGURES` / `SAVE_METRICS`; in `261123_FIGURE_novaseq49_share_v1.ipynb` the backbone heatmaps are left inline (`SAVE_FIGURES = False`) while the metric figures and tables are exported.
- Manuscript methods are provided in `Methods_nonredundant_mRNA.md` and `Methods_miRNA_seed_supplementary.md`.

## Troubleshooting

- **Missing `logs/` directory:** create it with `mkdir -p logs`; the Snakefile writes run logs there.
- **Bowtie2 index not found:** each sample requires an index named `{sample}_ref.*` in `BOWTIE2INDEX`.
- **Wrong adapter:** `ADAPTER_3P` is set for TruSeq small RNA libraries; change it for other kits.
- **Snakemake indentation error:** the `remove_multiple_mapped_reads` rule contains a tab before `shell:`. If your Python/Snakemake version reports inconsistent indentation, replace that tab with four spaces.

## Citation

If you use this code, please cite the associated manuscript.
