# Small RNA-seq preprocessing, non-redundant mRNA reference construction, and miRNA target-site analysis

This repository contains the code used to:

1. preprocess small RNA-seq FASTQ files and quantify shRNA-derived reads,
2. construct a non-redundant human mRNA reference set from UCSC/GENCODE annotations, and
3. search mRNAs for miRNA seed and 3′-supplementary target sites.

## Repository contents

| File                                                                                                               | Purpose                                                                                                                                                |
| ------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [`Snakefile_miseq4_v5-no-collapsing`](Snakefile_miseq4_v5-no-collapsing)                                           | Snakemake workflow for small RNA-seq preprocessing, Bowtie2 mapping, SAM filtering, BAM/BED generation, and strand-specific counting.                  |
| [`remove_multiple_mapped_reads.py`](remove_multiple_mapped_reads.py)                                               | SAM filter called by the Snakemake workflow. Removes multi-mapped reads according to the shRNA-specific FLAG rules used here.                          |
| [`260630_FIGURE_250114_salmon_index_knownCanonical.ipynb`](260630_FIGURE_250114_salmon_index_knownCanonical.ipynb) | Builds the non-redundant human mRNA and 3′ UTR tables from UCSC knownCanonical and foldUtr3. The mRNA table is converted to FASTA for Salmon indexing. |
| [`260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb`](260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb)         | Scans the non-redundant mRNAs for miR-494 seed and 3′-supplementary target sites.                                                                      |

## Requirements

- Python 3.8 or later
- Jupyter Notebook or JupyterLab
- `pandas`, `numpy`, `biopython`
- `snakemake`, `cutadapt`, `bowtie2`, `samtools`, `bedtools`
- `seqkit`
- Optional: `salmon`, if you want to build a Salmon index from the generated FASTA

A conda environment can be created with:

```bash
conda create -n mirna -c conda-forge -c bioconda \
  python=3.8 pandas numpy biopython jupyter snakemake \
  cutadapt bowtie2 samtools bedtools seqkit
conda activate mirna
```

The SAM-filtering script uses only the Python standard library.

## Workflow overview

```text
Small RNA-seq FASTQ
  └─ Snakefile_miseq4_v5-no-collapsing
       └─ Bowtie2 BAM/BED files and per-sample count/position tables

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

Create the log directory before running:

```bash
mkdir -p logs
```

### Run

Run the commands from the repository root; the default `RAW_DATA_DIR = '..'` is relative to that directory.

```bash
# Run with 8 cores
snakemake -s Snakefile_miseq4_v5-no-collapsing --cores 8
```

### Main outputs

| Path                                           | Description                                           |
| ---------------------------------------------- | ----------------------------------------------------- |
| `data/{sample}.trima.fasta.gz`                 | Adapter-trimmed reads.                                |
| `data/{sample}.insertonly.fasta.gz`            | Reads after random-end trimming and length filtering. |
| `data/{sample}.shrna.bt2output.sam`            | Bowtie2 alignments.                                   |
| `data/{sample}.shrna.bt2output.uniq.sam`       | Alignments after the multi-mapping filter.            |
| `data/{sample}.shrna.bt2output.unsorted.bam`   | Unsorted BAM file.                                    |
| `data/{sample}.shrna.bt2output.sorted.bam`     | Sorted BAM file.                                      |
| `data/{sample}.shrna.bt2output.sorted.bam.bai` | BAM index.                                            |
| `annotated/{sample}.shrna.bed`                 | BED file of aligned reads.                            |
| `annotated/{sample}.shrna.count.txt`           | Per-reference read counts (`+` strand only).          |
| `annotated/{sample}.shrna.position.txt`        | Per-read positions (`+` strand only).                 |

### The `remove_multiple_mapped_reads.py` script

The script is called automatically by the `remove_multiple_mapped_reads` rule. It can also be run manually:

```bash
python3 remove_multiple_mapped_reads.py input.sam output.sam
```

The script preserves SAM headers and applies the following shRNA-specific rules; it is tailored to the Bowtie2 `-k 2` output produced by this workflow.

- keeps alignments with FLAG 0;
- keeps alignments with FLAG 256 only when the same read name also has a FLAG 16 alignment;
- removes other records, including FLAG 4, FLAG 16, and FLAG 272.

---

## 2. Build the non-redundant human mRNA reference

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

| File                             | Columns      |
| -------------------------------- | ------------ |
| `Canonical_hg38_260630_mrna.txt` | `id`, `mrna` |
| `Canonical_hg38_260630_utr3.txt` | `id`, `utr3` |
| `Canonical_hg38_260630_mrna.fa`  | mRNA FASTA   |

The mRNA table can be converted to FASTA with SeqKit:

```bash
seqkit tab2fx -p id Canonical_hg38_260630_mrna.txt > Canonical_hg38_260630_mrna.fa
```

The notebook name includes `salmon_index` because the FASTA is intended for Salmon indexing. The notebook itself does not build the index; this can be done separately:

```bash
salmon index -t Canonical_hg38_260630_mrna.fa -i salmon_hg38_canonical
```

---

## 3. Search miRNA seed and 3′-supplementary sites

Open `260630_FIGURE_250714-novaseq50-suppl_in_mRNA.ipynb` and run the cells in order after editing the input and output paths.

### Input

- `Canonical_hg38_260630_mrna.txt` from the previous notebook.
- `Canonical_hg38_260630_utr3.txt` is loaded by the notebook but is not used in the current version.

Two 19-nt variants of miR-494-3p are analysed:

| Name               | Sequence (5′→3′)      |
| ------------------ | --------------------- |
| `miR-494-7merSupp` | `TGAAACATTGTGGGGATTG` |
| `miR-494-7mer`     | `TGAAACATTGTGCCCTTTG` |

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
- The small RNA-seq workflow and the mRNA target-site notebooks are independent analyses.
- Manuscript methods are provided in `Methods_nonredundant_mRNA.md` and `Methods_miRNA_seed_supplementary.md`.

## Troubleshooting

- **Missing `logs/` directory:** create it with `mkdir -p logs`; the Snakefile writes run logs there.
- **Bowtie2 index not found:** each sample requires an index named `{sample}_ref.*` in `BOWTIE2INDEX`.
- **Wrong adapter:** `ADAPTER_3P` is set for TruSeq small RNA libraries; change it for other kits.
- **Snakemake indentation error:** the `remove_multiple_mapped_reads` rule contains a tab before `shell:`. If your Python/Snakemake version reports inconsistent indentation, replace that tab with four spaces.

## Documentation and verification

This README was written by DeepSeek and checked by S. Chul Kwon.
