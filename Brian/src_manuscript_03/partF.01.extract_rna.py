#!/usr/bin/env python3
"""
partF.01.extract_rna.py
解析 MANE GenBank 原始文件，提取包含 CDS 的人类转录本。
"""

import os
import sys
import gzip
import csv
import time
from datetime import datetime

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("FATAL: 'biopython' library is required (pip install biopython).")

# ================= 配置区 =================
BASE_DIR = "/data/cai801/data/HKUcas/result/08.manuscript_03/partF"
GBFF_FILE = os.path.join(BASE_DIR, "00.raw_data", "MANE.GRCh38.v1.5.refseq_rna.gbff.gz")
OUT_CSV = os.path.join(BASE_DIR, "00.raw_data", "transcript_list.csv")

def log_info(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] INFO: {msg}", flush=True)

def process_mane_genbank(input_gz, output_csv):
    if not os.path.exists(input_gz):
        sys.exit(f"ERROR: Input file not found: {input_gz}")

    fieldnames = [
        "gene_name", "gene_definition", "transcript_accession", 
        "transcript_version", "species", "transcript_length", 
        "CDS_start", "CDS_end", "transcript_sequence"
    ]

    total_records = 0
    retained_records = 0
    start_time = time.time()

    log_info(f"Extracting from: {input_gz}")
    with open(output_csv, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()

        with gzip.open(input_gz, "rt", encoding="utf-8") as in_fh:
            for record in SeqIO.parse(in_fh, "genbank"):
                total_records += 1
                if total_records % 2000 == 0:
                    log_info(f"Processed {total_records} records...")

                has_cds, cds_start, cds_end, gene_name = False, None, None, ""

                for feature in record.features:
                    if feature.type == "gene":
                        gene_name = feature.qualifiers.get("gene", [""])[0]
                    elif feature.type == "CDS":
                        has_cds = True
                        cds_start = int(feature.location.start) + 1
                        cds_end = int(feature.location.end)
                        if not gene_name:
                            gene_name = feature.qualifiers.get("gene", [""])[0]

                if not has_cds: continue

                writer.writerow({
                    "gene_name": gene_name,
                    "gene_definition": record.description,
                    "transcript_accession": record.name,
                    "transcript_version": record.id,
                    "species": record.annotations.get("organism", "Unknown"),
                    "transcript_length": len(record.seq),
                    "CDS_start": cds_start,
                    "CDS_end": cds_end,
                    "transcript_sequence": str(record.seq).upper()
                })
                retained_records += 1

    log_info(f"Done! Retained {retained_records}/{total_records} transcripts.")
    log_info(f"Time elapsed: {time.time() - start_time:.2f} s -> {output_csv}")

if __name__ == "__main__":
    process_mane_genbank(GBFF_FILE, OUT_CSV)