#!/bin/bash
# 02.run_tiger2023.sh — TIGER 2023 CasRx prediction on NM_000371.4 at 23 nt guide length
# conda env: casrx_tiger2023
# Note: TIGER CNN input shape is fixed at 208 (=23nt guide). Cannot run at other lengths.
#
# Usage: nohup bash 02.run_tiger2023.sh > 02.run_tiger2023.log 2>&1 &

set -e

PROJECT_DIR="/data/cai801/data/HKUcas"
SOFTWARE_DIR="${PROJECT_DIR}/software/tiger/hugging_face"
DATA_DIR="${PROJECT_DIR}/data"
RESULT_DIR="${PROJECT_DIR}/result/08.manuscript_03/partA/02.original_pred"
ACC="NM_000371.4"
LEN=23
TMP_DIR=$(mktemp -d)

mkdir -p "${RESULT_DIR}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log "02. TIGER 2023: ${ACC} (${LEN} nt guide)"
log "    Note: TIGER CNN input is fixed at 208 dims (=23nt guide). Only 23nt is supported."

# Copy only the target FASTA to a temp directory (TIGER scans all .fasta files in the directory)
cp "${DATA_DIR}/${ACC}.fasta" "${TMP_DIR}/"
log "    Copied ${ACC}.fasta to temp directory"

cd "${SOFTWARE_DIR}"
source /data/cai801/miniconda3/etc/profile.d/conda.sh
conda activate casrx_tiger2023

python tiger.py --mode all --fasta_path "${TMP_DIR}"

if [ -f "on_target.csv" ]; then
    mv on_target.csv "${RESULT_DIR}/TIGER2023_${LEN}_${ACC}.csv"
    log "    Saved: TIGER2023_${LEN}_${ACC}.csv"
else
    log "    ERROR: on_target.csv not generated"
    rm -rf "${TMP_DIR}"
    exit 1
fi

# Cleanup
rm -rf "${TMP_DIR}"
log "    Cleaned up temp directory"

log "02. TIGER 2023: Done."
