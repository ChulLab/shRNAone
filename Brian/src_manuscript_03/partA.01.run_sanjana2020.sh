#!/bin/bash
# 01.run_sanjana2020.sh — Sanjana 2020 CasRx prediction on NM_000371.4 at multiple guide lengths
# conda env: casrx_sanjana2020
# Guide lengths: 23, 30 nt (18nt unsupported — internal features require >=23nt)
#
# Strategy: 
#   1. ALWAYS force-reset R script to default (23) at the start of each iteration
#   2. Replace hardcoded "23" at exact line numbers via Python
#   3. Run, then restore to default (23)
#   4. Trap ensures restoration even on crash
#
# Lines to modify:
#   Line 79:  guideLength = 23
#   Line 499: Y=23)
#   Line 524: MA[23,]
#   Line 592: E=23)
#
# Usage: nohup bash 01.run_sanjana2020.sh > 01.run_sanjana2020.log 2>&1 &

set -e

PROJECT_DIR="/data/cai801/data/HKUcas"
RSRC="${PROJECT_DIR}/software/cas13/Cas13designGuidePredictor/scripts/RfxCas13d_GuideScoring.R"
SOFTWARE_DIR="${PROJECT_DIR}/software/cas13/Cas13designGuidePredictor"
DATA_DIR="${PROJECT_DIR}/data"
RESULT_DIR="${PROJECT_DIR}/result/08.manuscript_03/partA/02.original_pred"
FASTA="${DATA_DIR}/NM_000371.4.fasta"
ACC="NM_000371.4"
LENGTHS=(23 30)
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

export PATH=/data/cai801/miniconda3/envs/casrx_sanjana2020/bin:$PATH
mkdir -p "${RESULT_DIR}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Force-reset the R script to its original default values (all 23)
# This is the SAFETY NET: guarantees clean state regardless of previous crashes
reset_to_default() {
    local file="$1"
    python3 - "$file" << 'PYEOF'
import sys

file = sys.argv[1]

replacements = {
    79:  "guideLength = 23 # fixed length",
    499: "GetUnpairedProb <- function( x = Guide.df , MA = log10(UnpairedProbabilities.tranformed), W=50, X=40 , Y=23){",
    524: "  VEC = MA[23,]  # is similar to Y",
    592: "GetLetterProbs = function( x = Guide.df, S=1, E=23){",
}

with open(file, 'r') as f:
    lines = f.readlines()

changed = False
for lineno, expected in replacements.items():
    idx = lineno - 1
    stripped = lines[idx].rstrip('\n')
    if stripped != expected:
        lines[idx] = expected + '\n'
        changed = True

if changed:
    with open(file, 'w', newline='\n') as f:
        f.writelines(lines)
    print("RESET", end='')
else:
    print("CLEAN", end='')
PYEOF
}

# Replace specific lines in the R script using Python (line-number exact)
# Args: $1=file $2=target_length
hack_lines() {
    local file="$1"
    local target="$2"
    python3 - "$file" "$target" << 'PYEOF'
import sys

file = sys.argv[1]
target = sys.argv[2]

replacements = {
    79:  ("guideLength = 23", f"guideLength = {target}"),
    499: ("Y=23)", f"Y={target})"),
    524: ("MA[23,]", f"MA[{target},]"),
    592: ("E=23)", f"E={target})"),
}

with open(file, 'r') as f:
    lines = f.readlines()

modified = False
for lineno, (old, new) in replacements.items():
    idx = lineno - 1
    if old in lines[idx]:
        lines[idx] = lines[idx].replace(old, new, 1)
        modified = True

if modified:
    with open(file, 'w', newline='\n') as f:
        f.writelines(lines)
PYEOF
}

log "01. Sanjana 2020: Starting prediction on ${ACC}"
log "    Timestamp: ${TIMESTAMP}"

for LEN in "${LENGTHS[@]}"; do
    log "    =============================================="
    log "    Processing guide length: ${LEN} nt"

    # STEP 1: FORCE RESET to default (23) — the critical safety net
    reset_status=$(reset_to_default "${RSRC}")
    if [ "${reset_status}" = "RESET" ]; then
        log "    [SAFETY] R script was dirty — forced reset to default (23)"
    else
        log "    [OK] R script already clean (default 23)"
    fi

    # STEP 2: Hack to target length
    hack_lines "${RSRC}" "${LEN}"
    if [ "${LEN}" -eq 23 ]; then
        log "    No modification needed (target length equals default 23)"
    else
        log "    Modified lines 79, 499, 524, 592 in RfxCas13d_GuideScoring.R (23 -> ${LEN})"
    fi

    # Copy FASTA to software dir (required by R script)
    cp "${FASTA}" "${SOFTWARE_DIR}/${ACC}.fasta"

    # Set trap: if ANYTHING fails, force-reset to default before exiting
    trap 'log "    ERROR: force-resetting R script to default (23)..."; reset_to_default "${RSRC}"; exit 1' ERR

    # Run R script
    cd "${SOFTWARE_DIR}"
    Rscript ./scripts/RfxCas13d_GuideScoring.R \
        "${SOFTWARE_DIR}/${ACC}.fasta" \
        ./data/Cas13designGuidePredictorInput.csv \
        true

    # Success — clear trap
    trap - ERR

    # Collect and rename output
    for csv in ./*.csv; do
        base=$(basename "$csv")
        target="${RESULT_DIR}/Sanjana2020_${LEN}_${ACC}.csv"
        mv "$csv" "$target"
        log "    Saved: $(basename $target)"
    done

    # Cleanup
    rm -f "${SOFTWARE_DIR}/${ACC}.fasta" ./*.pdf 2>/dev/null

    # STEP 3: Force-reset back to default (23) after successful run
    reset_to_default "${RSRC}" > /dev/null
    log "    Force-restored R script to default (23) after completion"
    log "    Completed guide length ${LEN} nt successfully"
done

log "01. Sanjana 2020: All lengths processed. R script confirmed at default (23). Done."
