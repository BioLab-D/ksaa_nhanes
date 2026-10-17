#!/usr/bin/env bash
# Reproduce every statistical result, table and figure of the manuscript from the files in data/.
# Run from the repository root:  bash run_all.sh          all scripts in code/analysis, in numbered order
#                                bash run_all.sh 13 14    only the scripts with these numbers (earlier ones must have run)
# Requirements: R >= 4.3 with the packages in environment/R_packages.txt; for Figure 4 a Python environment with
# CLARITE 2.3.6, Pillow and PyMuPDF (set CLARITE_PY, default: python3).
# Console output of each script goes to results/logs/<script>.txt.
set -euo pipefail
cd "$(dirname "$0")"
CLARITE_PY=${CLARITE_PY:-python3}
LOG=results/logs; mkdir -p $LOG

run() {
  local f=$1 s; s=$(basename "${f%.*}")
  echo "$(date '+%F %T')  $s"
  case $f in
    *.R)  Rscript "$f" > "$LOG/$s.txt" 2>&1 ;;
    *.py) $CLARITE_PY "$f" . > "$LOG/$s.txt" 2>&1 ;;
  esac || { echo "FAILED: $s (see $LOG/$s.txt)"; exit 1; }
}

if [ $# -gt 0 ]; then
  for n in "$@"; do run code/analysis/${n}_*.*; done
  exit 0
fi

# copy the shipped inputs to where the scripts read them
mkdir -p results/work results/predictions results/analysis/phewas results/analysis/tables results/figures
cp data/prepared/* data/images/*.csv data/images/cdc_knee_seqn.txt results/work/
cp data/predictions/* results/predictions/
cp data/phewas/* results/analysis/phewas/

for f in code/analysis/[0-9][0-9]_*.*; do run "$f"; done
run code/analysis/verify_derived.R   # compare the regenerated results with data/derived/; stops on any difference
echo "$(date '+%F %T')  done. Tables: results/analysis/tables/; figures: results/figures/; logs: $LOG/"
