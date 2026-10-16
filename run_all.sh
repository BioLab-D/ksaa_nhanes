#!/usr/bin/env bash
# Reproduce every statistical result, table and figure of the manuscript from the files in data/.
# Run from the repository root:  bash run_all.sh            (all steps)
#                                bash run_all.sh manuscript_tables   (one step; earlier steps must have run)
# Requirements: R >= 4.3 with the packages in environment/R_packages.txt;
# a Python environment with CLARITE 2.3.6, Pillow and PyMuPDF for Figure 4 (set CLARITE_PY, default: python3).
# Each step writes its console output to results/logs/<step>.txt.
set -euo pipefail
cd "$(dirname "$0")"
CLARITE_PY=${CLARITE_PY:-python3}
LOG=results/logs; mkdir -p $LOG

stage_inputs() {   # copy the shipped inputs to the locations the scripts read from
  mkdir -p results/work results/predictions results/analysis/phewas results/analysis/tables results/figures
  cp data/prepared/* results/work/
  cp data/images/boxes_rot.csv data/images/boxes_lowconf.csv data/images/metal_scan.csv data/images/cdc_knee_seqn.txt results/work/
  cp data/predictions/* results/predictions/
  cp data/phewas/* results/analysis/phewas/
}

run() {   # run <script> [env assignments...]
  local s=$1; shift
  echo "$(date '+%F %T')  $s"
  env "$@" Rscript "code/analysis/$s.R" > "$LOG/$s.txt" 2>&1 || { echo "FAILED: $s (see $LOG/$s.txt)"; exit 1; }
}

STEPS=(
  survival_labels           # survival labels (NHANES III linked mortality; time from the examination, PERMTH_EXM)
  ksaa_ensemble             # five-backbone fold-wise weighted ensemble + training-fold bias correction -> KSAA
  build_analysis_data       # analysis data set: covariates, Phenotypic Age, outcomes
  phenoage_comparison fracture_function_vars quartile_characteristics   # Table S3 P values use the weighted design, which reads fxf_vars.csv
  weighted weighted_boot    # survey-weighted primary analyses (design on all adult examinees, cohort as a domain) -> results/analysis/weighted/
  phewas_audit_counts phewas_audit phewas phewas_table
  fig2_mortality joint_quartiles fig3_incremental fig1 fig4 figS1_biascorr
  kl_bmd_partial_cor        # partial correlations of KSAA with the K-L grade and proximal femur BMD
  manuscript_tables tableS1_performance manuscript_numbers
  verify_derived            # compare the regenerated results with data/derived/ (stops with an error on any difference)
)

step() {
  case $1 in
    stage)               stage_inputs ;;
    survival_labels)     run survival_labels ;;
    ksaa_ensemble)       run ksaa_ensemble; cp $LOG/ksaa_ensemble.txt results/ksaa_ensemble_log.txt ;;   # tableS1_performance.R cross-checks this log
    phewas_audit_counts) run phewas_clean; mv $LOG/phewas_clean.txt $LOG/phewas_audit_counts.txt ;;   # per-variable counts for the audit + cleaned phenotypes
    fig4)                run fig4_phewas                                                         # Figure 4: panel B and data in R, panel A and assembly with CLARITE
                         $CLARITE_PY code/analysis/fig4_phewas_clarite.py . > $LOG/fig4_clarite.txt 2>&1 || { echo "FAILED: Figure 4A (see $LOG/fig4_clarite.txt)"; exit 1; } ;;
    *)                   run $1 ;;
  esac
}

if [ $# -gt 0 ]; then for s in "$@"; do step $s; done; exit 0; fi
step stage
for s in "${STEPS[@]}"; do step $s; done
echo "$(date '+%F %T')  done. Tables: results/analysis/tables/*.html; figures: results/figures/; logs: $LOG/"
