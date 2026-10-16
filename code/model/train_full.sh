#!/bin/bash
# Same recipe as the five folds: each of the five backbones is trained once on all 4,875 participants and saved (10% internal validation for early stopping only, no held-out fold), for new radiographs.
# Manuscript results come from the five-fold OOF predictions, not from these models. Deployment: single-model prediction = mean expected age of the two knees' local heads; ensemble weights and bias correction use the values calibrated on the five-fold OOF predictions,
#   never re-calibrated on the final models' in-sample predictions (that would underestimate the bias). Each configuration is retried once if no checkpoint is written. Usage: bash train_full.sh
export TORCH_HOME=/home/data/torch PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True N_VIEWS=2 DROP=0 NO_ALBUMENTATIONS_UPDATE=1
PY=${PY:-python}; cd /tmp/mammo_age_official; D=/home/data/knee_yolo/v2; W=$D/weights; mkdir -p $W
F='grep --line-buffered -vE check_version|fetch_version|validated_self|UserWarning|data.=.fetch'
C="--bs 16 --seed 20260831 --nblock 4 --img 448 --full --save_dir $W"
BASE="--aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15"
Q=$D/p5_queue.log
run(){ t=$1; shift; L=$D/$t.log
  if [ -s $W/${t}_full.pt ]; then echo "$(date '+%F %T')  skip $t (checkpoint exists)" >> $Q; return; fi
  for try in 1 2; do
    echo "=== $t attempt $try $(date '+%F %T')  $* ===" >> $L
    $PY knee_v2.py --tag $t $C $BASE "$@" 2>&1 | $F >> $L
    if [ -s $W/${t}_full.pt ] && grep -q "full done" $L; then
      echo "$(date '+%F %T')  $t  $(tail -n 1 $D/${t}_folds.csv | tr -d '\r')  checkpoint $(du -h $W/${t}_full.pt | cut -f1)" >> $Q; return
    fi
    echo "$(date '+%F %T')  $t attempt $try wrote no checkpoint (see $t.log)" >> $Q; sleep 90
  done
}
echo "=== v2 final-model checkpoint queue start $(date '+%F %T') ===" > $Q
run m5_full_resnet18       --arch resnet18
run m5_full_efficientnetb0 --arch efficientnet_b0
run m5_full_resnet50       --arch resnet50
run m5_full_densenet121    --arch densenet121
run m5_full_convnexttiny   --arch convnext_tiny
echo "=== done $(date '+%F %T') ===" >> $Q
