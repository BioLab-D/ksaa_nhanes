#!/bin/bash
# Five-fold cross-fitting of the multi-view Mammo-AGE knee age model (N = 4 fusion blocks, 448 ROI); writes per-fold internal validation MAE and OOF predictions.
# Recipe adapted to knee films: knee augmentation + per-image P0.5–P99.5 float normalization + AdamW with cosine annealing (wd 1e-4).
# Section 1 is the five backbones used in the manuscript; sections 2 and 3 are ablations not reported in the manuscript. Usage: bash train_5fold.sh
export TORCH_HOME=/home/data/torch PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True N_VIEWS=2 DROP=0 NO_ALBUMENTATIONS_UPDATE=1
PY=${PY:-python}; cd /tmp/mammo_age_official; D=${D:-/home/data/knee_yolo/v2}; mkdir -p $D
F='grep --line-buffered -vE check_version|fetch_version|validated_self|UserWarning|data.=.fetch'
C="--bs 16 --seed 20260831 --folds 0,1,2,3,4 --nblock 4"
BASE=${BASE:-"--aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15"}
Q=$D/p4_queue.log
run(){ t=$1; shift; L=$D/$t.log
  [ -s $D/${t}_folds.csv ] && [ "$(($(wc -l < $D/${t}_folds.csv)-1))" -ge 5 ] && { echo "$(date '+%F %T')  skip $t (all five folds done)" >> $Q; return; }
  echo "=== $t $(date '+%F %T')  $* ===" > $L
  $PY knee_v2.py --tag $t $C $BASE "$@" 2>&1 | $F >> $L
  echo "=== end $(date '+%F %T') ===" >> $L
  if [ -s $D/${t}_folds.csv ]; then
    v=$(tail -n +2 $D/${t}_folds.csv | cut -d, -f2)
    echo "$(date '+%F %T')  $t  $(echo $v | tr '\n' ' ') | mean $(echo "$v" | awk '{s+=$1;q+=$1*$1;n++}END{printf "%.4f ± %.4f", s/n, sqrt((q-s*s/n)/(n-1))}')" >> $Q
  else echo "$(date '+%F %T')  $t no output (see $t.log)" >> $Q; fi
}
echo "=== Mammo-AGE-style comparison (N=4, five folds) start $(date '+%F %T')  BASE=$BASE ===" > $Q

# ── 1 five backbones (used in the manuscript)
echo "$(date '+%F %T')  --- section 1: five backbones ---" >> $Q
run m4_resnet18       --arch resnet18       --img 448
run m4_resnet50       --arch resnet50       --img 448
run m4_efficientnetb0 --arch efficientnet_b0 --img 448
run m4_convnexttiny   --arch convnext_tiny  --img 448
run m4_densenet121    --arch densenet121    --img 448

# ── 2 module ablations @ ResNet-18 (not reported in the manuscript)
echo "$(date '+%F %T')  --- section 2: module ablations ---" >> $Q
run m4_ab_wo_citf --arch resnet18 --img 448 --nblock 0     # w/o CI-Transformer
run m4_ab_wo_poe  --arch resnet18 --img 448 --no_poe       # w/o POE
run m4_ab_wo_mvl  --arch resnet18 --img 448 --no_mvl       # w/o MVL
run m4_ab_wo_mt   --arch resnet18 --img 448 --no_density   # w/o Multi-task (no auxiliary knee labels)

# ── 3 input size @ ResNet-18 (not reported in the manuscript; 448 is in section 1)
echo "$(date '+%F %T')  --- section 3: input size ---" >> $Q
run m4_sz112 --arch resnet18 --img 112
run m4_sz224 --arch resnet18 --img 224
run m4_sz320 --arch resnet18 --img 320
run m4_sz512 --arch resnet18 --img 512
echo "=== Mammo-AGE-style comparison done $(date '+%F %T') ===" >> $Q
