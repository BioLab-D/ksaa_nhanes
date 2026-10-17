#!/bin/bash
export TORCH_HOME=/home/data/torch PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True N_VIEWS=2 DROP=0 NO_ALBUMENTATIONS_UPDATE=1
PY=${PY:-python}; cd /tmp/mammo_age_official; D=${D:-/home/data/knee_yolo/v2}; mkdir -p $D
F='grep --line-buffered -vE check_version|fetch_version|validated_self|UserWarning|data.=.fetch'
C="--bs 16 --seed 20260831 --folds 0,1,2,3,4 --nblock 4"
BASE=${BASE:-"--aug knee --norm pct --sched cosine --wd 1e-4 --warmup 5 --epochs 80 --min_ep 60 --patience 15"}
Q=$D/cv_queue.log
run(){ t=$1; shift; L=$D/$t.log
  if [ -s $D/${t}_folds.csv ] && [ "$(($(wc -l < $D/${t}_folds.csv)-1))" -ge 5 ]; then
    echo "$(date '+%F %T')  skip $t (all five folds done)" >> $Q; return
  fi
  echo "=== $t $(date '+%F %T')  $* ===" > $L
  $PY knee_v2.py --tag $t $C $BASE "$@" 2>&1 | $F >> $L
  echo "=== end $(date '+%F %T') ===" >> $L
  if [ -s $D/${t}_folds.csv ]; then
    v=$(tail -n +2 $D/${t}_folds.csv | cut -d, -f2)
    ms=$(echo "$v" | awk '{s+=$1;q+=$1*$1;n++}END{printf "%.4f +/- %.4f", s/n, sqrt((q-s*s/n)/(n-1))}')
    echo "$(date '+%F %T')  $t  $(echo $v | tr '\n' ' ') | mean $ms" >> $Q
  else echo "$(date '+%F %T')  $t no output (see $t.log)" >> $Q; fi
}
echo "=== Mammo-AGE-style comparison (N=4, five folds) start $(date '+%F %T')  BASE=$BASE ===" > $Q

# 1. five backbones (manuscript)
echo "$(date '+%F %T')  --- section 1: five backbones ---" >> $Q
run cv_resnet18       --arch resnet18       --img 448
run cv_resnet50       --arch resnet50       --img 448
run cv_efficientnetb0 --arch efficientnet_b0 --img 448
run cv_convnexttiny   --arch convnext_tiny  --img 448
run cv_densenet121    --arch densenet121    --img 448

# 2. module ablations, ResNet-18 (not reported)
echo "$(date '+%F %T')  --- section 2: module ablations ---" >> $Q
run cv_ab_wo_citf --arch resnet18 --img 448 --nblock 0     # w/o CI-Transformer
run cv_ab_wo_poe  --arch resnet18 --img 448 --no_poe       # w/o POE
run cv_ab_wo_mvl  --arch resnet18 --img 448 --no_mvl       # w/o MVL
run cv_ab_wo_mt   --arch resnet18 --img 448 --no_density   # w/o density branch (no auxiliary labels)

# 3. input size, ResNet-18 (not reported; 448 px is in section 1)
echo "$(date '+%F %T')  --- section 3: input size ---" >> $Q
run cv_sz112 --arch resnet18 --img 112
run cv_sz224 --arch resnet18 --img 224
run cv_sz320 --arch resnet18 --img 320
run cv_sz512 --arch resnet18 --img 512
echo "=== Mammo-AGE-style comparison done $(date '+%F %T') ===" >> $Q
